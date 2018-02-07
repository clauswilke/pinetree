#include <cmath>
#include <fstream>
#include <iostream>

#include "choices.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include "tracker.hpp"

const static double AVAGADRO = double(6.0221409e+23);

SpeciesReaction::SpeciesReaction(double rate_constant, double volume,
                                 const std::vector<std::string> &reactants,
                                 const std::vector<std::string> &products)
    : rate_constant_(rate_constant),
      reactants_(reactants),
      products_(products) {
  // Error checking
  if (reactants_.size() > 2) {
    throw std::runtime_error(
        "Simulation does not support reactions with "
        "more than two reactant species.");
  }
  if (volume <= 0) {
    throw std::runtime_error("Reaction volume cannot be zero.");
  }
  // Rate constant has to be transformed from macroscopic to mesoscopic
  // TODO: make more sophisticated and support different stoichiometries
  if (reactants_.size() == 2) {
    rate_constant_ = rate_constant_ / (AVAGADRO * volume);
  }
}

double SpeciesReaction::CalculatePropensity() {
  double propensity = rate_constant_;
  for (const auto &reactant : reactants_) {
    propensity *= SpeciesTracker::Instance().species(reactant);
  }
  return propensity;
}

void SpeciesReaction::Execute() {
  for (const auto &reactant : reactants_) {
    SpeciesTracker::Instance().Increment(reactant, -1);
  }
  for (const auto &product : products_) {
    SpeciesTracker::Instance().Increment(product, 1);
  }
}

Bind::Bind(double rate_constant, double volume,
           const std::string &promoter_name, const Polymerase &pol_template)
    : rate_constant_(rate_constant),
      promoter_name_(promoter_name),
      pol_name_(pol_template.name()),
      pol_template_(pol_template) {
  // Check volume
  if (volume <= 0) {
    throw std::runtime_error("Reaction volume cannot be zero.");
  }
  auto &tracker = SpeciesTracker::Instance();
  rate_constant_ = rate_constant_ / (AVAGADRO * volume);
}

double Bind::CalculatePropensity() {
  auto &tracker = SpeciesTracker::Instance();
  return rate_constant_ * tracker.species(pol_name_) *
         tracker.species(promoter_name_);
}

void Bind::Execute() {
  auto &tracker = SpeciesTracker::Instance();
  auto weights = std::vector<double>();
  for (const auto &polymer : tracker.FindPolymers(promoter_name_)) {
    weights.push_back(double(polymer->uncovered(promoter_name_)));
  }
  Polymer::Ptr polymer =
      Random::WeightedChoice(tracker.FindPolymers(promoter_name_), weights);
  auto new_pol = std::make_shared<Polymerase>(pol_template_);
  polymer->Bind(new_pol, promoter_name_);
  tracker.propensity_signal_.Emit(polymer->index());
  tracker.Increment(promoter_name_, -1);
  tracker.Increment(pol_name_, -1);
}

Simulation::Simulation(double cell_volume)
    : time_(0), alpha_sum_(0), cell_volume_(cell_volume) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();
}

void Simulation::seed(int seed) { Random::seed(seed); }

void Simulation::Run(int stop_time, int time_step,
                     const std::string &output_prefix) {
  auto &tracker = SpeciesTracker::Instance();
  Initialize();
  // Set up file output streams
  std::ofstream countfile(output_prefix + "_counts.tsv", std::ios::trunc);
  countfile << "time\tspecies\tcount\ttranscript\tribo_density\n";
  int out_time = 0;
  std::map<std::string, std::vector<double>> output;
  while (time_ < stop_time) {
    if ((out_time - time_) < 0.001) {
      for (auto elem : tracker.species()) {
        output[elem.first].push_back(elem.second);
        output[elem.first].push_back(0);
        output[elem.first].push_back(0);
      }
      for (auto transcript : tracker.transcripts()) {
        if (output[transcript.first].size() == 3) {
          output[transcript.first][1] = transcript.second;
        } else {
          output[transcript.first].push_back(0);
          output[transcript.first].push_back(transcript.second);
          output[transcript.first].push_back(0);
        }
        if (tracker.ribo_per_transcript().find(transcript.first) !=
            tracker.ribo_per_transcript().end()) {
          output[transcript.first][2] =
              double(tracker.ribo_per_transcript(transcript.first)) /
              double(output[transcript.first][1]);
        }
      }
      for (auto row : output) {
        countfile << (std::to_string(time_) + "\t" + row.first + "\t" +
                      std::to_string(row.second[0]) + "\t" +
                      std::to_string(row.second[1]) + "\t" +
                      std::to_string(row.second[2]) + "\n");
      }
      output.clear();
      countfile.flush();
      out_time += time_step;
    }
    Execute();
  }
  countfile.close();
}

void Simulation::RegisterReaction(Reaction::Ptr reaction) {
  auto it = std::find(reactions_.begin(), reactions_.end(), reaction);
  if (it == reactions_.end()) {
    reaction->index(reactions_.size());
    double new_prop = reaction->CalculatePropensity();
    alpha_list_.push_back(new_prop);
    alpha_sum_ += new_prop;
    reactions_.push_back(reaction);
  }
}

void Simulation::AddReaction(double rate_constant,
                             const std::vector<std::string> &reactants,
                             const std::vector<std::string> &products) {
  auto rxn = std::make_shared<SpeciesReaction>(rate_constant, cell_volume_,
                                               reactants, products);
  auto &tracker = SpeciesTracker::Instance();
  for (const auto &reactant : reactants) {
    tracker.Add(reactant, rxn);
    tracker.Increment(reactant, 0);
  }
  for (const auto &product : products) {
    tracker.Add(product, rxn);
    tracker.Increment(product, 0);
  }
  RegisterReaction(rxn);
}

void Simulation::AddSpecies(const std::string &name, int copy_number) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Increment(name, copy_number);
}

void Simulation::AddPolymerase(const std::string &name, int footprint,
                               double mean_speed, int copy_number) {
  auto pol = Polymerase(name, footprint, mean_speed);
  polymerases_.push_back(pol);
  AddSpecies(name, copy_number);
}

void Simulation::RegisterPolymer(Polymer::Ptr polymer) {
  polymer->InitElements();
  for (auto &elem : polymer->elements()) {
    if (elem->type() == "promoter") {
      elem->uncover_signal_.ConnectMember(shared_from_this(),
                                          &Simulation::FreePromoter);
      elem->cover_signal_.ConnectMember(shared_from_this(),
                                        &Simulation::BlockPromoter);
      SpeciesTracker::Instance().Add(elem->name(), polymer);
    }
  }
  // Encapsulate polymer in Bridge reaction and add to reaction list
  auto bridge = std::make_shared<Bridge>(polymer);
  RegisterReaction(bridge);
  polymer->index(bridge->index());
}

void Simulation::RegisterGenome(Genome::Ptr genome) {
  RegisterPolymer(genome);
  genome->termination_signal_.ConnectMember(
      shared_from_this(), &Simulation::TerminateTranscription);
  genome->transcript_signal_.ConnectMember(shared_from_this(),
                                           &Simulation::RegisterTranscript);
  genomes_.push_back(genome);
}

void Simulation::RegisterTranscript(Transcript::Ptr transcript) {
  RegisterPolymer(transcript);
  transcript->termination_signal_.ConnectMember(
      shared_from_this(), &Simulation::TerminateTranslation);
}

void Simulation::Initialize() {
  auto &tracker = SpeciesTracker::Instance();
  if (time_ == 0) {
    tracker.propensity_signal_.ConnectMember(shared_from_this(),
                                             &Simulation::UpdatePropensity);
    InitBindReactions();
    InitPropensity();
  }
  if (genomes_.size() == 0) {
    std::cerr << "Warning: There are no Genome objects registered with "
                 "Simulation. Did you forget to register a Genome?"
              << std::endl;
  }
}

void Simulation::InitPropensity() {
  for (int i = 0; i < reactions_.size(); i++) {
    UpdatePropensity(i);
  }
  // Make sure prop sum is starting from 0
  alpha_sum_ = 0;
  for (const auto &alpha : alpha_list_) {
    alpha_sum_ += alpha;
  }
}

void Simulation::InitBindReactions() {
  for (Genome::Ptr genome : genomes_) {
    for (auto promoter_name : genome->bindings()) {
      for (auto pol : polymerases_) {
        if (promoter_name.second.count(pol.name()) != 0) {
          double rate_constant = promoter_name.second[pol.name()];
          Polymerase pol_template = Polymerase(pol);
          auto reaction = std::make_shared<Bind>(
              rate_constant, cell_volume_, promoter_name.first, pol_template);
          auto &tracker = SpeciesTracker::Instance();
          tracker.Add(promoter_name.first, reaction);
          tracker.Add(pol.name(), reaction);
          RegisterReaction(reaction);
        }
      }
    }
  }
}

void Simulation::UpdatePropensity(int index) {
  double new_prop = reactions_[index]->CalculatePropensity();
  double diff = new_prop - alpha_list_[index];
  alpha_sum_ += diff;
  alpha_list_[index] = new_prop;
}

void Simulation::Execute() {
  // Generate random number
  if (alpha_sum_ <= 0) {
    throw std::runtime_error("Propensity of system is 0.");
  }
  // InitPropensity();
  double random_num = Random::random();
  // Calculate tau, i.e. time until next reaction
  double tau = (1.0 / alpha_sum_) * std::log(1.0 / random_num);
  if (!std::isnormal(tau)) {
    throw std::runtime_error("Underflow error.");
  }
  time_ += tau;
  // Randomly select next reaction to execute, weighted by propensities
  auto next_reaction = Random::WeightedChoiceIndex(reactions_, alpha_list_);
  reactions_[next_reaction]->Execute();
  UpdatePropensity(next_reaction);
  iteration_++;
}

void Simulation::FreePromoter(const std::string &species_name) {
  SpeciesTracker::Instance().Increment(species_name, 1);
}
void Simulation::BlockPromoter(const std::string &species_name) {
  SpeciesTracker::Instance().Increment(species_name, -1);
}

void Simulation::TerminateTranscription(int polymer_index,
                                        const std::string &pol_name,
                                        const std::string &gene_name) {
  SpeciesTracker::Instance().Increment(pol_name, 1);
  UpdatePropensity(polymer_index);
  CountTermination("transcript");
}

void Simulation::TerminateTranslation(int polymer_index,
                                      const std::string &pol_name,
                                      const std::string &gene_name) {
  SpeciesTracker::Instance().Increment(pol_name, 1);
  SpeciesTracker::Instance().Increment(gene_name, 1);
  SpeciesTracker::Instance().IncrementRibo(gene_name, -1);
  UpdatePropensity(polymer_index);
  CountTermination(gene_name);
}

void Simulation::CountTermination(const std::string &name) {
  auto new_name = name + "_total";
  if (terminations_.count(name) == 0) {
    terminations_[new_name] = 1;
  } else {
    terminations_[new_name]++;
  }
}
