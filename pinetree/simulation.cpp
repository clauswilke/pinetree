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
    : rate_constant_(rate_constant), reactants_(reactants),
      products_(products) {
  // Error checking
  if (reactants_.size() > 2) {
    throw std::runtime_error("Simulation does not support reactions with "
                             "more than two reactant species.");
  }
  if (volume <= 0) {
    throw std::runtime_error("Reaction volume cannot be zero.");
  }
  // Rate constant has to be transformed from macroscopic to mesoscopic
  // TODO: make more sophisticated and support different stoichiometries
  auto &tracker = SpeciesTracker::Instance();
  if (reactants_.size() == 2) {
    rate_constant_ = rate_constant_ / (AVAGADRO * volume);
  }
}

void SpeciesReaction::InitCounts() {
  auto &tracker = SpeciesTracker::Instance();
  for (const auto &reactant : reactants_) {
    tracker.Add(reactant, shared_from_this());
    tracker.Increment(reactant, 0);
  }
  for (const auto &product : products_) {
    tracker.Add(product, shared_from_this());
    tracker.Increment(product, 0);
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
    : rate_constant_(rate_constant), promoter_name_(promoter_name),
      pol_name_(pol_template.name()), pol_template_(pol_template) {
  // Check volume
  if (volume <= 0) {
    throw std::runtime_error("Reaction volume cannot be zero.");
  }
  auto &tracker = SpeciesTracker::Instance();
  rate_constant_ = rate_constant_ / (AVAGADRO * volume);
}

void Bind::InitCounts() {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Add(promoter_name_, shared_from_this());
  tracker.Add(pol_name_, shared_from_this());
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

Simulation::Simulation()
    : time_(0), stop_time_(0), time_step_(0), alpha_sum_(0) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();
}

void Simulation::Run(const std::string &output_name) {
  auto &tracker = SpeciesTracker::Instance();
  if (time_ == 0) {
    tracker.propensity_signal_.ConnectMember(shared_from_this(),
                                             &Simulation::UpdatePropensity);
    InitPropensity();
  }
  // Set up file output streams
  std::ofstream countfile(output_name + "_counts.tsv", std::ios::trunc);
  std::ofstream ribofile(output_name + "_density.tsv", std::ios::trunc);
  int out_time = 0;
  while (time_ < stop_time_) {
    Execute();
    if ((out_time - time_) < 0.001) {
      for (auto elem : tracker.species()) {
        countfile << (std::to_string(time_) + "\t" + elem.first + "\t" +
                      std::to_string(elem.second))
                  << std::endl;
      }
      countfile.flush();
      for (auto elem : tracker.ribo_per_transcript()) {
        double density =
            double(elem.second) / double(tracker.transcripts(elem.first));
        ribofile << (std::to_string(time_) + "\t" + elem.first + "\t" +
                     std::to_string(elem.second) + "\t" +
                     std::to_string(tracker.transcripts(elem.first)) + "\t" +
                     std::to_string(density))
                 << std::endl;
      }
      ribofile.flush();
      out_time += time_step_;
    }
  }
  countfile.close();
  ribofile.close();
}

void Simulation::RegisterReaction(Reaction::Ptr reaction) {
  reaction->InitCounts();
  auto it = std::find(reactions_.begin(), reactions_.end(), reaction);
  if (it == reactions_.end()) {
    reaction->index(reactions_.size());
    double new_prop = reaction->CalculatePropensity();
    alpha_list_.push_back(new_prop);
    alpha_sum_ += new_prop;
    reactions_.push_back(reaction);
  }
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
}

void Simulation::RegisterTranscript(Transcript::Ptr transcript) {
  RegisterPolymer(transcript);
  transcript->termination_signal_.ConnectMember(
      shared_from_this(), &Simulation::TerminateTranslation);
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

void Simulation::UpdatePropensity(int index) {
  double new_prop = reactions_[index]->CalculatePropensity();
  double diff = new_prop - alpha_list_[index];
  alpha_sum_ += diff;
  alpha_list_[index] = new_prop;
}

void Simulation::Execute() {
  // Generate random number
  if (alpha_sum_ == 0) {
    throw std::runtime_error("Propensity of system is 0.");
  }
  // InitPropensity();
  double random_num = Random::random();
  // Calculate tau, i.e. time until next reaction
  double tau = (1.0 / alpha_sum_) * std::log(1.0 / random_num);
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