#include <cmath>

#include "choices.hpp"
#include "polymer.hpp"
#include "simulation.hpp"

const static double AVAGADRO = double(6.0221409e+23);
const static double CELL_VOLUME = double(8e-15);

SpeciesReaction::SpeciesReaction(double rate_constant,
                                 const std::vector<std::string> &reactants,
                                 const std::vector<std::string> &products)
    : rate_constant_(rate_constant), reactants_(reactants),
      products_(products) {
  // Error checking
  if (reactants_.size() > 2) {
    throw std::runtime_error("Simulation does not support reactions with "
                             "more than two reactant species.");
  }
  // Rate constant has to be transformed from macroscopic to mesoscopic
  if (reactants_.size() == 2) {
    rate_constant_ = rate_constant_ / (AVAGADRO * CELL_VOLUME);
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

Bind::Bind(double rate_constant, const std::string &promoter_name,
           const Polymerase &pol_template)
    : rate_constant_(rate_constant), promoter_name_(promoter_name),
      pol_name_(pol_template.name()), pol_template_(pol_template) {
  rate_constant_ = rate_constant_ / (AVAGADRO * CELL_VOLUME);
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
  tracker.propensity_signal_.ConnectMember(this, &Simulation::UpdatePropensity);
}

void Simulation::Run() {
  // coming soon
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

void Simulation::RegisterPolymer(Polymer::Ptr polymer) {
  for (auto &elem : polymer->elements()) {
    elem->uncover_signal_.ConnectMember(this, &Simulation::FreePromoter);
    elem->cover_signal_.ConnectMember(this, &Simulation::BlockPromoter);
    SpeciesTracker::Instance().Add(elem->name(), polymer);
  }
  // Encapsulate polymer in Bridge reaction and add to reaction list
  auto bridge = std::make_shared<Bridge>(polymer);
  RegisterReaction(bridge);
  polymer->index(bridge->index());
}

void Simulation::RegisterGenome(Genome::Ptr genome) {
  RegisterPolymer(genome);
  genome->termination_signal_.ConnectMember(
      this, &Simulation::TerminateTranscription);
  genome->transcript_signal_.ConnectMember(this,
                                           &Simulation::RegisterTranscript);
}

void Simulation::RegisterTranscript(Transcript::Ptr transcript) {
  RegisterPolymer(transcript);
  transcript->termination_signal_.ConnectMember(
      this, &Simulation::TerminateTranslation);
}

void Simulation::UpdatePropensity(int index) {
  double new_prop = reactions_[index]->CalculatePropensity();
  double diff = new_prop - alpha_list_[index];
  alpha_sum_ += diff;
  alpha_list_[index] = new_prop;
}

void Simulation::Execute() {
  // Generate random number
  double random_num = Random::random();
  // Calculate tau, i.e. time until next reaction
  double tau = (1.0 / alpha_sum_) * std::log(1.0 / random_num);
  time_ += tau;
  // Randomly select next reaction to execute, weighted by propensities
  auto next_reaction = Random::WeightedChoice(reactions_, alpha_list_);
  next_reaction->Execute();
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

SpeciesTracker &SpeciesTracker::Instance() {
  static SpeciesTracker instance;
  return instance;
}

void SpeciesTracker::Clear() {
  species_.clear();
  promoter_map_.clear();
  species_map_.clear();
}

void SpeciesTracker::Register(SpeciesReaction::Ptr reaction) {
  // Add this reaction to SpeciesTracker
  for (const auto &reactant : reaction->reactants()) {
    Add(reactant, reaction);
    // Initialize reactant counts (this usually gets set to some non-zero
    // value later)
    Increment(reactant, 0);
  }
  // Do the same for products
  for (const auto &product : reaction->products()) {
    Add(product, reaction);
    Increment(product, 0);
  }
}

void SpeciesTracker::Increment(const std::string &species_name,
                               int copy_number) {
  if (species_.count(species_name) != 0) {
    species_[species_name] += copy_number;
  } else {
    species_[species_name] = copy_number;
  }
  if (copy_number == 0) {
    return;
  }
  if (species_map_.count(species_name) != 0) {
    for (const auto &reaction : species_map_[species_name]) {
      propensity_signal_.Emit(reaction->index());
    }
  }
}

void SpeciesTracker::Add(const std::string &species_name,
                         Reaction::Ptr reaction) {
  if (species_map_.count(species_name) == 0) {
    species_map_[species_name] = Reaction::VecPtr{reaction};
  } else {
    // TODO: Maybe use a better data type here like a set?
    auto it = std::find(species_map_[species_name].begin(),
                        species_map_[species_name].end(), reaction);
    if (it != species_map_[species_name].end()) {
      species_map_[species_name].push_back(reaction);
    }
  }
}

void SpeciesTracker::Add(const std::string &promoter_name,
                         Polymer::Ptr polymer) {
  if (promoter_map_.count(promoter_name) == 0) {
    promoter_map_[promoter_name] = Polymer::VecPtr{polymer};
  } else {
    promoter_map_[promoter_name].push_back(polymer);
  }
}

const Reaction::VecPtr &
SpeciesTracker::FindReactions(const std::string &species_name) {
  return species_map_[species_name];
}

const Polymer::VecPtr &
SpeciesTracker::FindPolymers(const std::string &promoter_name) {
  return promoter_map_[promoter_name];
}

int SpeciesTracker::species(const std::string &reactant) {
  return species_[reactant];
}