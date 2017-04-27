#include "simulation.hpp"

const static double AVAGADRO = double(6.0221409e+23);
const static double CELL_VOLUME = double(8e-15);

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

SpeciesReaction::SpeciesReaction(SpeciesTracker::Ptr tracker,
                                 double rate_constant,
                                 const std::vector<std::string> &reactants,
                                 const std::vector<std::string> &products)
    : tracker_(tracker), rate_constant_(rate_constant), reactants_(reactants),
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
  // Add this reaction to SpeciesTracker
  for (const auto &reactant : reactants_) {
    tracker_->Add(reactant, Ptr(this));
    // Initialize reactant counts (this usually gets set to some non-zero value
    // later)
    tracker_->Increment(reactant, 0);
  }
  // Do the same for products
  for (const auto &product : products_) {
    tracker_->Add(product, Ptr(this));
    tracker_->Increment(product, 0);
  }
}

double SpeciesReaction::CalculatePropensity() {
  double propensity = rate_constant_;
  for (const auto &reactant : reactants_) {
    propensity *= tracker_->species(reactant);
  }
  return propensity;
}

void SpeciesReaction::Execute() {
  for (const auto &reactant : reactants_) {
    tracker_->Increment(reactant, -1);
  }
  for (const auto &product : products_) {
    tracker_->Increment(product, 1);
  }
}