#include "reaction.hpp"
#include "choices.hpp"
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
