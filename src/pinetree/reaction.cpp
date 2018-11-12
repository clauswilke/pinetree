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
  old_prop_ = 0;
  if (reactants_.size() > 2) {
    throw std::invalid_argument(
        "Pinetree does not support reactions with "
        "more than two reactant species.");
  }
  if (reactants_.size() == 0 && products_.size() == 0) {
    throw std::invalid_argument(
        "You must specify at least one product or reactant.");
  }
  if (rate_constant <= 0) {
    throw std::invalid_argument("Reaction rate constant cannot be zero.");
  }
  if (volume <= 0) {
    throw std::invalid_argument("Reaction volume cannot be zero.");
  }
  // Rate constant has to be transformed from macroscopic to mesoscopic
  // TODO: make more sophisticated and support different stoichiometries
  if (reactants_.size() == 2) {
    rate_constant_ = rate_constant_ / (AVAGADRO * volume);
  }
}

double SpeciesReaction::CalculatePropensity() {
  if (remove_ == true) {
    old_prop_ = 0;
  }
  double new_prop = rate_constant_;
  for (const auto &reactant : reactants_) {
    new_prop *= SpeciesTracker::Instance().species(reactant);
  }
  double prop_diff = new_prop - old_prop_;
  old_prop_ = new_prop;
  return prop_diff;
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
           const std::string &promoter_name)
    : rate_constant_(rate_constant), promoter_name_(promoter_name) {
  old_prop_ = 0;
  // Check volume
  if (volume <= 0) {
    throw std::runtime_error("Reaction volume cannot be zero.");
  }
}

Polymer::Ptr Bind::ChoosePolymer() {
  auto &tracker = SpeciesTracker::Instance();
  auto weights = std::vector<double>();
  for (const auto &polymer : tracker.FindPolymers(promoter_name_)) {
    weights.push_back(double(polymer->uncovered(promoter_name_)));
  }
  Polymer::Ptr polymer =
      Random::WeightedChoice(tracker.FindPolymers(promoter_name_), weights);
  return polymer;
}

BindPolymerase::BindPolymerase(double rate_constant, double volume,
                               const std::string &promoter_name,
                               const Polymerase &pol_template)
    : Bind(rate_constant, volume, promoter_name), pol_template_(pol_template) {
  rate_constant_ = rate_constant_ / (AVAGADRO * volume);
}

double BindPolymerase::CalculatePropensity() {
  auto &tracker = SpeciesTracker::Instance();
  double new_prop = rate_constant_ * tracker.species(pol_template_.name()) *
                    tracker.species(promoter_name_);
  double prop_diff = new_prop - old_prop_;
  old_prop_ = new_prop;
  return prop_diff;
}

void BindPolymerase::Execute() {
  auto polymer = ChoosePolymer();
  auto new_pol = std::make_shared<Polymerase>(pol_template_);
  polymer->Bind(new_pol, promoter_name_);
  SpeciesTracker::Instance().propensity_signal_.Emit(polymer->wrapper());
  // Polymer should handle decrementing promoter
  SpeciesTracker::Instance().Increment(new_pol->name(), -1);
}

BindRnase::BindRnase(double rate_constant, double volume,
                     const Rnase &rnase_template, const std::string &name)
    : Bind(rate_constant, volume, name), pol_template_(rnase_template) {}

void BindRnase::Execute() {
  auto polymer = ChoosePolymer();
  auto new_pol = std::make_shared<Rnase>(pol_template_);
  polymer->Bind(new_pol, promoter_name_);
  SpeciesTracker::Instance().propensity_signal_.Emit(polymer->wrapper());
}

double BindRnase::CalculatePropensity() {
  if (remove_ == true) {
    old_prop_ = 0;
  }
  auto &tracker = SpeciesTracker::Instance();
  double new_prop = rate_constant_ * tracker.species(promoter_name_);
  double prop_diff = new_prop - old_prop_;
  old_prop_ = new_prop;
  return prop_diff;
}

PolymerWrapper::PolymerWrapper(Polymer::Ptr polymer) : polymer_(polymer) {
  old_prop_ = 0;
  polymer_->Initialize();
}

PolymerWrapper::~PolymerWrapper() {
  polymer_->Unlink();
  // std::cout << "Destroying polymer wrapper." << std::endl;
}

void PolymerWrapper::index(int index) {
  index_ = index;
  polymer_->index(index);
}

void PolymerWrapper::Execute() {
  polymer_->Execute();
  if (polymer_->degrade() == true && polymer_->attached() == false) {
    remove_ = true;
    // std::cout << "Removing polymer wrapper...\n" << std::endl;
  }
}