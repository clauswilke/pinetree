#include "gillespie.hpp"
#include "choices.hpp"

void Gillespie::LinkReaction(Reaction::Ptr reaction) {
  auto it = std::find(reactions_.begin(), reactions_.end(), reaction);
  if (it == reactions_.end()) {
    reaction->index(reactions_.size());
    double new_prop = reaction->CalculatePropensity();
    alpha_list_.push_back(new_prop);
    alpha_sum_ += new_prop;
    reactions_.push_back(reaction);
  }
}

void Gillespie::DeleteReaction(int index) {
  // Update alpha sum
  alpha_sum_ -= alpha_list_[index];
  // Remove from alpha list
  alpha_list_.erase(alpha_list_.begin() + index);
  // Remove from reactions list
  reactions_.erase(reactions_.begin() + index);
}

void Gillespie::UpdatePropensity(int index) {
  if (index >= reactions_.size() || index >= alpha_list_.size() || index < 0) {
    throw std::range_error("Gillespie: Reaction index out of range.");
  }
  double new_prop = reactions_[index]->CalculatePropensity();
  double diff = new_prop - alpha_list_[index];
  alpha_sum_ += diff;
  alpha_list_[index] = new_prop;
}

void Gillespie::Iterate() {
  // Make sure propensities have been initialized
  if (initialized_ == false) {
    Initialize();
  }

  // Basic sanity checks
  if (alpha_sum_ <= 0) {
    throw std::runtime_error(
        "Gillespie: Propensity of system is 0. No reactions will execute.");
  }
  double random_num = Random::random();
  // Calculate tau, i.e. time until next reaction
  double tau = (1.0 / alpha_sum_) * std::log(1.0 / random_num);
  if (!std::isnormal(tau)) {
    throw std::underflow_error("Underflow error.");
  }
  time_ += tau;
  // Randomly select next reaction to execute, weighted by propensities
  auto next_reaction = Random::WeightedChoiceIndex(reactions_, alpha_list_);
  reactions_[next_reaction]->Execute();
  UpdatePropensity(next_reaction);
  if (reactions_[next_reaction]->remove() == true) {
    DeleteReaction(next_reaction);
  }
  iteration_++;
}

void Gillespie::Initialize() {
  // Check that propensities have not already been initialized.
  if (initialized_ == true) {
    throw std::runtime_error(
        "Gillespie: Reaction propensities have already been initialized.");
  }
  // Update all propensities
  for (int i = 0; i < reactions_.size(); i++) {
    UpdatePropensity(i);
  }
  // Make sure prop sum is starting from 0, sum over all propensities.
  alpha_sum_ = 0;
  for (const auto &alpha : alpha_list_) {
    alpha_sum_ += alpha;
  }
  initialized_ = true;
}
