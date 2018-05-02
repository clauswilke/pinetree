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
  if (index >= reactions_.size() || index < 0) {
    throw std::range_error(
        "Gillespie: Reaction index out of range for reaction deletion.");
  }
  // Update alpha sum
  alpha_sum_ -= reactions_[index]->CalculatePropensity();
  // Remove from alpha list
  alpha_list_.erase(alpha_list_.begin() + index);
  // Remove from reactions list
  reactions_.erase(reactions_.begin() + index);
}

void Gillespie::UpdatePropensity(Reaction::Ptr reaction) {
  // if (index >= reactions_.size() || index >= alpha_list_.size() || index < 0)
  // {
  //  throw std::range_error(
  //      "Gillespie: Reaction index out of range for propensity update.");
  // }
  // double new_prop = reactions_[index]->CalculatePropensity();
  // double diff = new_prop - alpha_list_[index];
  // alpha_sum_ += diff;
  // alpha_list_[index] = new_prop;
  double alpha_diff = reaction->CalculatePropensity();
  auto it = std::find(reactions_.begin(), reactions_.end(), reaction);
  if (it != reactions_.end()) {
    auto index = std::distance(reactions_.begin(), it);
    alpha_list_[index] += alpha_diff;
  } else {
    // Don't throw an error unless everything has been initialized
    if (initialized_ == true) {
      throw std::runtime_error(
          "Attempting to update propensity of invalid reaction.");
    }
  }
  alpha_sum_ += alpha_diff;
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
  // std::cout << std::to_string(alpha_list_[next_reaction]) << std::endl;
  UpdatePropensity(reactions_[next_reaction]);
  if (reactions_[next_reaction]->remove() == true) {
    // std::cout << std::to_string(alpha_list_[next_reaction]) << std::endl;
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
    UpdatePropensity(reactions_[i]);
  }
  // Make sure prop sum is starting from 0, sum over all propensities.
  // alpha_sum_ = 0;
  // for (const auto &alpha : alpha_list_) {
  //  alpha_sum_ += alpha;
  // }
  initialized_ = true;
}
