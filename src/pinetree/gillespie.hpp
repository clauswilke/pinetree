#ifndef SRC_GILLESPIE_HPP  // header guard
#define SRC_GILLESPIE_HPP

#include <vector>

#include "reaction.hpp"

class Gillespie {
 public:
  /**
   * Add Reaction object to reaction queue.
   */
  void LinkReaction(Reaction::Ptr reaction);
  /**
   * Remove Reaction object from reaction queue.
   */
  void DeleteReaction(int index);
  /**
   * Update propensity of a reaction.
   */
  void UpdatePropensity(Reaction::Ptr reaction);
  /**
   * Execute one iteration of the gillespie algorithm.
   */
  void Iterate();
  /**
   * Getters and setters.
   */
  double time() { return time_; }

 private:
  /**
   * True if Initialize() has been called.
   */
  bool initialized_ = false;
  /**
   * Current simulation time.
   */
  double time_ = 0;
  /**
   * Current simulation iteration.
   */
  int iteration_ = 0;
  /**
   * Vector of individual reaction propensities in same order as reactions_.
   */
  std::vector<double> alpha_list_;
  /**
   * Running total of propensities.
   */
  double alpha_sum_ = 0;
  /**
   * Vector of all reactions.
   */
  Reaction::VecPtr reactions_;
  /**
   * Compute all propensities after all reactions have been added.
   */
  void Initialize();
};

#endif  // header guard