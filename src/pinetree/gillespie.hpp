#ifndef SRC_GILLESPIE_HPP  // header guard
#define SRC_GILLESPIE_HPP

#include <vector>
#include "simulation.hpp"

class Gillespie {
 public:
  /**
   * Add Reaction object to reaction queue.
   */
  void LinkReaction(Reaction::Ptr reaction);
  /**
   * Compute all propensities after all reactions have been added.
   */
  void Initialize();
  /**
   * Update propensity of a reaction at a given index.
   */
  void UpdatePropensity(int index);
  /**
   * Execute one iteration of the gillespie algorithm.
   */
  void Iterate();

 private:
  /**
   * True if Initialize() has been called.
   */
  bool initialized_;
  /**
   * Current simulation time.
   */
  double time_;
  /**
   * Current simulation iteration.
   */
  int iteration_;
  /**
   * Vector of individual reaction propensities in same order as reactions_.
   */
  std::vector<double> alpha_list_;
  /**
   * Running total of propensities.
   */
  double alpha_sum_;
  /**
   * Vector of all reactions.
   */
  Reaction::Ptr reactions_;
};

#endif  // header guard