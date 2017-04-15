/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_FEATURE_HPP_  // header guard
#define SRC_FEATURE_HPP_

#include <string>
#include <vector>

#include "my_signal.h"
/**
 * A generic feature in or on `Polymer`. Designed to be extended
 * by `Terminator`, `Promoter`, `Polymerase`, etc.
 *
 * TODO: make abstract?
 */
class Feature {
 public:
  /**
   * Default constructor of Feature
   * @param name a (unique?) name for this feature
   * @param start start position of this feature within polymer
   * @param stop stop position of this feature within polymer
   * @param interactions list of names of other features that this feature
   *                     interacts with
   */
  Feature(const std::string &name, int start, int stop,
    const std::vector<std::string> &interactions);
  /**
   * Getters and setters.
   */
  int get_start() { return start_; }
  int get_stop() { return stop_; }
  void set_start(int start) { start_ = start; }
  void set_stop(int stop) { stop_ = stop; }
  /**
   * Check to see if some other feature interacts with this feature.
   * @param  name name of a feature object
   * @return      TRUE if features interact
   */
  bool check_interaction(const std::string &name);

 protected:
  std::string name_;
  /**
   * The start site of the feature. Usually the most upstream site position.
   */
  int start_;
  /**
   * The stop site of the feature. Usually the most downstream site position.
   */
  int stop_;
  /**
   * Vector of names of other features/polymerases that this feature interacts
   * with.
   */
  std::vector<std::string> interactions_;
};

/**
 * A molecule that binds to `Polymer` and moves.
 */
class Polymerase {
 public:
  /**
   * The only constructor for Polymerase.
   * @param name name of polymerase (unique?)
   * @param footprint polymerase footprint
   * @param speed speed of polymerase
   */
  Polymerase(const std::string &name, int footprint, int speed);
  /**
   * Getters and setters.
   */
  int get_start() { return start_; }
  int get_stop() { return stop_; }
  /**
   * Move one position forward.
   */
  void move();
  /**
   * Move one positioin back.
   */
  void move_back();

  // Gallant::Signal1<std::string> move_signal;
  // Signal release_signal;

 private:
  /**
   * Name of polymerase (for determining interactions)
   */
  std::string name_;
  /**
   * Start position of polymerase.
   */
  int start_;
  /**
   * End position of polymerase.
   */
  int stop_;
  /**
   * Foot print in base pairs of polymerase.
   */
  int footprint_;
  /**
   * Speed in bp/s.
   */
  int speed_;
  /**
   * Index of left-most interacting element; used for increasing efficiency of
   * polymerase movement.
   */
  int left_most_element_;
  /**
   * Where did this polymerase bind?
   */
  int bound_;
  /**
   * Type of polymerase. (Is this used?)
   */
  std::string type_;
  /**
   * Reading frame of polymerase (0, 1, or 2).
   */
  int reading_frame_;
};

/**
 * A pseudo-feature that tracks which portion of a genome or polymer are not
 * yet accessible. For example, as the genome is entering the cell, or as a
 * transcript is being synthesized.
 */
class Mask : public Feature {
 public:
  /**
   * Only constructor for Mask.
   */
  Mask(const std::string &name, int start, int stop,
    const std::vector<std::string> &interactions);
  /**
   * Shift mask backwards one position.
   */
  void recede() { start_++; }
};

/**
 * A fixed feature in the polymer that can be covered or uncovered.
 */
class Element : public Feature {
 public:
  /**
   * Only constructor of Element.
   */
  Element(const std::string &name, int start, int stop,
    const std::vector<std::string> &interactions);
  /**
   * Save covering state.
   */
  void save_state() { old_covered_ = covered_; }
  /**
   * Was this element just uncovered?
   * @return True if element was just uncovered.
   */
  bool was_uncovered() { return old_covered_ >= 1; }
  /**
   * Was this element just covered?
   * @return True if element was just covered.
   */
  bool was_covered() { return old_covered_ == 0 && covered_ > 0; }
  /**
   * Cover this element. Elements can be covered by multiple features.
   */
  void cover() { covered_++; }
  /**
   * Uncover element.
   */
  void uncover() { if (covered_ > 0) covered_--; }
  /**
   * Is this element covered at all?
   * @return True if at least one feature is covering element.
   */
  bool is_covered() { return covered_ > 0; }
  /**
   * Check for change in state and react appropriately. (Should be overridden
   * by children).
   */
  void check_state();

 private:
  /**
   * Count of how many features are currently covering this element.
   */
  int covered_;
  /**
   * Used to cache old covering count to then test for changes in state.
   */
  int old_covered_;
};

#endif  // SRC_FEATURE_HPP_
