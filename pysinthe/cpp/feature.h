/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef CPP_FEATURE_H_  // header guard
#define CPP_FEATURE_H_

#include <string>
#include <vector>
/**
 * A generic feature in or on `Polymer`. Designed to be extended
 * by `Terminator`, `Promoter`, `Polymerase`, etc.
 *
 * TODO: make abstract?by `Terminator`, `Promoter`, `Polymerase`, etc.
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
    std::vector<std::string> interactions);
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

 private:
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


#endif  // CPP_FEATURE_H_
