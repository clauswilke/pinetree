/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_FEATURE_HPP_  // header guard
#define SRC_FEATURE_HPP_

#include <iostream>
#include <string>
#include <vector>

#include "event_signal.hpp"

/**
 * Abstract class from which all fixed elements on a polymer inherit. These
 * include promoters, terminators, ribosome binding sites, and stop codons.
 * They all share a common interface for tracking whether they're covered
 * or uncovered.
 */
class FixedElement : public std::enable_shared_from_this<FixedElement> {
 public:
  /**
   * Only constructor for FixedElement.
   *
   * @param name name of the element
   * @param start start position of element on polymer
   * @param stop stop position of element on polymer
   * @param interactions map of names of MobileElements that this element
   *                     interacts with, and the rate constant of that
   *                     interaction
   */
  FixedElement(const std::string &name, int start, int stop,
               const std::map<std::string, double> &interactions);
  /**
   * Forces child classes to define a destructor.
   */
  virtual ~FixedElement() = 0;
  /**
   * Save covering state.
   */
  void ResetState() { old_covered_ = covered_; }
  /**
   * Was this element just uncovered?
   * @return True if element was just uncovered.
   */
  bool WasUncovered() { return old_covered_ >= 1 and covered_ == 0; }
  /**
   * Was this element just covered?
   * @return True if element was just covered.
   */
  bool WasCovered() { return old_covered_ == 0 && covered_ > 0; }
  /**
   * Cover this element. Elements can be covered by multiple features.
   */
  void Cover() { covered_++; }
  /**
   * Uncover element.
   */
  void Uncover() {
    if (covered_ > 0) {
      covered_ = covered_ - 1;
    }
  }
  /**
   * Is this element covered at all?
   * @return True if at least one feature is covering element.
   */
  bool IsCovered() { return covered_ > 0; }
  /**
   * Getters and setters
   */
  const std::string &gene() const { return gene_; }
  void gene(const std::string &gene) { gene_ = gene; }
  std::string const &name() const { return name_; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  int reading_frame() const { return reading_frame_; }
  void set_reading_frame(int reading_frame) { reading_frame_ = reading_frame; }

 protected:
  /**
   * Name of this feature.
   */
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
  std::map<std::string, double> interactions_;
  /**
   * Name of gene that this terminator terminates on. This is the value that
   * will get reported to the species tracker.
   */
  std::string gene_;
  /**
   * Count of how many features are currently covering this element.
   */
  int covered_;
  /**
   * Used to cache old covering count to then test for changes in state.
   */
  int old_covered_;
  /**
   * Reading frame for terminator.
   */
  int reading_frame_;
};

/**
 * A promoter class
 */
class Promoter : public FixedElement {
 public:
  /**
   * The only constructor for Promoter.
   *
   * @param name name of promoter
   * @param start start position of promoter
   * @param stop stop position of promoter (also transcription start site)
   * @param interactions vector of polymerases that interact with this promoter
   */
  Promoter(const std::string &name, int start, int stop,
           const std::map<std::string, double> &interactions);
  /**
   * Promoter does not create or accept any new resources (i.e. pointers).
   */
  ~Promoter(){};
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<Promoter> Ptr;
  typedef std::vector<std::shared_ptr<Promoter>> VecPtr;
  /**
   * Check to see if covering state has changed and fire appropriate signals.
   */
  Promoter::Ptr Clone() const;

  bool CheckInteraction(const std::string &name);

  bool first_exposure() { return first_exposure_; }
  void first_exposure(bool first_exposure) { first_exposure_ = first_exposure; }

 private:
  /**
   * Has the site been exposed before?
   */
  bool first_exposure_;
};

class Terminator : public FixedElement {
 public:
  /**
   * Only constructor for Terminator.
   *
   * @param name name of terminator
   * @param start start position of terminator
   * @param stop stop position of terminator
   * @param interactions list of features that this terminator interacts with
   */
  Terminator(const std::string &name, int start, int stop,
             const std::map<std::string, double> &interactions);
  /**
   * Terminator does not create or accept any new resources (i.e. pointers).
   */
  ~Terminator(){};
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<Terminator> Ptr;
  typedef std::vector<std::shared_ptr<Terminator>> VecPtr;

  Terminator::Ptr Clone() const;
  /**
   * Check to see if feature interacts with this terminator and is in the
   * correct reading frame.
   *
   * @param name name of other feature
   * @param reading_frame reading frame of of interacting feature
   *
   * @return bool true if feature interacts with terminator
   */
  bool CheckInteraction(const std::string &name, int reading_frame);
  /**
   * Getters and setters
   */
  bool readthrough() const { return readthrough_; }
  void set_readthrough(bool readthrough) { readthrough_ = readthrough; }
  double efficiency(const std::string &pol_name) {
    return interactions_[pol_name];
  }

 private:
  /**
   * Readthrough state of terminator. True if a polymerase is reading through it
   * and false otherwise.
   */
  bool readthrough_;
};

class MobileElement : public std::enable_shared_from_this<MobileElement> {
 public:
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<MobileElement> Ptr;
  typedef std::vector<std::shared_ptr<MobileElement>> VecPtr;

  MobileElement(const std::string &name, int footprint, int speed);
  std::string const &name() const { return name_; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  void set_start(int start) { start_ = start; }
  void set_stop(int stop) { stop_ = stop; }
  double speed() const { return speed_; }
  int footprint() const { return footprint_; }
  /**
   * Getters and setters.
   */
  int reading_frame() const { return reading_frame_; }
  void set_reading_frame(int reading_frame) { reading_frame_ = reading_frame; }
  /**
   * Move one position forward.
   */
  virtual void Move() = 0;
  /**
   * Move one positioin back.
   */
  virtual void MoveBack() = 0;

 protected:
  /**
   * Name of this feature.
   */
  std::string name_;
  std::string type_;
  /**
   * The start site of the feature. Usually the most upstream site position.
   */
  int start_;
  /**
   * The stop site of the feature. Usually the most downstream site position.
   */
  int stop_;
  /**
   * Foot print in base pairs of polymerase.
   */
  int footprint_;
  /**
   * Speed in bp/s.
   */
  double speed_;
  /**
   * Reading frame of polymerase (0, 1, or 2).
   */
  int reading_frame_;
};

/**
 * A molecule that binds to `Polymer` and moves.
 */
class Polymerase : public MobileElement {
 public:
  /**
   * The only constructor for Polymerase.
   * @param name name of polymerase (unique?)
   * @param footprint polymerase footprint
   * @param speed speed of polymerase
   */
  Polymerase(const std::string &name, int footprint, int speed);
  /**
   * Some typedefs to make code less verbose
   */
  typedef std::shared_ptr<Polymerase> Ptr;
  typedef std::vector<std::shared_ptr<Polymerase>> VecPtr;
  /**
   * Move one position forward.
   */
  void Move();
  /**
   * Move one positioin back.
   */
  void MoveBack();
};

/**
 * A pseudo-feature that tracks which portion of a genome or polymer are not
 * yet accessible. For example, as the genome is entering the cell, or as a
 * transcript is being synthesized.
 */
class Mask : public MobileElement {
 public:
  /**
   * Only constructor for Mask.
   */
  Mask(int start, int stop, const std::map<std::string, double> &interactions);
  /**
   * Shift mask backwards one position.
   */
  void Move() { start_++; }
  void MoveBack() { start_--; }

  bool CheckInteraction(const std::string &name);

 private:
  /**
   * Vector of names of other features/polymerases that this feature interacts
   * with.
   */
  std::map<std::string, double> interactions_;
};

/**
 * A polymerase-like object that degrades RNA.
 */
class Rnase : public MobileElement {
 public:
  Rnase(int footprint, int speed);
  void Move() { stop_++; }
  void MoveBack() { stop_--; }
};

#endif  // SRC_FEATURE_HPP_
