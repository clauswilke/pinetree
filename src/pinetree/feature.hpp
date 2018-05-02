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
  void reading_frame(int reading_frame) { reading_frame_ = reading_frame; }
  bool first_exposure() const { return first_exposure_; }
  void first_exposure(bool first_exposure) { first_exposure_ = first_exposure; }

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
   * Name of gene associated with this FixedElement. This is the value that
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
   * Reading frame for FixedElement.
   */
  int reading_frame_;
  /**
   * Has the site been exposed before?
   */
  bool first_exposure_;
};

/**
 * A BindingSite class for both promoters and ribosome binding sites
 */
class BindingSite : public FixedElement {
 public:
  /**
   * The only constructor for BindingSite.
   *
   * @param name name of BindingSite
   * @param start start position of BindingSite
   * @param stop stop position of BindingSite (also transcription/translation
   *             start site)
   * @param interactions name--binding-strength map of MobileElements that
   *                     interact with this BindingSite
   */
  BindingSite(const std::string &name, int start, int stop,
              const std::map<std::string, double> &interactions);
  /**
   * BindingSite does not create or accept any new resources (i.e. pointers).
   */
  ~BindingSite(){};
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<BindingSite> Ptr;
  typedef std::vector<std::shared_ptr<BindingSite>> VecPtr;
  /**
   * Create a deep copy of BindingSite. Used by Polymer when creating
   * transcripts from a transcript template.
   *
   * @return std::shared_ptr<BindingSite> pointer to deep copy of BindingSite
   */
  BindingSite::Ptr Clone() const;
  /**
   * Check to see if BindingSite interacts with the MobileElement (name).
   *
   * @param name name of MobileElement
   *
   * @return bool true if MobileElement interacts with BindingSite
   */
  bool CheckInteraction(const std::string &name);
  /**
   * Mark this site as degraded.
   */
  void Degrade();
  bool degraded() { return degraded_; }

 private:
  /**
   * Has this site been degraded? (I.e. covered by RNAse)
   */
  bool degraded_ = false;
};

/**
 * ReleaseSite class for terminators and stop codons
 */
class ReleaseSite : public FixedElement {
 public:
  /**
   * Only constructor for ReleaseSite.
   *
   * @param name name of ReleaseSite
   * @param start start position of ReleaseSite
   * @param stop stop position of ReleaseSite
   * @param interactions list of features that this ReleaseSite interacts with
   */
  ReleaseSite(const std::string &name, int start, int stop,
              const std::map<std::string, double> &interactions);
  /**
   * ReleaseSite does not create or accept any new resources (i.e. pointers).
   */
  ~ReleaseSite(){};
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<ReleaseSite> Ptr;
  typedef std::vector<std::shared_ptr<ReleaseSite>> VecPtr;
  /**
   * Create a deep copy of ReleaseSite. Used by Polymer when creating
   * transcripts from a transcript template.
   *
   * @return std::shared_ptr<BindingSite> pointer to deep copy of BindingSite
   */
  ReleaseSite::Ptr Clone() const;
  /**
   * Check to see if feature interacts with this ReleaseSite and is in the
   * correct reading frame.
   *
   * @param name name of other feature
   * @param reading_frame reading frame of of interacting feature
   *
   * @return bool true if feature interacts with ReleaseSite
   */
  bool CheckInteraction(const std::string &name, int reading_frame);
  /**
   * Getters and setters
   */
  bool readthrough() const { return readthrough_; }
  void readthrough(bool readthrough) { readthrough_ = readthrough; }
  double efficiency(const std::string &pol_name) {
    return interactions_[pol_name];
  }

 private:
  /**
   * Readthrough state of ReleaseSite. True if a polymerase is reading through
   * it and false otherwise.
   */
  bool readthrough_;
};

/**
 * Abstract class for all elements capable of moving on a polymer. This includes
 * polymerase, ribosomes, masks, and RNases.
 */
class MobileElement : public std::enable_shared_from_this<MobileElement> {
 public:
  /**
   * The only constructor for MobileElement.
   *
   * @param name name of this mobile element
   * @param footprint footpring of this mobile element in base pairs
   * @param speed average speed of this polymerase in bp/s
   */
  MobileElement(const std::string &name, int footprint, int speed);
  /**
   * Force child classes to explicitly define a destructor.
   */
  virtual ~MobileElement() = 0;
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<MobileElement> Ptr;
  typedef std::vector<std::shared_ptr<MobileElement>> VecPtr;
  /**
   * Move one position forward.
   */
  virtual void Move() = 0;
  /**
   * Move one positioin back.
   */
  virtual void MoveBack() = 0;
  /**
   * Getters and setters.
   */
  std::string const &name() const { return name_; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  void start(int start) { start_ = start; }
  void stop(int stop) { stop_ = stop; }
  double speed() const { return speed_; }
  int footprint() const { return footprint_; }
  int reading_frame() const { return reading_frame_; }
  void reading_frame(int reading_frame) { reading_frame_ = reading_frame; }
  std::string gene_bound() const { return gene_bound_; }
  void gene_bound(std::string gene) { gene_bound_ = gene; }

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
  /**
   * Which gene did this MobileElement originally bind to?
   * (used for ribosomes)
   */
  std::string gene_bound_ = "";
};

/**
 * A molecule that binds to `Polymer` and moves. This may be either a ribosome
 * or a polymerase.
 */
class Polymerase : public MobileElement {
 public:
  /**
   * The only constructor for Polymerase.
   *
   * @param name name of polymerase (unique?)
   * @param footprint polymerase footprint
   * @param speed speed of polymerase
   */
  Polymerase(const std::string &name, int footprint, int speed);
  /**
   * Polymerase does not create or accept any new resources (i.e. pointers)
   */
  ~Polymerase(){};
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
   * Mask does not create or accept any new resources (i.e. pointers)
   */
  ~Mask(){};
  /**
   * Shift start position of mask forwards one position (uncovering more of
   * the polymer).
   */
  void Move() {
    start_++;
    footprint_--;
  }
  /**
   * Shift start position of mask backward one position, covering more of the
   * polymer
   */
  void MoveBack();
  /**
   * Does this polymerase interact with this mask?
   *
   * @param name name of polymerase
   *
   * @return bool true if elements interact
   */
  bool CheckInteraction(const std::string &name);

 private:
  /**
   * Vector of names of other features/polymerases that this feature interacts
   * with.
   */
  std::map<std::string, double> interactions_;
};

/**
 * A polymerase-like object that degrades RNA from 5' to 3' end.
 */
class Rnase : public MobileElement {
 public:
  /**
   * The only constructor for Rnase
   *
   * @param footprint footprint of Rnase
   * @param speed average speed that Rnase degrades polymer
   */
  Rnase(int footprint, int speed);
  /**
   * Rnase does not create or accept any new resources (i.e. pointers)
   */
  ~Rnase(){};
  /**
   * Move the Rnase one step forward (i.e. degrade more polymer)
   */
  void Move() {
    stop_++;
    footprint_++;
  }
  /**
   * Move Rnase back.
   */
  void MoveBack() {
    stop_--;
    footprint_--;
  }
};

#endif  // SRC_FEATURE_HPP_
