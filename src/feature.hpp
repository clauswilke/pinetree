/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_FEATURE_HPP_ // header guard
#define SRC_FEATURE_HPP_

#include <string>
#include <vector>

#include "event_signal.hpp"
/**
 * A generic feature in or on `Polymer`. Designed to be extended
 * by `Terminator`, `Promoter`, `Polymerase`, etc.
 *
 * TODO: make abstract?
 */
class Feature
{
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
  std::string const &name() const { return name_; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  void set_start(int start) { start_ = start; }
  void set_stop(int stop) { stop_ = stop; }
  /**
   * Check to see if some other feature interacts with this feature.
   * @param  name name of a feature object
   * @return      TRUE if features interact
   */
  bool CheckInteraction(const std::string &name);

protected:
  /**
   * Name of this feature.
   */
  std::string name_;
  /**
   * The type of this feature (e.g. "promoter", "terminator", etc.)
   */
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
   * Vector of names of other features/polymerases that this feature interacts
   * with.
   */
  std::vector<std::string> interactions_;
};

/**
 * A molecule that binds to `Polymer` and moves.
 */
class Polymerase
{
public:
  /**
   * The only constructor for Polymerase.
   * @param name name of polymerase (unique?)
   * @param footprint polymerase footprint
   * @param speed speed of polymerase
   */
  Polymerase(const std::string &name, int footprint, int speed);
  typedef std::shared_ptr<Polymerase> Ptr;
  bool operator<(Polymerase const &other_pol) const;
  /**
   * Getters and setters.
   */
  int start() const { return start_; }
  int stop() const { return stop_; }
  double speed() const { return speed_; }
  /**
   * Move one position forward.
   */
  void Move();
  /**
   * Move one positioin back.
   */
  void MoveBack();

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
  double speed_;
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
class Mask : public Feature
{
public:
  /**
   * Only constructor for Mask.
   */
  Mask(const std::string &name, int start, int stop,
       const std::vector<std::string> &interactions);
  /**
   * Shift mask backwards one position.
   */
  void Recede() { start_++; }
};

/**
 * A fixed feature in the polymer that can be covered or uncovered.
 */
class Element : public Feature
{
public:
  /**
   * Only constructor of Element.
   */
  Element(const std::string &name, int start, int stop,
          const std::vector<std::string> &interactions);
  typedef std::shared_ptr<Element> Ptr;
  /**
   * Save covering state.
   */
  void SaveState() { old_covered_ = covered_; }
  /**
   * Was this element just uncovered?
   * @return True if element was just uncovered.
   */
  bool WasUncovered() { return old_covered_ >= 1; }
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
  void Uncover()
  {
    if (covered_ > 0)
      covered_--;
  }
  /**
   * Is this element covered at all?
   * @return True if at least one feature is covering element.
   */
  bool IsCovered() { return covered_ > 0; }
  /**
   * Check for change in state and react appropriately. (Should be overridden
   * by children).
   */
  virtual void CheckState() {}
  /**
   * Signal to fire when element changes state from uncovered to covered
   */
  Signal<const std::string &> cover_signal_;
  /**
   * Signal to fire when element changes state from covered to uncovered
   */
  Signal<const std::string &> uncover_signal_;

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

/**
 * A promoter class
 */
class Promoter : public Element
{
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
           const std::vector<std::string> &interactions);
  typedef std::shared_ptr<Promoter> Ptr;
  /**
   * Check to see if covering state has changed and fire appropriate signals.
   */
  void CheckState();
};

class Terminator : public Element
{
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
             const std::vector<std::string> &interactions);
  typedef std::shared_ptr<Terminator> Ptr;
  /**
   * Check for changes in covering state and fire appropriate signals.
   */
  void CheckState();
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
  int reading_frame() { return reading_frame_; }
  void set_reading_frame(int reading_frame) { reading_frame_ = reading_frame; }
  bool readthrough() { return readthrough_; }
  void set_readthrough(bool readthrough) { readthrough_ = readthrough; }

private:
  /**
   * Name of gene that this terminator terminates on. This is the value that 
   * will get reported to the species tracker.
   */
  std::string gene_;
  /**
   * Readthrough state of terminator. True if a polymerase is reading through it
   * and false otherwise.
   */
  bool readthrough_;
  /**
   * Reading frame for terminator.
   */
  int reading_frame_;
};

#endif // SRC_FEATURE_HPP_
