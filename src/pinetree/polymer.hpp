/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_POLYMER_HPP_  // header guard
#define SRC_POLYMER_HPP_

#include <map>
#include <string>
#include <vector>

#include "IntervalTree.h"
#include "feature.hpp"

/**
 * Hack-y forward declaration.
 */
class Polymer;
class PolymerWrapper;
class Reaction;

/**
 * Manages all MobileElements (e.g., polymerases and ribosomes) on a Polymer.
 * MobileElements are maintained in order. It also tracks any polymers
 * (transcripts) that may be attached to a polymerase. This class also tracks
 * the total propensity of all the MobileElements.
 */
class MobileElementManager {
 public:
  /**
   * Only constructor of MobileElementManager
   *
   * @param weights Base-pair specific movement weights.
   */
  MobileElementManager(const std::vector<double> &weights);
  /**
   * Insert an MobileElement-Polymer pair while maintaining order of
   * MobileElements
   *
   * @param pol polymerase
   * @param polymer polymer that may be attached to MobileElement
   */
  void Insert(std::shared_ptr<MobileElement> pol,
              std::shared_ptr<Polymer> polymer);
  /**
   * Delete a polymerase by index.
   *
   * @param index Index of polymerase-polymer pair
   */
  void Delete(int index);
  /**
   * Return a randomly selected MobileElement, weighted by speed and base-pair
   * specific weights.
   *
   * @return index of MobileElement-Polymer pair
   */
  int Choose();
  /**
   * Is this a valid index?
   *
   * @param index Index of MobileElement-Polymer pair
   * @return true if index is valid
   */
  bool ValidIndex(int index) { return index < polymerases_.size(); };
  /**
   * Get a Polymer at a given index.
   *
   * @return pointer to MobileElement
   */
  std::shared_ptr<MobileElement> GetPol(int index);
  /**
   * Get an attached Polymer at a given index.
   *
   * @param index of MobileElement-Polymer pair
   * @return pointer to Polymer
   */
  std::shared_ptr<Polymer> GetAttached(int index);
  /**
   * Update movement propensity of MobileElement at a given index.
   *
   * @param index Index of MobileElement-Polymer pair
   */
  void UpdatePropensity(int index);
  /**
   * Getters and setters.
   */
  double prop_sum() { return prop_sum_; }
  int pol_count() { return pol_count_; }

 private:
  /**
   * Total propensity sum of all MobileElements being tracked
   */
  double prop_sum_ = 0;
  /**
   * Total number of polymerases (anything except RNases). Used to make sure
   * that there are no active polymerases left on the polymer before degrading
   * polymer.
   */
  int pol_count_ = 0;
  /**
   * Vector of propensities corresponding to each MobileElement-Polymer pair
   */
  std::vector<double> prop_list_;
  /**
   * MobileElement-Polymer pairs
   */
  std::vector<
      std::pair<std::shared_ptr<MobileElement>, std::shared_ptr<Polymer>>>
      polymerases_;
  /**
   * Base-pair specific movement weights.
   */
  std::vector<double> weights_;
};

/**
 * Track element objects, polymerase objects, and collisions on a single
 * polymer. Move polymerase objects along the polymer. Handle logic for
 * covering and uncovering of elements. This class contains the core of the
 * single-molecule tracking used for both genomes (transcription) and
 * transcripts (translation).
 *
 * The terms polymer, polymerase, promoter, and terminator are all used
 * generically in this class. Each term could refer to a different biological
 *  definition in the context of transcription and translation.
 *
 * * polymer: genome, transcript
 * * polymerase: RNA polymerase, ribosome, any object that binds to polymer
 * * promoter: promoter, ribosome binding site, any site on a polymer in which
 *      a protein (polymerase) can bind
 * * terminator: terminator, stop codon, any site on polymer that ends
 *      polymerization
 */
class Polymer : public std::enable_shared_from_this<Polymer> {
 public:
  /**
   * The most used constructor for Polymer.
   *
   * @param name name of this polymer (should it be unique?)
   * @param length length of polymer (used purely for debugging, may not be
   *     needed)
   * @param elements all elements on this polymer, including promoters,
   *     terminators, etc.
   * @param mask mask object which determines which portions of the polymer
   *     are currently inaccessible
   */
  Polymer(const std::string &name, int start, int stop);
  /**
   * Remove from promoter-polymer lap. Error checking to make sure this
   * polymer is no longer linked to a polymerase. Make sure there are no
   * polymerases on here.
   */
  ~Polymer();
  /**
   * De-register polymer from promoter-polymer map in SpeciesTracker.
   */
  void Unlink();
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<Polymer> Ptr;
  typedef std::vector<std::shared_ptr<Polymer>> VecPtr;
  /**
   * Create interval trees.
   * Make sure elements covered by mask have correct state.
   */
  virtual void Initialize();
  /**
   * Bind a polymerase object to the polymer. Randomly select an open
   * promoter with which to bind and update the polymerases position to the
   * position of that promoter.
   *
   * @param pol polymerase object (pointer)
   * @param promoter_name the name of a promoter that pol will bind
   */
  virtual void Bind(MobileElement::Ptr pol, const std::string &promoter_name);
  /**
   * Select a polymerase to move next and deal with terminations.
   */
  void Execute();
  /**
   * Move polymerase and deal with collisions and covering/uncovering of
   * elements.
   *
   * This method is optimized to take advantage of the fact that we only ever
   * need to look at the elements that this polymerase is covering, and then
   * one element ahead on the polymer, and one element behind after the
   * polymerase has moved.
   *
   * It probably shouldn't be public, but the logic in here is complex and
   * needs to be unit tested.
   *
   * @param pol polymerase to move
   */
  void Move(int pol_index);
  /**
   * Shift mask by 1 base-pair and check for uncovered elements.
   */
  virtual void ShiftMask();

  /**
   * Getters and setters. There are two getters for prop_sum_... mostly to
   * maintain the interface that Model expects.
   */
  void index(int index) { index_ = index; }
  int index() { return index_; }
  double prop_sum() { return polymerases_.prop_sum(); }
  int uncovered(const std::string &name) { return uncovered_[name]; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  bool degrade() { return degrade_; }
  bool attached() { return attached_; }
  void attached(bool attached) { attached_ = attached; }
  void wrapper(std::shared_ptr<PolymerWrapper> wrapper) { wrapper_ = wrapper; }
  std::shared_ptr<PolymerWrapper> wrapper() { return wrapper_.lock(); }

  /**
   * Signal to fire when a polymerase terminates.
   */
  Signal<std::shared_ptr<PolymerWrapper>, const std::string &,
         const std::string &>
      termination_signal_;

 protected:
  std::weak_ptr<PolymerWrapper> wrapper_;
  int index_;
  /**
   * Name of polymer
   */
  std::string name_;
  /**
   * Start position.
   */
  int start_;
  /**
   * End position
   */
  int stop_;
  /**
   * Vector of polymerases currently on this polymer.
   */
  MobileElementManager polymerases_;
  /**
   * Total count of promoters/terminators on polymer.
   */
  int total_elements_;
  /**
   * Running count of promoters or terminators that have been marked as
   * degraded.
   */
  int degraded_elements_;
  /**
   * Should this polymer be degraded?
   */
  bool degrade_ = false;
  /**
   * Is this polymer still attached to whatever polymerase is generating it?
   */
  bool attached_ = false;

  /**
   * Vector of binding site intervals (start/stop positions)
   */
  std::vector<Interval<BindingSite::Ptr>> binding_intervals_;
  /**
   * Vector of release site intervals
   */
  std::vector<Interval<ReleaseSite::Ptr>> release_intervals_;
  /**
   * Interval tree of binding sites
   */
  IntervalTree<BindingSite::Ptr> binding_sites_;
  /**
   * Interval tree of release sites
   */
  IntervalTree<ReleaseSite::Ptr> release_sites_;
  /**
   * Mask corresponding to this polymer. Controls which elements are hidden.
   */
  Mask mask_ = Mask(0, 0, std::map<std::string, double>());
  /**
   * Cached count of uncovered elements on this polymer, used by Model.
   */
  std::map<std::string, int> uncovered_;
  /**
   * Vector of the same length as this polymer, containing weights for different
   * positions along the polymer. When a polymerase passes over a given position
   * in the genome, the weight * speed of polymerase will determine the
   * propensity for the next movement of that polymerase.
   */
  std::vector<double> weights_;
  /**
   * Finding which binding site (promoter) that the polymerase should bind to.
   *
   * @param pol Pointer to polymerase object
   * @param promoter_name Name of promoter as a string
   */
  BindingSite::Ptr FindBindingSite(MobileElement::Ptr pol,
                                   const std::string &promoter_name);
  /**
   * Attach a polymerase to the polymer.
   *
   * @param pol Pointer to polymerase
   */
  virtual void Attach(MobileElement::Ptr pol);
  /**
   * Check downstream of polymerase for any interactions and respond
   * accordingly.
   *
   * @param old_stop Prior stop position of polymerase
   * @param new_stop Current stop position of polymerase
   */
  void CheckAhead(int old_stop, int new_stop);
  /**
   * Check downstream of Rnase for potential interactions.
   *
   * @param old_stop Prior stop position of Rnase
   * @param new_stop Current stop position of Rnase
   */
  void CheckAheadRnase(int old_stop, int new_stop);
  /**
   * Check upstream of polymerase for any changes/newly uncovered elements
   *
   * @param old_start Prior start position of polymerase
   * @param new_start Current start position polymerase
   */
  void CheckBehind(int old_start, int new_start);
  /**
   * Determine if polymerase should terminate upon interacting with a terminator
   * or if it should read through the terminator.
   *
   * @param pol pointer to polymerase object
   * @param element pointer to terminator object
   *
   * @return true if polymerase is terminating
   */
  bool CheckTermination(int pol_index);
  /**
   * Check for collisions between polymerase and this polymer's mask.
   *
   * @param pol pointer to polymerase object
   *
   * @return true if this pol will collide with mask (but not shift mask)
   */
  bool CheckMaskCollisions(MobileElement::Ptr pol);
  /**
   * Check for collisions between polymerases.
   *
   * @param pol pointer to polymerase object
   *
   * @return true if polymerases will collide
   */
  bool CheckPolCollisions(int pol_index);
  /**
   * Update the cached count of uncovered promoters/elements.
   *
   * @param species_name name of species to cover
   */
  void LogCover(const std::string &species_name);
  /**
   * Update the cached count of uncovered promoters/elements.
   *
   * @param species_name name of species to uncover
   */
  void LogUncover(const std::string &species_name);
};

/**
 * An mRNA transcript. Tracks ribosomes and protein production. Only differs
 * from Polymer in capability to receive signasl from a moving polymerase on a
 * Genome and uncover the appropriate, "newly-synthesized" RNA elements.
 */
class Transcript : public Polymer {
 public:
  /**
   * The only constructor of Transcript.
   *
   * @param name name of transcript
   * @param start start position of transcript (in genomic coordinates)
   * @param stop stop position of transcript (in genomeic coordinates)
   * @param rbs_intervals vector of RBS intervals (stop/start positions)
   * @param stop_site_invervals vector of stop codon intervals (start/stop
   *        positions)
   * @param mask mask object that gets shifted as polymerase "synthesizes" more
   *  of the transcript
   */
  Transcript(const std::string &name, int start, int stop,
             const std::vector<Interval<BindingSite::Ptr>> &rbs_intervals,
             const std::vector<Interval<ReleaseSite::Ptr>> &stop_site_intervals,
             const Mask &mask, const std::vector<double> &weights);
  /**
   * Convenience typdefs.
   */
  typedef std::shared_ptr<Transcript> Ptr;
  /**
   * Bind ribosome to this transcript.
   *
   * @param pol polymerase (ribosome) object pointer
   * @param promoter_name name of RBS to bind (typically just "rbs")
   */
  void Bind(MobileElement::Ptr pol, const std::string &promoter_name);
  /**
   * Hack-y redeclaration of ShiftMask so that Signal doesn't complain about
   * mismatched types.
   */
  void ShiftMask() { Polymer::ShiftMask(); }
};

/**
 * Track polymeraes on DNA, deal with collisions, promoters, terminators, and
 * constructing transcripts. Inherits from Polymer. Unlike Polymer, Genome must
 * construct a transcript upon promoter binding.
 */
class Genome : public Polymer {
 public:
  /**
   * Genome's only constructor.
   *
   * @param name name of this genome
   * @param length length of genome (do we still need this ?)
   * @param elements DNA elements
   * @param transcript_template vector of all possible elements that could
   *  be producted by this genome (i. e. the largest possible polycistronic)
   *  transcript)
   * @param mask polymer mask (i.e. protion of genome that has not yet entered
   *  the cell)
   */
  Genome(const std::string &name, int length,
         double transcript_degradation_rate,
         double transcript_degradation_rate_ext, double rnase_speed,
         double rnase_footprint);
  void Initialize();
  void AddMask(int start, const std::vector<std::string> &interactions);
  void AddPromoter(const std::string &name, int start, int stop,
                   const std::map<std::string, double> &interactions);
  void AddTerminator(const std::string &name, int start, int stop,
                     const std::map<std::string, double> &efficiency);
  void AddGene(const std::string &name, int start, int stop, int rbs_start,
               int rbs_stop, double rbs_strength);
  void AddRnaseSite(int start, int stop);
  void AddWeights(const std::vector<double> &transcript_weights);
  const std::map<std::string, std::map<std::string, double>> &bindings();
  const double &transcript_degradation_rate() {
    return transcript_degradation_rate_;
  }
  const double &transcript_degradation_rate_ext() {
    return transcript_degradation_rate_ext_;
  }
  const double &rnase_speed() { return rnase_speed_; }
  int rnase_footprint() { return rnase_footprint_; }
  /**
   * Convenience typedefs
   */
  typedef std::shared_ptr<Genome> Ptr;
  typedef std::vector<std::shared_ptr<Genome>> VecPtr;
  /**
   * Bind a polymerase to genome and construct new transcript.
   *
   * @param pol pointer to polymerase to bind
   * @param promoter name of promoter to which this polymerase binds
   */
  void Attach(MobileElement::Ptr pol);
  Signal<Transcript::Ptr> transcript_signal_;

 private:
  std::vector<Interval<BindingSite::Ptr>> transcript_rbs_intervals_;
  std::vector<Interval<ReleaseSite::Ptr>> transcript_stop_site_intervals_;
  IntervalTree<BindingSite::Ptr> transcript_rbs_;
  IntervalTree<ReleaseSite::Ptr> transcript_stop_sites_;
  std::vector<double> transcript_weights_;
  std::map<std::string, std::map<std::string, double>> bindings_;
  double transcript_degradation_rate_ = 0.0;
  double transcript_degradation_rate_ext_ = 0.0;
  double rnase_speed_ = 0.0;
  int rnase_footprint_ = 0;
  /**
   * Build a transcript object corresponding to start and stop positions
   * within this genome.
   *
   * NOTE: Assumes that elements are already ordered by start position.
   *
   * @param start start position of transcript within genome
   * @param stop stop position of transcript within genome
   *
   * @returns pointer to Transcript object
   */
  Transcript::Ptr BuildTranscript(int start, int stop);
};

#endif  // SRC_POLYMER_HPP_