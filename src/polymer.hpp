/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_POLYMER_HPP_ // header guard
#define SRC_POLYMER_HPP_

#include <map>
#include <string>
#include <vector>

#include "feature.hpp"

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
class Polymer {
public:
  /**
    * The only constructor for Polymer.
    *
    * @param name name of this polymer (should it be unique?)
    * @param length length of polymer (used purely for debugging, may not be
    *     needed)
    * @param elements all elements on this polymer, including promoters,
    *     terminators, etc.
    * @param mask mask object which determines which portions of the polymer
    *     are currently inaccessible
    */
  Polymer(const std::string &name, int start, int stop,
          const Element::VecPtr &elements, const Mask &mask);
  /**
   * Bind a polymerase object to the polymer. Randomly select an open
   * promoter with which to bind and update the polymerases position to the
   * position of that promoter.
   *
   * @param pol polymerase object (pointer)
   * @param promoter_name the name of a promoter that pol will bind
   */
  void Bind(Polymerase::Ptr pol, const std::string &promoter_name);
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
  void Move(Polymerase::Ptr pol);
  /**
   * Shift mask by 1 base-pair and check for uncovered elements.
   */
  void ShiftMask();
  /**
   * Terminate polymerization reaction, fire the appropriate signals,
   * reduce the total propensity of this polymer, and then delete the
   * polymerase object itself.
   *
   * @param pol polymerase object to be deleted
   * @param last_gene last gene polymerase passed over before terminating
   */
  void Terminate(Polymerase::Ptr pol, const std::string &last_gene);
  /**
   * Update the cached count of uncovered promoters/elements.
   * TODO: Make private and refactor covering signals
   *
   * @param species_name name of species to cover
   */
  void CoverElement(const std::string &species_name);
  /**
   * Update the cached count of uncovered promoters/elements.
   * TODO: Make private and refactor covering signals
   *
   * @param species_name name of species to uncover
   */
  void UncoverElement(const std::string &species_name);
  /**
   * Getters and setters
   */
  int index() { return index_; }
  /**
   * There are two getters for prop_sum_... mostly to maintain the interface
   * that Simulation expects.
   */
  double CalculatePropensity() { return prop_sum_; }
  double prop_sum() { return prop_sum_; }
  int uncovered(const std::string &name) { return uncovered_[name]; }
  int start() const { return start_; }
  int stop() const { return stop_; }
  /**
   * Signal to fire when a polymerase terminates.
   */
  Signal<int, std::string, std::string> termination_signal_;

protected:
  /**
   * An index for this polymer, used by Simulation.
   */
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
  Polymerase::VecPtr polymerases_;
  /**
   * Vector of elements on this polymer (Promoters, Terminators, etc.)
   */
  Element::VecPtr elements_;

public:
  /**
   * Mask corresponding to this polymer. Controls which elements are hidden.
   */
  Mask mask_;

protected:
  /**
   * Cached total propensity of this polymerase, i.e. the sum of all of the
   * polymerase speeds.
   */
  double prop_sum_;
  /**
   * List of all propensities for polymerases on this polymer.
   */
  std::vector<double> prop_list_;
  /**
   * Cached count of uncovered elements on this polymer, used by Simulation.
   */
  std::map<std::string, int> uncovered_;
  /**
   * Add a polymerase to polymerase list, while maintaining the
   * order in which polymerases currently on the polymer. Higher
   * indices correspond to downstream polymerases, and lower
   * indices correspond to upstream polymerases.
   *
   * @param pol pointer to polymerase object
   */
  void Insert(Polymerase::Ptr pol);
  /**
   * Randomly select next polymerase to move, weighted by propensity
   * (i.e. speed).
   *
   * @return pointer to selected polymerase object
   */
  Polymerase::Ptr Choose();
  /**
   * Temporarily uncover all elements covered by a given polymerase.
   *
   * @param pol pointer to polymerase object
   */
  void UncoverElements(Polymerase::Ptr pol);
  /**
   * Recover all elements that a given polymerase should be covering,
   * and trigger actions if there has been a change in state (for example,
   * from covered to uncovered).
   *
   * @param pol pointer to polymerase object
   */
  void RecoverElements(Polymerase::Ptr pol);
  /**
   * Determine if polymerase should terminate upon interacting with a terminator
   * or if it should read through the terminator.
   *
   * @param pol pointer to polymerase object
   * @param element pointer to terminator object
   *
   * @return true if polymerase is terminating
   */
  bool ResolveTermination(Polymerase::Ptr pol, Element::Ptr element);
  /**
   * Check for collisions between polymerase and this polymer's mask.
   *
   * @param pol pointer to polymerase object
   *
   * @return true if this pol will collide with mask (but not shift mask)
   */
  bool ResolveMaskCollisions(Polymerase::Ptr pol);
  /**
   * Check for collisions between polymerases.
   *
   * @param pol pointer to polymerase object
   *
   * @return true if polymerases will collide
   */
  bool ResolveCollisions(Polymerase::Ptr pol);
  /**
   * Do two elements intersect?
   *
   * @param elem1 an element
   * @param elem2 another element
   *
   * @return true if elements intersect
   */
  bool Intersect(const Feature &elem1, const Feature &elem2);
};

class Transcript : public Polymer {
public:
  Transcript(const std::string &name, int start, int stop,
             const Element::VecPtr &elements, const Mask &mask);
  typedef std::shared_ptr<Transcript> Ptr;
  void Bind(Polymerase::Ptr pol, const std::string &promoter_name);
  void ShiftMask() { Polymer::ShiftMask(); }
};

/**
 * Track polymeraes on DNA, deal with collisions, promoters, terminators, and
 * constructing transcripts. Inherits from Polymer. Unlike Polymer, Genome must
 * * construct a transcript upon promoter binding.
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
  Genome(const std::string &name, int length, const Element::VecPtr &elements,
         const Element::VecPtr &transcript_template, const Mask &mask);
  /**
   * Bind a polymerase to genome and construct new transcript.
   *
   * @param pol pointer to polymerase to bind
   * @param promoter name of promoter to which this polymerase binds
   */
  void Bind(Polymerase::Ptr pol, const std::string &promoter_name);
  Signal<Transcript::Ptr> transcript_signal_;

private:
  Element::VecPtr transcript_template_;
  /**
   * Build a transcript object corresponding to start and stop positions within
   * this genome.
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

#endif // SRC_POLYMER_HPP_