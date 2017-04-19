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
 * single-moleculre tracking used for both genomes (transcription) and
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
  void Execute();
  void ShiftMask();
  void Terminate();
  void CoverElement(const std::string &species_name) {}
  void UncoverElement(const std::string &species_name) {}
  double CalculatePropensity();
  /**
   * Getters and setters
   */
  int index() { return index_; }
  double prop_sum() { return prop_sum_; }
  int uncovered(const std::string &name) { return uncovered_[name]; }
  Signal<std::string, std::string> termination_signal_;

private:
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
  /**
   * Mask corresponding to this polymer. Controls which elements are hidden.
   */
  Mask mask_;
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

  void Insert(Polymerase::Ptr pol);
  Polymerase Choose();
  void Move(const Polymerase &pol);
  void UncoverElements(const Polymerase &pol);
  void RecoverElements(const Polymerase &pol);
  bool ResolveTermination(const Polymerase &pol, const Element &elem);
  bool ResolveMaskCollisions(const Polymerase &pol);
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

#endif // SRC_POLYMER_HPP_