/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_POLYMER_HPP_ // header guard
#define SRC_POLYMER_HPP_

#include <vector>
#include <map>
#include <string>

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
class Polymer
{
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
            const std::vector<Element> &elements, const Mask &mask);
    void bind_polymerase();
    void execute();
    void shift_mask();
    void terminate();
    int count_uncovered(const std::string &species_name);
    void cover_element(const std::string &species_name) {}
    void uncover_element(const std::string &species_name);
    double calculate_propensity();
    // Gallant::Signal2<std::string, std::string> termination_signal;

    int get_index() { return index_; }

  private:
    int index_;
    std::string name_;
    int start_;
    int stop_;
    std::vector<Polymerase> polymerases_;
    std::vector<Element> elements_;
    Mask mask_;
    double prop_sum_;
    std::vector<double> prop_list_;
    std::map<std::string, int> uncovered_;

    void insert_polymerase_(const Polymerase &pol);
    Polymerase &choose_polymerase_();
    void move_polymerase_(const Polymerase &pol);
    void uncover_elements_(const Polymerase &pol);
    void recover_elements_(const Polymerase &pol);
    bool resolve_termination_(const Polymerase &pol, const Element &elem);
    bool resolve_mask_collisions_(const Polymerase &pol);
    bool elements_intersect_(const Element &elem1, const Element &elem2);
};

#endif // SRC_POLYMER_HPP_