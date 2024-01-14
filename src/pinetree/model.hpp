/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_SIMULATION_HPP  // header guard
#define SRC_SIMULATION_HPP

#include <memory>

#include "gillespie.hpp"
#include "polymer.hpp"
#include "reaction.hpp"

/**
 * Coordinate polymers and species-level reactions.
 */
class Model : public std::enable_shared_from_this<Model> {
 public:
  /**
   * Construct a simulation
   */
  Model(double cell_volume);
  /**
   * Run the simulation until the given time point and write output to a file.
   *
   * @param prefix for output files
   */
  void Simulate(int time_limit, double time_step, const std::string &output);
  /**
   * Set a seed for random number generator.
   */
  void seed(int seed);
  /**
   * Run the simulation with dynamic tRNAs. 
   * 
   * @param codons the codon-anticodon pairs and initial tRNA copy numbers
   * @param rate_constant rate for tRNA charging
   */ 
  // void AddtRNA(const CodonMap &codons, double rate_constant);
  void AddtRNA(std::map<std::string, std::map<std::string, std::map<std::string, int>>> &codons, double rate_constant);
  /**
   * Alternative method for adding dynamic tRNAs to the simulation. Should be used when tRNA species- specific charging
   * rates are needed and/or if considering wobble base pairings. 
   * 
   * @param codon_map the codon-anticodon pairs
   * @param counts Intial species counts for charged and uncharged tRNAs
   * @param rate_constants charging rate for each tRNA species
   */ 
  void AddtRNA(std::map<std::string, std::vector<std::string>> &codon_map, 
               std::map<std::string, std::pair<int, int>> &counts, 
               std::map<std::string, double> &rate_constants);
  /**
   * Add species to simulation.
   *
   * @param name species name
   * @param copy_number copy number of species
   */
  void AddSpecies(const std::string &name, int copy_number);
  /**
   * Add a polymerase to simulation.
   *
   * @param name polymerase name
   * @param footprint footprint (in basepairs) of polymerase on DNA
   * @param mean_speed mean speed of polymerase
   * @param copy_number copy number of free polymerase
   */
  void AddPolymerase(const std::string &name, int footprint, double speed,
                     int copy_number);
  /**
   * Add a ribosome to simulation.
   *
   * @param footprint footprint (in basepairs) of ribosome on RNA
   * @param speed mean speed of ribosome
   * @param copy_number copy number of free ribosomes
   */
  void AddRibosome(int footprint, double mean_speed, int copy_number);
  /**
   * Add a species reaction to simulation.
   *
   * @param rate_constant macroscopic rate constant of reaction
   * @param reactants vector of reactant names
   * @param products vector of product names
   */
  void AddReaction(double rate_constant,
                   const std::vector<std::string> &reactants,
                   const std::vector<std::string> &products);
  /**
   * Add a genome to the list of reactions.
   *
   * @param genome pointer to Genome object
   */
  void RegisterGenome(Genome::Ptr genome);
  /**
   * Add a transcript to the list of reactions.
   *
   * @param pointer to Transcript object
   */
  void RegisterTranscript(Transcript::Ptr transcript);
  void Initialize();
  /**
   * Record when a polymerase reaches a terminator so that we can track total
   * quantity of a species synthesized during simulation.
   *
   * TODO: Move to species tracker.
   */
  void CountTermination(const std::string &name);

 private:
  /**
   * Gillespie object
   */
  Gillespie gillespie_;
  /**
   * Vector of all genomes in this ismulation
   */
  Genome::VecPtr genomes_;
  /**
   * Vector of predefined transcripts in simulation
   */
  Transcript::VecPtr transcripts_;
  /**
   * Vector of all polymerases in simulation
   */
  std::vector<Polymerase> polymerases_;
  /**
   * Cell volume
   */
  double cell_volume_;
  /**
   * Has this model been initialized?
   */
  bool initialized_ = false;
  /**
   * Map of terminations.
   */
  std::map<std::string, int> terminations_;
  /**
   * Add a generic polymer to the list of reactions.
   *
   * @param polymer pointer to Polymer object
   */
  void RegisterPolymer(Polymer::Ptr polymer);
};

#endif  // header guard
