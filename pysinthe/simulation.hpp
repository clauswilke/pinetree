/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_SIMULATION_HPP // header guard
#define SRC_SIMULATION_HPP

#include <memory>

#include "polymer.hpp"
/**
 * An abstract reaction class. Propensity refers to the reaction propensity, or
 * the probability that this reaction will occur in the next time step. It is a
 * concept taken directly from the Gillepsie Algorithm.
 */
class Reaction : public std::enable_shared_from_this<Reaction> {
public:
  /**
   * Some convenience typedefs.
   */
  typedef std::shared_ptr<Reaction> Ptr;
  typedef std::vector<std::shared_ptr<Reaction>> VecPtr;
  /**
   * Return the propensity of this reaction.
   *
   * @return propensity of the reaction
   */
  virtual double CalculatePropensity() = 0;
  /**
   * Execute the reaction.
   */
  virtual void Execute() = 0;
  /**
   * Return the index of the reaction.
   */
  int index() const { return index_; }
  void index(int index) { index_ = index; }
  void virtual InitCounts() {}

protected:
  /**
   * The index of this reaction in the reaction list maintained by Simulation.
   */
  int index_;
};

/**
 * A generic class for a species-level reaction. It currently only supports 2 or
 * fewer reactants.
 */
class SpeciesReaction : public Reaction {
public:
  /**
   * The only constructor of SpeciesReaction.
   *
   * @param tracker pointer to SpeciesTracker object
   * @param rate_constant mesoscopic rate constant (this is different, but
   *  similar in concept, to rate constants used in ODE models)
   * @param reactants vector of reactant names
   * @param products vector of product names
   *
   * TODO: Change SpeciesTracker be a static class instead of passing it to
   * constructor?
   */
  SpeciesReaction(double rate_constant, double volume,
                  const std::vector<std::string> &reactants,
                  const std::vector<std::string> &products);
  /**
   * Convenience typedefs.
   */
  typedef std::shared_ptr<SpeciesReaction> Ptr;
  /**
   * Calculate propensity of this reaction.
   *
   * @return propensity of reaction
   */
  double CalculatePropensity();
  /**
   * Execute the reaction. Decrement reactants and increment products.
   */
  void Execute();
  /**
   * Getters and setters.
   */
  const std::vector<std::string> &reactants() const { return reactants_; }
  const std::vector<std::string> &products() const { return products_; }
  void InitCounts();

private:
  /**
   * Rate constant of reaction.
   */
  double rate_constant_;
  /**
   * Vector of reactant names.
   */
  const std::vector<std::string> reactants_;
  /**
   * Vector of product names.
   */
  const std::vector<std::string> products_;
};

/**
 * Bind a polymerase to a polymer.
 */
class Bind : public Reaction {
public:
  /**
   * The only constructor of Bind.
   *
   * @param rate_constant rate constant of the binding reaction
   * @param promoter_name name of promoter involved in this reaction
   * @param pol_template Polymerase object that will get copied and bound to
   *  Polymer upon execution of this reaction
   */
  Bind(double rate_constant, double volume, const std::string &promoter_name,
       const Polymerase &pol_template);
  /**
   * Calculate propensity of binding reaction.
   *
   * @return propensity of this reaction
   */
  double CalculatePropensity();
  /**
   * Decrement reactants, choose polymer to bind, construct a new polymerase,
   * and bind the polymerase to the polymer.
   */
  void Execute();
  void InitCounts();

private:
  /**
   * Rate constant of this reaction.
   */
  double rate_constant_;
  /**
   * Name of promoter involved in this binding reaction.
   */
  const std::string promoter_name_;
  /**
   * Name of polymerase involved in this binding reaction.
   */
  const std::string pol_name_;
  /**
   * Polymerase object to be copied and bound to Polymer upon execution.
   */
  const Polymerase pol_template_;
};

/**
 * A thin wrapper for Polymer so it can participate in species-level reaction
 * processing.
 */
class Bridge : public Reaction {
public:
  /**
   * Only constructor for Bridge.
   *
   * @param polymer pointer to polymer object that this reaction is
   *  encapsulating
   */
  Bridge(Polymer::Ptr polymer) : polymer_(polymer) {}
  /**
   * Retrieve total propensity of all reactions tha tmay occur within this
   * polymer.
   *
   * @return total propensity of reactions within polymer
   */
  double CalculatePropensity() { return polymer_->CalculatePropensity(); }
  /**
   * Execute reaction within polymer (e.g. typically moving a polymerase)
   */
  void Execute() { polymer_->Execute(); }

private:
  /**
   * Pointer to polymer object that this reaction encapsulates.
   */
  Polymer::Ptr polymer_;
};

/**
 * Coordinate polymers and species-level reactions.
 */
class Simulation : public std::enable_shared_from_this<Simulation> {
public:
  Simulation();
  /**
   * Run the simulation until the given time point and write output to a file.
   *
   * @param prefix for output files
   */
  void Run(const std::string &output_name);
  /**
   * Add a SpeciesReaction object to the list of reactions.
   *
   * @param reaction pointer to SpeciesReaction object
   */
  void RegisterReaction(Reaction::Ptr reaction);
  /**
   * Add a generic polymer to the list of reactions.
   *
   * @param polymer pointer to Polymer object
   */
  void RegisterPolymer(Polymer::Ptr polymer);
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
  /**
   * Initialize propensities. Must be called before simulation is run.
   *
   * TODO: Make private?
   */
  void InitPropensity();
  /**
   * Update propensity of a reaction at a given index.
   *
   * @index index of polymer in reaction list
   */
  void UpdatePropensity(int index);
  /**
   * Execute one iteration of gillespie algorithm.
   */
  void Execute();
  /**
   * Increment promoter count at species level.
   *
   * TODO: Possibley remove this function and call SpeciesTracker directly
   * elsewhere?
   *
   * @param species_name name of promoter to free
   */
  void FreePromoter(const std::string &species_name);
  /**
   * Decrement promoter count at species level.
   *
   * TODO: Possibly remove this function and call SpeciesTracker directly.
   *
   * @param species_name name of promoter to block
   */
  void BlockPromoter(const std::string &species_name);
  /**
   * Update propensities and species counts after transcription has terminated.
   *
   * @param polymer_index index of genome in reaction list
   * @param pol_name name of polymerase completing transcription
   * @param gene_name name of last gene on the polymerase encountered (not
   *  currently used, but may be used in future)
   */
  void TerminateTranscription(int polymer_index, const std::string &pol_name,
                              const std::string &gene_name);
  /**
  * Update propensities and species counts after translation has terminated.
  *
  * @param polymer_index index of transcript in reaction list
  * @param pol_name name of ribosome completing transcription
  * @param protein_name name of newly-synthesized protein
  */
  void TerminateTranslation(int polymer_index, const std::string &pol_name,
                            const std::string &protein_name);
  /**
   * Record when a polymerase reaches a terminator so that we can track total
   * quantity of a species synthesized during simulation.
   *
   * TODO: Move to species tracker.
   */
  void CountTermination(const std::string &name);
  /**
   * Getters and setters
   */
  void stop_time(double stop_time) { stop_time_ = stop_time; }
  double stop_time() { return stop_time_; }
  void time_step(double time_step) { time_step_ = time_step; }
  double time_step() { return time_step_; }
  double alpha_sum() { return alpha_sum_; }

private:
  /**
   * Current simulation time.
   */
  double time_;
  /**
   * Time at which simulation should stop.
   */
  double stop_time_;
  /**
   * Time step at which to report data.
   */
  double time_step_;
  /**
   * Simulation iteraction counter
   */
  int iteration_;
  /**
   * Vector of all reactions in this simulation.
   */
  Reaction::VecPtr reactions_;
  /**
   * Vector of propensity values for reactions in this simulation.
   */
  std::vector<double> alpha_list_;
  /**
   * Total simulation propensity.
   */
  double alpha_sum_;
  /**
   * Map of terminations.
   */
  std::map<std::string, int> terminations_;
};

#endif // header guard
