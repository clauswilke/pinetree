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
class Reaction {
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
  SpeciesReaction(double rate_constant,
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
  Bind(double rate_constant, const std::string &promoter_name,
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
 * Tracks species' copy numbers and maintains promoter-to-polymer and species-
 * to-reaction maps to easily look up which polymers contain a given promoter
 * and which reactions involve a given species. These maps are needed to cache
 * propensities and increase the performance of the simulation.
 *
 * Follows a singleton pattern:
 * http://stackoverflow.com/a/1008289/5360633
 *
 * TODO: Move propensity cache from Simulation into this class?
 */
class SpeciesTracker {
public:
  /**
   * Retrieve the only instance of this object.
   */
  static SpeciesTracker &Instance();
  /**
   * Clear all data in the tracker.
   */
  void Clear();
  /**
   * Delete copy constructor and assignment operator.
   */
  SpeciesTracker(SpeciesTracker const &) = delete;
  void operator=(SpeciesTracker const &) = delete;
  /**
   * Register a SpeciesReaction with the species tracker.
   */
  void Register(SpeciesReaction::Ptr reaction);
  /**
   * Change a species count by a given value (positive or negative).
   *
   * @param species_name name of species to change count
   * @param copy_number number to add to current copy number count
   */
  void Increment(const std::string &species_name, int copy_number);
  /**
   * Add a species-reaction pair to species-reaction map.
   *
   * @param species_name name of species
   * @param reaction reaction object (pointer) that involves species
   */
  void Add(const std::string &species_name, Reaction::Ptr reaction);
  /**
   * Add a promoter-polymer pair to promoter-polymer map.
   *
   * @param promter_name of promoter
   * @param polymer polymer object that contains the named promoter (pointer)
   */
  void Add(const std::string &promoter_name, Polymer::Ptr polymer);
  /**
   * Get polymers that contain a given promoter.
   *
   * @param promoter_name name of promoter
   *
   * @return vector of pointers to Polymer objects that contain promoter_name
   */
  const Polymer::VecPtr &FindPolymers(const std::string &promoter_name);
  /**
   * Get reactions that involve a given species.
   *
   * @param species_name name of species
   *
   * @return vector of pointers to Reaction objects that involve species_name
   */
  const Reaction::VecPtr &FindReactions(const std::string &species_name);
  /**
   * Getters and setters
   */
  int species(const std::string &reactant);
  /**
   * Signal to fire when propensity needs to be updated.
   */
  Signal<int> propensity_signal_;

private:
  /**
   * Private constructor for singleton
   */
  SpeciesTracker() {}
  /**
   * Species counts.
   */
  std::map<std::string, int> species_;
  /**
   * Promoter-to-polymer map.
   */
  std::map<std::string, Polymer::VecPtr> promoter_map_;
  /**
   * Species-to-reaction map.
   */
  std::map<std::string, Reaction::VecPtr> species_map_;
};

#endif // header guard
