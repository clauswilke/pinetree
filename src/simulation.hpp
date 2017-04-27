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
 * Tracks species' copy numbers and maintains promoter-to-polymer and species-
 * to-reaction maps to easily look up which polymers contain a given promoter
 * and which reactions involve a given species. These maps are needed to cache
 * propensities and increase the performance of the simulation.
 *
 * TODO: Move propensity cache from Simulation into this class?
 */
namespace SpeciesTracker {
/**
 * Some convenience typedefs.
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
 * Species counts.
 */
static std::map<std::string, int> species_;
/**
 * Promoter-to-polymer map.
 */
static std::map<std::string, Polymer::VecPtr> promoter_map_;
/**
 * Species-to-reaction map.
 */
static std::map<std::string, Reaction::VecPtr> species_map_;
/**
 * Getters and setters
 */
int species(const std::string &reactant);
/**
 * Signal to fire when propensity needs to be updated.
 */
static Signal<int> propensity_signal_;
};

#endif // header guard
