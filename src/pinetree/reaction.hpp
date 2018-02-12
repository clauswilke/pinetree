/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_REACTION_HPP  // header guard
#define SRC_REACTION_HPP

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
   * Some getters and setters.
   */
  virtual int index() const { return index_; }
  virtual void index(int index) { index_ = index; }

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
   * @param rate_constant macroscopic rate constant (SpeciesReaction will
   *                      convert rate constant from macroscopic to mesoscopic)
   * @param reactants vector of reactant names
   * @param products vector of product names
   * @param volume the volume in which these reactions will occur
   *
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
  Bridge(Polymer::Ptr polymer);
  /**
   * Retrieve total propensity of all reactions tha tmay occur within this
   * polymer.
   *
   * @return total propensity of reactions within polymer
   */
  double CalculatePropensity() { return polymer_->prop_sum(); }
  /**
   * Execute reaction within polymer (e.g. typically moving a polymerase)
   */
  void Execute();

  void index(int index);
  int index() const { return index_; }

 private:
  /**
   * Pointer to polymer object that this reaction encapsulates.
   */
  Polymer::Ptr polymer_;
  void ReportChanges();
};

#endif  // header guard