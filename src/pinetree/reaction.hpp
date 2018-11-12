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
   * Should this reaction be removed from the reaction queue?
   *
   * @return true if reaction should be removed
   */
  bool remove() { return remove_; }
  /**
   * Some getters and setters.
   */
  virtual int index() const { return index_; }
  virtual void index(int index) { index_ = index; }

 protected:
  /**
   * The index of this reaction in the reaction list maintained by Model.
   */
  int index_;

  double old_prop_ = 0;
  /**
   * Flag to mark reaction for removal.
   */
  bool remove_ = false;
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
   * Calculate *change* in propensity of this reaction.
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
 * Bind a mobile element to a polymer.
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
  Bind(double rate_constant, double volume, const std::string &promoter_name);
  /**
   * Calculate propensity of binding reaction.
   *
   * @return propensity of this reaction
   */
  virtual double CalculatePropensity() = 0;
  /**
   * Decrement reactants, choose polymer to bind, construct a new polymerase,
   * and bind the polymerase to the polymer.
   */
  virtual void Execute() = 0;
  /**
   * Randomly select a polymer to bind that has an open promoter/binding site.
   *
   * @return pointer to polymer
   */
  Polymer::Ptr ChoosePolymer();

 protected:
  /**
   * Rate constant of this reaction.
   */
  double rate_constant_;
  /**
   * Name of promoter involved in this binding reaction.
   */
  const std::string promoter_name_;
};

/**
 * Bind a Polymerase to a polymer.
 */
class BindPolymerase : public Bind {
 public:
  /**
   * Only constructor for BindPolymerase.
   *
   * @param rate_constant binding rate constant
   * @param volume volume that in which reaction occurs
   * @param pol_template polymerase object to construct upon binding
   */
  BindPolymerase(double rate_constant, double volume,
                 const std::string &promoter_name,
                 const Polymerase &pol_template);
  /**
   * Bind the polymerase.
   */
  void Execute();
  /**
   * Calculate the propensity of binding occurring.
   */
  double CalculatePropensity();

 private:
  /**
   * Polymerase object to be copied and bound to Polymer upon execution.
   */
  const Polymerase pol_template_;
};

/**
 * Bind an RNase to a polymer.
 */
class BindRnase : public Bind {
 public:
  /**
   * Only constructor of BindRnase.
   *
   * @param rate_constant rate constant of binding reaction
   * @param volume volume in which reaction occurs
   * @param rnase_template Rnase to construct upon binding
   */
  BindRnase(double rate_constant, double volume, const Rnase &rnase_template,
            const std::string &name);
  /**
   * Bind Rnase to open binding site.
   */
  void Execute();
  /**
   * Calculate propensity of binding reaction occurring.
   */
  double CalculatePropensity();

 private:
  /**
   * Polymerase object to be copied and bound to Polymer upon execution.
   */
  const Rnase pol_template_;
};

/**
 * A thin wrapper for Polymer so it can participate in species-level reaction
 * processing.
 */
class PolymerWrapper : public Reaction {
 public:
  /**
   * Only constructor for PolymerWrapper.
   *
   * @param polymer pointer to polymer object that this reaction is
   *  encapsulating
   */
  PolymerWrapper(Polymer::Ptr polymer);
  /**
   * Make sure polymer that this object is wrapping is properly destroyed.
   */
  ~PolymerWrapper();
  /**
   * Retrieve total propensity of all reactions tha tmay occur within this
   * polymer.
   *
   * @return total propensity of reactions within polymer
   */
  double CalculatePropensity() {
    if (remove_ == true) {
      old_prop_ = 0;
    }
    double new_prop = polymer_->prop_sum() - old_prop_;
    old_prop_ = polymer_->prop_sum();
    return new_prop;
  }
  /**
   * Execute reaction within polymer (e.g. typically moving a polymerase)
   */
  void Execute();
  /**
   * Getters and setters
   */
  void index(int index);
  int index() const { return index_; }

 private:
  /**
   * Pointer to polymer object that this reaction encapsulates.
   */
  Polymer::Ptr polymer_;
};

#endif  // header guard