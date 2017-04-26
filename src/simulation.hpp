/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#ifndef SRC_SIMULATION_HPP // header guard
#define SRC_SIMULATION_HPP

#include <memory>

#include "polymer.hpp"

class SpeciesTracker {
public:
  SpeciesTracker();
  typedef std::shared_ptr<SpeciesTracker> Ptr;
  void Increment(const std::string &species_name, int copy_number);
  void Add(const std::string &species_name, Reaction::Ptr);
  void Add(const std::string &promoter_name, Polymer::Ptr);
  const Polymer::VecPtr &FindPolymers(const std::string &promoter_name);
  const Reaction::VecPtr &FindReactions(const std::string &species_name);
  Signal<> propensity_signal_;

private:
  std::map<std::string, int> species_;
  std::map<std::string, Polymer::Ptr> promoter_map_;
  std::map<std::string, Reaction::Ptr> species_map_;
};

class Reaction {
public:
  typedef std::shared_ptr<Reaction> Ptr;
  typedef std::vector<std::shared_ptr<Reaction>> VecPtr;
  Reaction();
  virtual double CalculatePropensity() = 0;
  virtual void Execute() = 0;
};

class SpeciesReaction : public Reaction {
public:
  SpeciesReaction(SpeciesTracker::Ptr tracker, double rate_constant,
                  const std::vector<std::string> &reactants,
                  const std::vector<std::string> &products);
  double CalculatePropensity();
  void Execute();
};

#endif // header guard
