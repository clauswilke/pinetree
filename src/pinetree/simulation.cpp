#include <cmath>
#include <fstream>
#include <iostream>

#include "choices.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include "tracker.hpp"

Simulation::Simulation(double cell_volume) : cell_volume_(cell_volume) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();
  gillespie_ = Gillespie();
  tracker.propensity_signal_.ConnectMember(&gillespie_,
                                           &Gillespie::UpdatePropensity);
}

void Simulation::seed(int seed) { Random::seed(seed); }

void Simulation::Run(int stop_time, int time_step,
                     const std::string &output_prefix) {
  auto &tracker = SpeciesTracker::Instance();
  Initialize();
  // Set up file output streams
  std::ofstream countfile(output_prefix + "_counts.tsv", std::ios::trunc);
  // Output header
  countfile << "time\tspecies\tcount\ttranscript\tribo_density\n";
  int out_time = 0;
  while (gillespie_.time() < stop_time) {
    if ((out_time - gillespie_.time()) < 0.001) {
      countfile << tracker.GatherCounts(gillespie_.time());
      countfile.flush();
      out_time += time_step;
    }
    gillespie_.Iterate();
  }
  countfile.close();
}

void Simulation::AddReaction(double rate_constant,
                             const std::vector<std::string> &reactants,
                             const std::vector<std::string> &products) {
  auto rxn = std::make_shared<SpeciesReaction>(rate_constant, cell_volume_,
                                               reactants, products);
  auto &tracker = SpeciesTracker::Instance();
  for (const auto &reactant : reactants) {
    tracker.Add(reactant, rxn);
  }
  for (const auto &product : products) {
    tracker.Add(product, rxn);
  }
  gillespie_.LinkReaction(rxn);
}

void Simulation::AddSpecies(const std::string &name, int copy_number) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Increment(name, copy_number);
}

void Simulation::AddPolymerase(const std::string &name, int footprint,
                               double mean_speed, int copy_number) {
  auto pol = Polymerase(name, footprint, mean_speed);
  polymerases_.push_back(pol);
  AddSpecies(name, copy_number);
}

void Simulation::RegisterPolymer(Polymer::Ptr polymer) {
  // Encapsulate polymer in Bridge reaction and add to reaction list
  auto bridge = std::make_shared<Bridge>(polymer);
  gillespie_.LinkReaction(bridge);
}

void Simulation::RegisterGenome(Genome::Ptr genome) {
  RegisterPolymer(genome);
  genome->termination_signal_.ConnectMember(
      &SpeciesTracker::Instance(), &SpeciesTracker::TerminateTranscription);
  genome->transcript_signal_.ConnectMember(shared_from_this(),
                                           &Simulation::RegisterTranscript);
  genomes_.push_back(genome);
}

void Simulation::RegisterTranscript(Transcript::Ptr transcript) {
  RegisterPolymer(transcript);
  transcript->termination_signal_.ConnectMember(
      &SpeciesTracker::Instance(), &SpeciesTracker::TerminateTranslation);
}

void Simulation::Initialize() {
  if (genomes_.size() == 0) {
    std::cerr << "Warning: There are no Genome objects registered with "
                 "Simulation. Did you forget to register a Genome?"
              << std::endl;
  }
  // Create Bind reactions for each promoter-polymerase pair
  for (Genome::Ptr genome : genomes_) {
    for (auto promoter_name : genome->bindings()) {
      for (auto pol : polymerases_) {
        if (promoter_name.second.count(pol.name()) != 0) {
          double rate_constant = promoter_name.second[pol.name()];
          Polymerase pol_template = Polymerase(pol);
          auto reaction = std::make_shared<Bind>(
              rate_constant, cell_volume_, promoter_name.first, pol_template);
          auto &tracker = SpeciesTracker::Instance();
          tracker.Add(promoter_name.first, reaction);
          tracker.Add(pol.name(), reaction);
          gillespie_.LinkReaction(reaction);
        }
      }
    }
  }
}

void Simulation::CountTermination(const std::string &name) {
  auto new_name = name + "_total";
  if (terminations_.count(name) == 0) {
    terminations_[new_name] = 1;
  } else {
    terminations_[new_name]++;
  }
}
