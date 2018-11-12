#include <cmath>
#include <fstream>
#include <iostream>

#include "choices.hpp"
#include "model.hpp"
#include "polymer.hpp"
#include "tracker.hpp"

Model::Model(double cell_volume) : cell_volume_(cell_volume) {
  auto &tracker = SpeciesTracker::Instance();
  tracker.Clear();
  gillespie_ = Gillespie();
  tracker.propensity_signal_.ConnectMember(&gillespie_,
                                           &Gillespie::UpdatePropensity);
}

void Model::seed(int seed) { Random::seed(seed); }

void Model::Simulate(int time_limit, int time_step,
                     const std::string &output = "counts.tsv") {
  auto &tracker = SpeciesTracker::Instance();
  Initialize();
  // Set up file output streams
  std::ofstream countfile(output, std::ios::trunc);
  // Output header
  countfile << "time\tspecies\tprotein\ttranscript\tribo_density\n";
  int out_time = 0;
  while (gillespie_.time() < time_limit) {
    if ((out_time - gillespie_.time()) < 0.001) {
      countfile << tracker.GatherCounts(gillespie_.time());
      countfile.flush();
      out_time += time_step;
    }
    gillespie_.Iterate();
  }
  countfile.close();
}

void Model::AddReaction(double rate_constant,
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

void Model::AddSpecies(const std::string &name, int copy_number) {
  if (name.substr(0, 2) == "__") {
    throw std::invalid_argument(
        "Names prefixed with '__' (double underscore) are reserved for "
        "internal use.");
  }
  auto &tracker = SpeciesTracker::Instance();
  tracker.Increment(name, copy_number);
}

void Model::AddPolymerase(const std::string &name, int footprint,
                          double mean_speed, int copy_number) {
  auto pol = Polymerase(name, footprint, mean_speed);
  polymerases_.push_back(pol);
  auto &tracker = SpeciesTracker::Instance();
  tracker.Increment(name, copy_number);
}

void Model::AddRibosome(int footprint, double mean_speed, int copy_number) {
  auto pol = Polymerase("__ribosome", footprint, mean_speed);
  polymerases_.push_back(pol);
  auto &tracker = SpeciesTracker::Instance();
  tracker.Increment("__ribosome", copy_number);
}

void Model::RegisterPolymer(Polymer::Ptr polymer) {
  // Encapsulate polymer in PolymerWrapper reaction and add to reaction list
  auto wrapper = std::make_shared<PolymerWrapper>(polymer);
  polymer->wrapper(wrapper);
  gillespie_.LinkReaction(wrapper);
}

void Model::RegisterGenome(Genome::Ptr genome) {
  RegisterPolymer(genome);
  genome->termination_signal_.ConnectMember(
      &SpeciesTracker::Instance(), &SpeciesTracker::TerminateTranscription);
  genome->transcript_signal_.ConnectMember(shared_from_this(),
                                           &Model::RegisterTranscript);
  genomes_.push_back(genome);
}

void Model::RegisterTranscript(Transcript::Ptr transcript) {
  RegisterPolymer(transcript);
  transcript->termination_signal_.ConnectMember(
      &SpeciesTracker::Instance(), &SpeciesTracker::TerminateTranslation);
}

void Model::Initialize() {
  if (genomes_.size() == 0) {
    std::cerr << "Warning: There are no Genome objects registered with "
                 "Model. Did you forget to register a Genome?"
              << std::endl;
  }
  // Create Bind reactions for each promoter-polymerase pair
  for (Genome::Ptr genome : genomes_) {
    for (auto promoter_name : genome->bindings()) {
      for (auto pol : polymerases_) {
        if (promoter_name.second.count(pol.name()) != 0) {
          double rate_constant = promoter_name.second[pol.name()];
          Polymerase pol_template = Polymerase(pol);
          auto reaction = std::make_shared<BindPolymerase>(
              rate_constant, cell_volume_, promoter_name.first, pol_template);
          auto &tracker = SpeciesTracker::Instance();
          tracker.Add(promoter_name.first, reaction);
          tracker.Add(pol.name(), reaction);
          gillespie_.LinkReaction(reaction);
        }
      }
    }
    if (genome->transcript_degradation_rate() != 0.0) {
      // TODO: user defined Rnase speed
      // auto rnase_template = Rnase(10, 30);
      auto rnase_template =
          Rnase(genome->rnase_footprint(), genome->rnase_speed());
      auto reaction = std::make_shared<BindRnase>(
          genome->transcript_degradation_rate(), cell_volume_, rnase_template,
          "__rnase_site");
      auto &tracker = SpeciesTracker::Instance();
      tracker.Add("__rnase_site", reaction);
      // tracker.Add("__rnase", reaction);
      gillespie_.LinkReaction(reaction);

      auto rnase_template_ext =
          Rnase(genome->rnase_footprint(), genome->rnase_speed());
      auto reaction_ext = std::make_shared<BindRnase>(
          genome->transcript_degradation_rate_ext(), cell_volume_,
          rnase_template_ext, "__rnase_site_ext");
      tracker.Add("__rnase_site_ext", reaction_ext);
      // tracker.Add("__rnase", reaction_ext);
      gillespie_.LinkReaction(reaction_ext);
    }
  }
}

void Model::CountTermination(const std::string &name) {
  auto new_name = name + "_total";
  if (terminations_.count(name) == 0) {
    terminations_[new_name] = 1;
  } else {
    terminations_[new_name]++;
  }
}
