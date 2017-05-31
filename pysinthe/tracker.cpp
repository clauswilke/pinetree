#include <algorithm>

#include "tracker.hpp"

SpeciesTracker &SpeciesTracker::Instance() {
  static SpeciesTracker instance;
  return instance;
}

void SpeciesTracker::Clear() {
  species_.clear();
  promoter_map_.clear();
  species_map_.clear();
  propensity_signal_ = Signal<int>();
}

void SpeciesTracker::Register(SpeciesReaction::Ptr reaction) {
  // Add this reaction to SpeciesTracker
  for (const auto &reactant : reaction->reactants()) {
    Add(reactant, reaction);
    // Initialize reactant counts (this usually gets set to some non-zero
    // value later)
    Increment(reactant, 0);
  }
  // Do the same for products
  for (const auto &product : reaction->products()) {
    Add(product, reaction);
    Increment(product, 0);
  }
}

void SpeciesTracker::Increment(const std::string &species_name,
                               int copy_number) {
  if (species_.count(species_name) == 0) {
    species_[species_name] = copy_number;
  } else {
    species_[species_name] += copy_number;
  }
  if (species_map_.count(species_name) > 0) {
    for (const auto &reaction : species_map_[species_name]) {
      propensity_signal_.Emit(reaction->index());
    }
  }
  if (species_[species_name] < 0) {
    throw std::runtime_error("Species count less than 0." + species_name);
  }
}

void SpeciesTracker::IncrementRibo(const std::string &transcript_name,
                                   int copy_number) {
  if (ribo_per_transcript_.count(transcript_name) == 0) {
    ribo_per_transcript_[transcript_name] = copy_number;
  } else {
    ribo_per_transcript_[transcript_name] += copy_number;
  }
  if (ribo_per_transcript_[transcript_name] < 0) {
    throw std::runtime_error("Ribosome count less than 0." + transcript_name);
  }
}

void SpeciesTracker::IncrementTranscript(const std::string &transcript_name,
                                         int copy_number) {
  if (transcripts_.count(transcript_name) == 0) {
    transcripts_[transcript_name] = copy_number;
  } else {
    transcripts_[transcript_name] += copy_number;
  }
  if (transcripts_[transcript_name] < 0) {
    throw std::runtime_error("Transcript count less than 0." + transcript_name);
  }
}

void SpeciesTracker::Add(const std::string &species_name,
                         Reaction::Ptr reaction) {
  if (species_map_.count(species_name) == 0) {
    species_map_[species_name] = Reaction::VecPtr{reaction};
  } else {
    // TODO: Maybe use a better data type here like a set?
    auto it = std::find(species_map_[species_name].begin(),
                        species_map_[species_name].end(), reaction);
    if (it == species_map_[species_name].end()) {
      species_map_[species_name].push_back(reaction);
    }
  }
}

void SpeciesTracker::Add(const std::string &promoter_name,
                         Polymer::Ptr polymer) {
  if (promoter_map_.count(promoter_name) == 0) {
    promoter_map_[promoter_name] = Polymer::VecPtr{polymer};
  } else {
    promoter_map_[promoter_name].push_back(polymer);
  }
}

const Reaction::VecPtr &
SpeciesTracker::FindReactions(const std::string &species_name) {
  return species_map_[species_name];
}

const Polymer::VecPtr &
SpeciesTracker::FindPolymers(const std::string &promoter_name) {
  return promoter_map_[promoter_name];
}

int SpeciesTracker::species(const std::string &reactant) {
  return species_[reactant];
}

int SpeciesTracker::transcripts(const std::string &transcript_name) {
  return transcripts_[transcript_name];
}