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
  transcripts_.clear();
  ribo_per_transcript_.clear();
  propensity_signal_.DisconnectAll();
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
      propensity_signal_.Emit(reaction);
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
  Increment(species_name, 0);
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

void SpeciesTracker::Remove(const std::string &promoter_name,
                            Polymer::Ptr polymer) {
  if (promoter_map_.count(promoter_name) != 0) {
    auto it = std::find(promoter_map_[promoter_name].begin(),
                        promoter_map_[promoter_name].end(), polymer);
    if (it != promoter_map_[promoter_name].end()) {
      promoter_map_[promoter_name].erase(it);
    }
  }
}

void SpeciesTracker::TerminateTranscription(
    std::shared_ptr<PolymerWrapper> wrapper, const std::string &pol_name,
    const std::string &gene_name) {
  Increment(pol_name, 1);
  propensity_signal_.Emit(wrapper);
  // CountTermination("transcript");
}

void SpeciesTracker::TerminateTranslation(
    std::shared_ptr<PolymerWrapper> wrapper, const std::string &pol_name,
    const std::string &gene_name) {
  Increment(pol_name, 1);
  Increment(gene_name, 1);
  IncrementRibo(gene_name, -1);
  propensity_signal_.Emit(wrapper);
  // CountTermination(gene_name);
}

const Reaction::VecPtr &SpeciesTracker::FindReactions(
    const std::string &species_name) {
  if (species_map_.count(species_name) == 0) {
    throw std::runtime_error("Species not found in tracker.");
  }
  return species_map_[species_name];
}

const Polymer::VecPtr &SpeciesTracker::FindPolymers(
    const std::string &promoter_name) {
  if (promoter_map_.count(promoter_name) == 0) {
    throw std::runtime_error("Species not found in tracker.");
  }
  return promoter_map_[promoter_name];
}

int SpeciesTracker::species(const std::string &reactant) {
  if (species_.count(reactant) == 0) {
    throw std::runtime_error("Species not found in tracker.");
  }
  return species_[reactant];
}

int SpeciesTracker::transcripts(const std::string &transcript_name) {
  return transcripts_[transcript_name];
}

int SpeciesTracker::ribo_per_transcript(const std::string &transcript_name) {
  return ribo_per_transcript_[transcript_name];
}

const std::string SpeciesTracker::GatherCounts(double time_stamp) {
  std::map<std::string, std::vector<double>> output;
  for (auto elem : species_) {
    output[elem.first].push_back(elem.second);
    output[elem.first].push_back(0);
    output[elem.first].push_back(0);
  }
  for (auto transcript : transcripts_) {
    if (output[transcript.first].size() == 3) {
      output[transcript.first][1] = transcript.second;
    } else {
      output[transcript.first].push_back(0);
      output[transcript.first].push_back(transcript.second);
      output[transcript.first].push_back(0);
    }
    if (ribo_per_transcript_.find(transcript.first) !=
        ribo_per_transcript_.end()) {
      output[transcript.first][2] =
          double(ribo_per_transcript_[transcript.first]) /
          double(output[transcript.first][1]);
    }
  }
  std::string out_string;
  for (auto row : output) {
    out_string = out_string + (std::to_string(time_stamp) + "\t" + row.first +
                               "\t" + std::to_string(row.second[0]) + "\t" +
                               std::to_string(row.second[1]) + "\t" +
                               std::to_string(row.second[2]) + "\n");
  }
  return out_string;
}