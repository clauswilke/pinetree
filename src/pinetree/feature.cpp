/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "feature.hpp"
#include "tracker.hpp"

Feature::Feature(const std::string &name, int start, int stop,
                 const std::map<std::string, double> &interactions)
    : name_(name), start_(start), stop_(stop), interactions_(interactions) {}

bool Feature::CheckInteraction(const std::string &name) {
  return interactions_.count(name);
}

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : Feature(name, 0, footprint - 1, std::map<std::string, double>()),
      footprint_(footprint),
      speed_(speed) {
  left_most_element_ = 0;
  bound_ = 0;
  type_ = "polymerase";
  reading_frame_ = -1;
}

void Polymerase::Move() {
  start_++;
  stop_++;
}

void Polymerase::MoveBack() {
  start_--;
  stop_--;
}

Mask::Mask(const std::string &name, int start, int stop,
           const std::map<std::string, double> &interactions)
    : Feature(name, start, stop, interactions) {}

Element::Element(const std::string &name, int start, int stop,
                 const std::map<std::string, double> &interactions)
    : Feature(name, start, stop, interactions), covered_(0), old_covered_(0) {
  first_exposure_ = false;
}

Promoter::Promoter(const std::string &name, int start, int stop,
                   const std::map<std::string, double> &interactions)
    : Element(name, start, stop, interactions) {
  type_ = "promoter";
}

void Promoter::CheckState() {
  if (WasCovered()) {
    cover_signal_.Emit(name_);
  } else if (WasUncovered()) {
    uncover_signal_.Emit(name_);
    if (!first_exposure_ && interactions().count("ribosome") == 1) {
      SpeciesTracker::Instance().IncrementTranscript(gene_, 1);
      first_exposure_ = true;
    }
  }
}

Promoter::Ptr Promoter::Clone() const {
  return std::make_shared<Promoter>(*this);
}

Terminator::Terminator(const std::string &name, int start, int stop,
                       const std::map<std::string, double> &interactions)
    : Element(name, start, stop, interactions) {
  readthrough_ = false;
  reading_frame_ = -1;
  type_ = "terminator";
}

void Terminator::CheckState() {
  if (WasUncovered()) {
    readthrough_ = false;
    uncover_signal_.Emit(name_);
  } else if (WasCovered()) {
    cover_signal_.Emit(name_);
  }
}

Terminator::Ptr Terminator::Clone() const {
  return std::make_shared<Terminator>(*this);
}

bool Terminator::CheckInteraction(const std::string &name, int reading_frame) {
  return (reading_frame == reading_frame_ && Element::CheckInteraction(name));
}
