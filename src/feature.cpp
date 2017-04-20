/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "feature.hpp"

Feature::Feature(const std::string &name, int start, int stop,
                 const std::vector<std::string> &interactions)
    : name_(name), start_(start), stop_(stop), interactions_(interactions) {}

bool Feature::CheckInteraction(const std::string &name) {
  return std::find(interactions_.begin(), interactions_.end(), name) !=
         interactions_.end();
}

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : Feature(name, 0, footprint - 1, std::vector<std::string>()),
      footprint_(footprint), speed_(speed) {
  left_most_element_ = 0;
  bound_ = 0;
  type_ = "polymerase";
  reading_frame_ = 0;
}

bool Polymerase::operator<(const Polymerase &other_pol) const {
  return start_ < other_pol.start();
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
           const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions) {}

Element::Element(const std::string &name, int start, int stop,
                 const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions), covered_(0), old_covered_(0) {}

Promoter::Promoter(const std::string &name, int start, int stop,
                   const std::vector<std::string> &interactions)
    : Element(name, start, stop, interactions) {
  type_ = "promoter";
}

void Promoter::CheckState() {
  if (WasCovered()) {
    cover_signal_.Emit(name_);
  } else if (WasUncovered()) {
    uncover_signal_.Emit(name_);
  }
}

Terminator::Terminator(const std::string &name, int start, int stop,
                       const std::vector<std::string> &interactions,
                       const std::map<std::string, double> &efficiency)
    : Element(name, start, stop, interactions), efficiency_(efficiency) {
  readthrough_ = false;
  reading_frame_ = -1;
}

void Terminator::CheckState() {
  if (WasUncovered()) {
    readthrough_ = false;
    uncover_signal_.Emit(name_);
  } else if (WasCovered()) {
    cover_signal_.Emit(name_);
  }
}

bool Terminator::CheckInteraction(const std::string &name, int reading_frame) {
  return (reading_frame == reading_frame_ && Element::CheckInteraction(name));
}
