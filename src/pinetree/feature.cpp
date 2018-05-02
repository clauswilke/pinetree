/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "feature.hpp"
#include "tracker.hpp"

FixedElement::FixedElement(const std::string &name, int start, int stop,
                           const std::map<std::string, double> &interactions)
    : name_(name),
      start_(start),
      stop_(stop),
      interactions_(interactions),
      covered_(0),
      old_covered_(0),
      reading_frame_(-1) {
  if (start_ < 0 || stop_ < 0) {
    throw std::invalid_argument(
        "Fixed element '" + name_ +
        "' has a negative start and/or stop coordinate.");
  }
}

FixedElement::~FixedElement(){};

BindingSite::BindingSite(const std::string &name, int start, int stop,
                         const std::map<std::string, double> &interactions)
    : FixedElement(name, start, stop, interactions) {
  first_exposure_ = false;
  for (auto const &item : interactions) {
    if (item.second < 0) {
      throw std::invalid_argument(
          "Binding site '" + name_ +
          "' must have non-negative interaction rate constants.");
    }
  }
}

bool BindingSite::CheckInteraction(const std::string &name) {
  return interactions_.count(name);
}

BindingSite::Ptr BindingSite::Clone() const {
  return std::make_shared<BindingSite>(*this);
}

void BindingSite::Degrade() {
  if (covered_ == 0) {
    std::runtime_error(
        "Attempting to mark an uncovered binding site for degradation.");
  } else {
    degraded_ = true;
  }
}

ReleaseSite::ReleaseSite(const std::string &name, int start, int stop,
                         const std::map<std::string, double> &interactions)
    : FixedElement(name, start, stop, interactions) {
  first_exposure_ = false;
  readthrough_ = false;
  for (auto const &item : interactions) {
    if (item.second < 0 || item.second > 1) {
      throw std::invalid_argument(
          "Release site '" + name_ +
          "' must have efficiency values between 0.0 and 1.0.");
    }
  }
}

bool ReleaseSite::CheckInteraction(const std::string &name, int reading_frame) {
  if (interactions_.count(name) == 1) {
    if (reading_frame_ == -1) {
      return true;
    }
    if (reading_frame == reading_frame_) {
      return true;
    }
  }
  return false;
}

ReleaseSite::Ptr ReleaseSite::Clone() const {
  return std::make_shared<ReleaseSite>(*this);
}

MobileElement::MobileElement(const std::string &name, int footprint, int speed)
    : name_(name), footprint_(footprint), speed_(speed), reading_frame_(-1) {
  start_ = 0;
  stop_ = start_ + footprint_;
  if (footprint_ < 0) {
    throw std::invalid_argument("Mobile element '" + name_ +
                                "' has a negative footprint size.");
  }
  if (speed_ < 0) {
    throw std::invalid_argument("Mobile element '" + name_ +
                                "' has a negative average speed.");
  }
}

MobileElement::~MobileElement(){};

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : MobileElement(name, footprint, speed) {
  reading_frame_ = -1;
}

void Polymerase::Move() {
  start_++;
  stop_++;
}

void Polymerase::MoveBack() {
  if (start_ > 0) {
    start_--;
    stop_--;
  } else {
    throw std::runtime_error(
        "Attempting to assign negative start position to Polymerase object '" +
        name_ + "'.");
  }
}

Mask::Mask(int start, int stop,
           const std::map<std::string, double> &interactions)
    : MobileElement("__mask", stop - start + 1, 0),
      interactions_(interactions) {
  start_ = start;
  stop_ = stop;
}

/**
 * TODO: refactor to make vector rather than map
 */
bool Mask::CheckInteraction(const std::string &name) {
  return interactions_.count(name);
}

void Mask::MoveBack() {
  if (start_ > 0) {
    start_--;
    footprint_++;
  } else {
    throw std::runtime_error(
        "Attempting to assign negative start position to Mask object '" +
        name_ + "'.");
  }
}

Rnase::Rnase(int footprint, int speed)
    : MobileElement("__rnase", footprint, speed) {}