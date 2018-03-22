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

Promoter::Promoter(const std::string &name, int start, int stop,
                   const std::map<std::string, double> &interactions)
    : FixedElement(name, start, stop, interactions) {
  first_exposure_ = false;
}

bool Promoter::CheckInteraction(const std::string &name) {
  return interactions_.count(name);
}

Promoter::Ptr Promoter::Clone() const {
  return std::make_shared<Promoter>(*this);
}

Terminator::Terminator(const std::string &name, int start, int stop,
                       const std::map<std::string, double> &interactions)
    : FixedElement(name, start, stop, interactions) {
  readthrough_ = false;
}

bool Terminator::CheckInteraction(const std::string &name, int reading_frame) {
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

Terminator::Ptr Terminator::Clone() const {
  return std::make_shared<Terminator>(*this);
}

MobileElement::MobileElement(const std::string &name, int footprint, int speed)
    : name_(name), footprint_(footprint), speed_(speed), reading_frame_(-1) {}

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : MobileElement(name, footprint, speed) {
  type_ = "__polymerase";
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

Rnase::Rnase(int footprint, int speed)
    : MobileElement("__rnase", footprint, speed) {}

Mask::Mask(int start, int stop,
           const std::map<std::string, double> &interactions)
    : MobileElement("__mask", 0, 0), interactions_(interactions) {
  start_ = start;
  stop_ = stop;
}

bool Mask::CheckInteraction(const std::string &name) {
  return interactions_.count(name);
}