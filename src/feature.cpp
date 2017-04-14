/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "feature.hpp"

Feature::Feature(const std::string &name, int start, int stop,
  const std::vector<std::string> &interactions)
    : name_(name),
      start_(start),
      stop_(stop),
      interactions_(interactions) {
}

bool Feature::check_interaction(const std::string &name) {
  return std::find(interactions_.begin(), interactions_.end(), name) !=
    interactions_.end();
}

Polymerase::Polymerase(const std::string &name, int footprint, int speed)
    : name_(name),
      footprint_(footprint),
      speed_(speed) {
  start_ = 0;
  stop_ = footprint - 1;
  left_most_element_ = 0;
  bound_ = 0;
  type_ = "polymerase";
  reading_frame_ = 0;
}

void Polymerase::move() {
  start_++;
  stop_++;
}

void Polymerase::move_back() {
  start_--;
  stop_--;
}

Mask::Mask(const std::string &name, int start, int stop,
  const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions) {
}

Element::Element(const std::string &name, int start, int stop,
  const std::vector<std::string> &interactions)
    : Feature(name, start, stop, interactions),
      covered_(0),
      old_covered_(0) {
}
