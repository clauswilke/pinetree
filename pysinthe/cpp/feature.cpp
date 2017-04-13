/* Copyright (c) 2017 Benjamin Jack All Rights Reserved. */

#include <algorithm>
#include <string>
#include <vector>

#include "./feature.h"

Feature::Feature(const std::string &name, int start, int stop,
  std::vector<std::string> interactions)
    : name_(name),
      start_(start),
      stop_(stop),
      interactions_(interactions) {
}

bool Feature::check_interaction(const std::string &name) {
  return std::find(interactions_.begin(), interactions_.end(), name) !=
    interactions_.end();
}
