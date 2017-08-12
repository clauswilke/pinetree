#include "choices.hpp"

void Random::seed(int seed) {
  gen_.seed(seed);
  seeded_ = true;
}

double Random::random() {
  if (!seeded_) {
    std::random_device rd;
    gen_.seed(rd());
    seeded_ = true;
  }
  return dis_(gen_);
}