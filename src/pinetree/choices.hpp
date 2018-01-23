#ifndef SRC_CHOICES_HPP_ // header guard
#define SRC_CHOICES_HPP_

#include <algorithm>
#include <numeric>
#include <random>

namespace Random {
static bool seeded_ = false;
static std::mt19937 gen_;
static std::uniform_real_distribution<> dis_(0, 1);
void seed(int seed);
double random();
template <typename T>
int WeightedChoiceIndex(const std::vector<T> &population,
                        const std::vector<double> &weights) {
  double random_num = random();
  // Calculate cumulative sums for weights
  std::vector<double> cum_weights(weights.size());
  std::partial_sum(weights.begin(), weights.end(), cum_weights.begin());
  // Bisect cumulative weights vector
  std::vector<double>::iterator upper;
  upper = std::upper_bound(cum_weights.begin(), cum_weights.end(),
                           random_num * cum_weights.back());
  // Calculate index and return list item at index
  int index = (upper - cum_weights.begin());
  return index;
}
template <typename T>
T WeightedChoice(const std::vector<T> &population,
                 const std::vector<double> &weights) {
  int index = WeightedChoiceIndex(population, weights);
  return population[index];
}

template <typename T> T WeightedChoice(const std::vector<T> &population) {
  // Grab random number from python's random module
  double random_num = random();
  int index = random_num * population.size();
  return population[index];
}
}

#endif // SRC_CHOICES_HPP_