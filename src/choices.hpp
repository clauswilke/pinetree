#ifndef SRC_CHOICES_HPP_ // header guard
#define SRC_CHOICES_HPP_

#include <algorithm>
#include <numeric>
#include <random>

class Random
{
  public:
    Random();
    void seed(int seed);
    double random();
    template <class T>
    int WeightedChoice(const std::vector<T> &population,
                       const std::vector<double> &weights)
    {
        double random_num = random();
        // Calculate cumulative sums for weights
        std::vector<double> cum_weights(weights.size());
        std::partial_sum(weights.begin(), weights.end(), cum_weights.begin());
        // Bisect cumulative weights vector
        std::vector<double>::iterator upper;
        upper = std::upper_bound(cum_weights.begin(),
                                 cum_weights.end(),
                                 random_num * cum_weights.back());
        // Calculate index and return list item at index
        int index = (upper - cum_weights.begin());
        return index;
    }

  private:
    static std::mt19937 gen_;
    static std::uniform_real_distribution<> dis_;
};

// py::object Choice(py::list population) {
//     // Grab random number from python's random module
//     py::object randfunc = py::module::import("random").attr("random");
//     double random_num = randfunc().cast<double>();
//     int index = random_num * population.size();
//     return population[index];
// }

#endif // SRC_CHOICES_HPP_