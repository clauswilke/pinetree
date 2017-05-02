
/* <% setup_pybind11(cfg) %> */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <numeric>

namespace py = pybind11;

py::object weighted_choice(py::list population,
                           const std::vector<double> &weights) {
    // Grab random number from python's random module
    py::object randfunc = py::module::import("random").attr("random");
    double random_num = randfunc().cast<double>();
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
    return population[index];
}

py::object choice(py::list population) {
    // Grab random number from python's random module
    py::object randfunc = py::module::import("random").attr("random");
    double random_num = randfunc().cast<double>();
    int index = random_num * population.size();
    return population[index];
}

PYBIND11_PLUGIN(c_choices) {
    py::module m("c_choices", "weighted_choice rewritten in C++");
    m.def("weighted_choice", &weighted_choice,
          py::arg("population").noconvert(), py::arg("weights"));
    m.def("weighted_choice", &choice, py::arg("population").noconvert());
    return m.ptr();
}
