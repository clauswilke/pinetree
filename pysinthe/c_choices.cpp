<%
setup_pybind11(cfg)
%>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <algorithm>

namespace py = pybind11;

py::object weighted_choice(py::list population,
                           const std::vector<double> &weights) {
    // Grab random number from python's random module
    py::object randfunc = py::module::import("random").attr("random");
    double random_num = randfunc().cast<double>();
    // Calculate cumulative sums for weights
    double sum = 0;
    std::vector<double> cum_weights = std::vector<double>(weights.size());
    for (std::vector<double>::size_type i = 0; i < cum_weights.size(); i++) {
      sum = weights[i] + sum;
      cum_weights[i] = sum;
    }
    // Bisect cumulative weights vector
    std::vector<double>::iterator upper;
    upper = std::upper_bound(cum_weights.begin(),
                             cum_weights.end(),
                             random_num * sum);
    // Calculate index and return list item at index
    int index = (upper - cum_weights.begin());
    return population[index];
}

PYBIND11_PLUGIN(c_choices) {
    py::module m("c_choices", "weighted_choice rewritten in C++");
    m.def("weighted_choice", &weighted_choice);
    return m.ptr();
}
