#include "simulation.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(simulation) {
  py::module m("simulation", "pybind11 example plugin");

  py::class_<Simulation>(m, "Simulation")
      .def(py::init<>())
      .def_property(
          "stop_time", (void (Simulation::*)(double)) & Simulation::stop_time,
          (double (Simulation::*)()) & Simulation::stop_time, "get time");

  return m.ptr();
}