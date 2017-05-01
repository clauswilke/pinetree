#include "choices.hpp"
#include "simulation.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_PLUGIN(pysinthe) {
  py::module m("pysinthe", "pybind11 example plugin");

  m.def("seed", &Random::seed, "set a global seed for the simulation");

  py::class_<Simulation>(m, "Simulation")
      .def(py::init<>())
      .def_property("stop_time",
                    (void (Simulation::*)(double)) & Simulation::stop_time,
                    (double (Simulation::*)()) & Simulation::stop_time,
                    "stop time of simulation")
      .def_property("time_step",
                    (void (Simulation::*)(double)) & Simulation::time_step,
                    (double (Simulation::*)()) & Simulation::time_step,
                    "time step at which to output data");

  return m.ptr();
}