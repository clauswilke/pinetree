#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_PLUGIN(core) {
  py::module m("core", "pybind11 example plugin");

  m.def("seed", &Random::seed, "set a global seed for the simulation");

  py::class_<SpeciesTracker>(m, "SpeciesTracker")
      .def("get_instance", SpeciesTracker::Instance,
           py::return_value_policy::reference)
      .def("increment", &SpeciesTracker::Increment);

  py::class_<Simulation>(m, "Simulation")
      .def(py::init<>())
      .def_property("stop_time",
                    (double (Simulation::*)()) & Simulation::stop_time,
                    (void (Simulation::*)(double)) & Simulation::stop_time,
                    "stop time of simulation")
      .def_property("time_step",
                    (double (Simulation::*)()) & Simulation::time_step,
                    (void (Simulation::*)(double)) & Simulation::time_step,
                    "time step at which to output data")
      .def("register_reaction", &Simulation::RegisterReaction,
           "register a species-level reaction")
      .def("register_genome", &Simulation::RegisterGenome, "register a genome");

  // Binding for abtract Reaction so pybind11 doesn't complain when doing
  // conversions between Reaction and its child classes
  py::class_<Reaction, Reaction::Ptr>(m, "Reaction");
  py::class_<SpeciesReaction, SpeciesReaction::Ptr, Reaction>(m,
                                                              "SpeciesReaction")
      .def(py::init<double, const std::vector<std::string> &,
                    const std::vector<std::string> &>());
  py::class_<Bind, std::shared_ptr<Bind>, Reaction>(m, "Bind").def(
      py::init<double, const std::string &, const Polymerase &>());
  py::class_<Bridge, std::shared_ptr<Bridge>, Reaction>(m, "Bridge")
      .def(py::init<Polymer::Ptr>());

  // Features and elements
  py::class_<Feature>(m, "Feature")
      .def_property("start", (int (Feature::*)()) & Feature::start,
                    (void (Feature::*)(int)) & Feature::start)
      .def_property("stop", (int (Feature::*)()) & Feature::stop,
                    (void (Feature::*)(int)) & Feature::stop);
  py::class_<Mask, Feature>(m, "Mask").def(
      py::init<const std::string &, int, int,
               const std::vector<std::string> &>());
  py::class_<Element, Feature>(m, "Element");
  py::class_<Promoter, Element>(m, "Promoter")
      .def(py::init<const std::string &, int, int,
                    const std::vector<std::string> &>());
  py::class_<Terminator, Element>(m, "Terminator")
      .def(py::init<const std::string &, int, int,
                    const std::vector<std::string> &,
                    const std::map<std::string, double> &>());

  // Polymers, genomes, and transcripts
  py::class_<Polymer>(m, "Polymer");
  py::class_<Transcript, Polymer>(m, "Transcript")
      .def(py::init<const std::string &, int, int, const Element::VecPtr &,
                    const Mask &>());
  py::class_<Genome, Polymer>(m, "Genome")
      .def(py::init<const std::string &, int, const Element::VecPtr &,
                    const Element::VecPtr &, const Mask &>());

  return m.ptr();
}