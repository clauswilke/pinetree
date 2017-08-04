#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include "tracker.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_PLUGIN(core) {
  py::module m("core", R"doc(
    Python module
    -----------------------
    .. currentmodule:: pysinthe.core
    .. autosummary::
        :toctree: _generate

        Simulation
        Reaction
        SpeciesReaction
  )doc");

  m.def("seed", &Random::seed, "set a global seed for the simulation");

  py::class_<SpeciesTracker>(m, "SpeciesTracker")
      .def("get_instance", SpeciesTracker::Instance,
           py::return_value_policy::reference)
      .def("increment", &SpeciesTracker::Increment);

  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation")
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
      .def("register_genome", &Simulation::RegisterGenome, "register a genome")
      .def("run", &Simulation::Run, "run the simulation");

  // Binding for abtract Reaction so pybind11 doesn't complain when doing
  // conversions between Reaction and its child classes
  py::class_<Reaction, Reaction::Ptr>(m, "Reaction");
  py::class_<SpeciesReaction, Reaction, SpeciesReaction::Ptr>(m,
                                                              "SpeciesReaction")
      .def(py::init<double, double, const std::vector<std::string> &,
                    const std::vector<std::string> &>());
  py::class_<Bind, Reaction, std::shared_ptr<Bind>>(m, "Bind").def(
      py::init<double, double, const std::string &, const Polymerase &>());
  py::class_<Bridge, Reaction, std::shared_ptr<Bridge>>(m, "Bridge")
      .def(py::init<Polymer::Ptr>());

  // Features and elements
  py::class_<Feature, std::shared_ptr<Feature>>(m, "Feature")
      .def_property("start", (int (Feature::*)()) & Feature::start,
                    (void (Feature::*)(int)) & Feature::start)
      .def_property("stop", (int (Feature::*)()) & Feature::stop,
                    (void (Feature::*)(int)) & Feature::stop)
      .def_property_readonly(
          "type", (const std::string &(Feature::*)() const) & Feature::type)
      .def_property_readonly(
          "name", (const std::string &(Feature::*)() const) & Feature::name);
  py::class_<Mask, Feature, std::shared_ptr<Mask>>(m, "Mask").def(
      py::init<const std::string &, int, int,
               const std::vector<std::string> &>());
  py::class_<Element, Feature, Element::Ptr>(m, "Element")
      .def_property("gene",
                    (const std::string &(Element::*)() const) & Element::gene,
                    (void (Element::*)(const std::string &)) & Element::gene);
  py::class_<Promoter, Element, Promoter::Ptr>(m, "Promoter")
      .def(py::init<const std::string &, int, int,
                    const std::vector<std::string> &>());
  py::class_<Terminator, Element, Terminator::Ptr>(m, "Terminator")
      .def(py::init<const std::string &, int, int,
                    const std::vector<std::string> &,
                    const std::map<std::string, double> &>())
      .def_property(
          "reading_frame", (int (Terminator::*)()) & Terminator::reading_frame,
          (void (Terminator::*)(int)) & Terminator::set_reading_frame);
  py::class_<Polymerase, Feature, Polymerase::Ptr>(m, "Polymerase")
      .def(py::init<const std::string &, int, int>());

  // Polymers, genomes, and transcripts
  py::class_<Polymer, Polymer::Ptr>(m, "Polymer");
  py::class_<Genome, Polymer, Genome::Ptr>(m, "Genome")
      .def(py::init<const std::string &, int, const Element::VecPtr &,
                    const Element::VecPtr &, const Mask &>())
      .def(py::init<const std::string &, int, const Element::VecPtr &,
                    const Element::VecPtr &, const Mask &,
                    const std::vector<double> &>());

  return m.ptr();
}