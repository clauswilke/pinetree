#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include "tracker.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(pinetree, m) {
  m.doc() = (R"doc(
    Python module
    -----------------------
    .. currentmodule:: pysinthe.core
  )doc");

  py::class_<Simulation, std::shared_ptr<Simulation>>(
      m, "Simulation", "Set up and run a gene expression simulation.")
      .def(py::init<double>(), "cell_volume"_a, R"doc(
             Define a new gene expression simulation.
             
             1. **cell_volume** The volume, in liters, of the system being 
                simulated.

           )doc")
      .def("seed", &Simulation::seed,
           "set a seed for the random number generator")
      .def("add_reaction", &Simulation::AddReaction,
           "add a species-level reaction")
      .def("register_genome", &Simulation::RegisterGenome, "register a genome")
      .def("add_species", &Simulation::AddSpecies, "add species")
      .def("add_polymerase", &Simulation::AddPolymerase, "name"_a,
           "footprint"_a, "speed"_a, "copy_number"_a, "add a polymerase")
      .def("run", &Simulation::Run, "stop_time"_a, "time_step"_a,
           "output_name"_a,
           R"doc(
            
            run the simulation

            1. **stop_time** Simulated time, in seconds at which this 
                simulation should 
                stop executing reactions. Note that this *simulated* time 
                and not real time. The real time that it takes for the 
                simulation to complete depends on the number of reactions and species (genomes, transcripts, proteins, etc) in the system.
             
             2. **time_step** Time interval, in seconds, that species counts 
                are reported.
            
             3. **output_name** Prefix for output files.

          )doc");

  // Polymers, genomes, and transcripts
  py::class_<Polymer, Polymer::Ptr>(m, "Polymer");
  py::class_<Genome, Polymer, Genome::Ptr>(m, "Genome")
      .def(py::init<const std::string &, int>(), "name"_a, "length"_a)
      .def("add_mask", &Genome::AddMask, "start"_a, "interactions"_a)
      .def("add_weights", &Genome::AddWeights, "weights"_a)
      .def("add_promoter", &Genome::AddPromoter, "name"_a, "start"_a, "stop"_a,
           "interactions"_a)
      .def("add_terminator", &Genome::AddTerminator, "name"_a, "start"_a,
           "stop"_a, "efficiency"_a)
      .def("add_gene", &Genome::AddGene, "name"_a, "start"_a, "stop"_a,
           "rbs_start"_a, "rbs_stop"_a, "rbs_strength"_a);
}