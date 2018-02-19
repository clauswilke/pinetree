#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"
#include "simulation.hpp"
#include "tracker.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(core, m) {
  m.doc() = (R"doc(
    Python module
    -----------------------
    .. currentmodule:: pinetree
  )doc");

  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation", R"doc(
            
            Define and run a gene expression simulation.
            
            Args:
                cell_volume (float): The volume, in liters, of the system being simulated.
             
            Examples:

                >>> import pinetree.pinetree as pt
                >>> sim = pt.Simulation(cell_volume=8e-16) # Approximate volume of E. coli cell

           )doc")
      .def(py::init<double>(), "cell_volume"_a)
      .def("seed", &Simulation::seed,
           R"doc(
             
             Set a seed for reproducible simulations.
             
             Args:
                seed (int): a seed for the random number generator

             )doc")
      .def("add_reaction", &Simulation::AddReaction, "rate_constant"_a,
           "reactants"_a, "products"_a, R"doc(
             
            Define a reaction between species, which may include free 
            ribosomes and polymerases.

            Args:
                rate_constant (float): Macroscopic rate constant of the 
                    reaction. This will be converted into a stochastic 
                    mesoscopic rate constant automatically.
                reactants (list): List of reactants, which may be species, 
                    ribosomes, or polymerases
                products (list): List of products which may be species, 
                    ribosomes, or polymerases
        
            Note:
                Reaction rate constants should be given as macroscopic rate 
                constants, the same constants used in differential 
                equation-based models. The simulation will automatically 
                convert these rate constants to mesoscopic constants required 
                for a stochastic simulation.
            
            Example:

                >>> coming soon
             
           )doc")
      .def("add_species", &Simulation::AddSpecies, "name"_a, "copy_number"_a,
           R"doc(
           
           Defines individual chemical species not specified by either ``add_ribosome()`` or ``add_polymerase()``.

           Args:
              name (str): Name of chemical species which can be referred to in    reactions added with ``add_reaction()``.
              copy_number (int): Initial number of copies of the chemical species 

           )doc")
      .def("add_polymerase", &Simulation::AddPolymerase, "name"_a,
           "footprint"_a, "speed"_a, "copy_number"_a, R"doc(

           Add a polymerase to the simulation. Each type of polymerase should define the following fields:
           
           Args:
              name (str): Name of the polymerase which can be referred to in 
                  ``add_reaction()`` and ``add_promoter()``.
              copy_number (int): Initial number of copies of the polymerase
              speed (int): Speed, in base pairs per second, at which the 
                  polymerase transcribes
              footprint (int): Footprint, in base pairs, of the polymerase on 
                  the genome

           )doc")
      .def("register_genome", &Simulation::RegisterGenome, R"doc(
        
        Register a genome with the simulation.

        Args:
            genome (Genome): a pinetree ``Genome`` object.
        
        )doc")
      .def("run", &Simulation::Run, "stop_time"_a, "time_step"_a,
           "output_prefix"_a, R"doc(
            
            Run the simulation.

            Args:
                stop_time (int): Simulated time, in seconds at which this 
                    simulation should stop executing reactions. Note that this 
                    *simulated* time and not real time. The real time that it 
                    takes for the simulation to complete depends on the number 
                    of reactions and species (genomes, transcripts, proteins, 
                    etc) in the system.
                time_step (int): Time interval, in seconds, that species counts 
                    are reported.
                output_prefix (str): Prefix for output files.

          )doc");

  // Polymers, genomes, and transcripts
  py::class_<Polymer, Polymer::Ptr>(m, "Polymer");
  py::class_<Genome, Polymer, Genome::Ptr>(m, "Genome")
      .def(py::init<const std::string &, int, double>(), "name"_a, "length"_a,
           "transcript_degradation_rate"_a = 0.0)
      .def("add_mask", &Genome::AddMask, "start"_a, "interactions"_a)
      .def("add_weights", &Genome::AddWeights, "weights"_a)
      .def("add_promoter", &Genome::AddPromoter, "name"_a, "start"_a, "stop"_a,
           "interactions"_a)
      .def("add_terminator", &Genome::AddTerminator, "name"_a, "start"_a,
           "stop"_a, "efficiency"_a)
      .def("add_gene", &Genome::AddGene, "name"_a, "start"_a, "stop"_a,
           "rbs_start"_a, "rbs_stop"_a, "rbs_strength"_a);
}