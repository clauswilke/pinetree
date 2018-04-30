#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "choices.hpp"
#include "feature.hpp"
#include "polymer.hpp"
#include "reaction.hpp"
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

  py::class_<BindingSite, std::shared_ptr<BindingSite>>(m, "BindingSite",
                                                        R"doc(
            BindingSite class that corresponds to both promoters and ribosome 
            binding sites. For internal use only.

            )doc")
      .def(py::init<std::string, int, int, std::map<std::string, double>>())
      .def("reset_state", &BindingSite::ResetState)
      .def("was_uncovered", &BindingSite::WasUncovered)
      .def("was_covered", &BindingSite::WasCovered)
      .def("cover", &BindingSite::Cover)
      .def("uncover", &BindingSite::Uncover)
      .def("is_covered", &BindingSite::IsCovered)
      .def("clone", &BindingSite::Clone)
      .def("check_interaction", &BindingSite::CheckInteraction)
      .def_property(
          "first_exposure",
          (bool (BindingSite::*)(void) const) & BindingSite::first_exposure,
          (void (BindingSite::*)(bool)) & BindingSite::first_exposure);

  py::class_<ReleaseSite, std::shared_ptr<ReleaseSite>>(m, "ReleaseSite",
                                                        R"doc(
            ReleaseSite class that corresponds to both terminators and stop codons. For internal use only.

            )doc")
      .def(py::init<std::string, int, int, std::map<std::string, double>>())
      .def("reset_state", &ReleaseSite::ResetState)
      .def("was_uncovered", &ReleaseSite::WasUncovered)
      .def("was_covered", &ReleaseSite::WasCovered)
      .def("cover", &ReleaseSite::Cover)
      .def("uncover", &ReleaseSite::Uncover)
      .def("is_covered", &ReleaseSite::IsCovered)
      .def("clone", &ReleaseSite::Clone)
      .def("check_interaction", &ReleaseSite::CheckInteraction)
      .def_property(
          "readthrough",
          (bool (ReleaseSite::*)(void) const) & ReleaseSite::readthrough,
          (void (ReleaseSite::*)(bool)) & ReleaseSite::readthrough)
      .def("efficiency", &ReleaseSite::efficiency);

  py::class_<Polymerase, std::shared_ptr<Polymerase>>(m, "Polymerase",
                                                      R"doc(
            Polymerase class that corresponds to both biological polymerases and ribosomes. For internal use only.

            )doc")
      .def(py::init<std::string, int, int>())
      .def("move", &Polymerase::Move)
      .def("move_back", &Polymerase::MoveBack)
      .def_property("start",
                    (int (Polymerase::*)(void) const) & Polymerase::start,
                    (void (Polymerase::*)(int)) & Polymerase::start)
      .def_property("stop",
                    (int (Polymerase::*)(void) const) & Polymerase::stop,
                    (void (Polymerase::*)(int)) & Polymerase::stop)
      .def_property_readonly(
          "speed", (int (Polymerase::*)(void) const) & Polymerase::speed)
      .def_property_readonly("footprint", (int (Polymerase::*)(void) const) &
                                              Polymerase::footprint)
      .def_property(
          "reading_frame",
          (int (Polymerase::*)(void) const) & Polymerase::reading_frame,
          (void (Polymerase::*)(int)) & Polymerase::reading_frame);

  py::class_<Mask, std::shared_ptr<Mask>>(m, "Mask",
                                          R"doc(
            Mask class that corresponds to polymers that are still undergoing 
            synthesis. For internal use only.

            )doc")
      .def(py::init<int, int, std::map<std::string, double>>())
      .def("move", &Mask::Move)
      .def("move_back", &Mask::MoveBack)
      .def_property("start", (int (Mask::*)(void) const) & Mask::start,
                    (void (Mask::*)(int)) & Mask::start)
      .def_property("stop", (int (Mask::*)(void) const) & Mask::stop,
                    (void (Mask::*)(int)) & Mask::stop)
      .def_property_readonly("speed", (int (Mask::*)(void) const) & Mask::speed)
      .def_property_readonly("footprint",
                             (int (Mask::*)(void) const) & Mask::footprint)
      .def_property("reading_frame",
                    (int (Mask::*)(void) const) & Mask::reading_frame,
                    (void (Mask::*)(int)) & Mask::reading_frame)
      .def("check_interaction", &Mask::CheckInteraction);

  py::class_<Rnase, std::shared_ptr<Rnase>>(m, "Rnase",
                                            R"doc(
            Rnase class that corresponds to polymers that are being degraded 
            from 5' to 3' end. For internal use only.

            )doc")
      .def(py::init<int, int>())
      .def("move", &Rnase::Move)
      .def("move_back", &Rnase::MoveBack)
      .def_property("start", (int (Rnase::*)(void) const) & Rnase::start,
                    (void (Rnase::*)(int)) & Rnase::start)
      .def_property("stop", (int (Rnase::*)(void) const) & Rnase::stop,
                    (void (Rnase::*)(int)) & Rnase::stop)
      .def_property_readonly("speed",
                             (int (Rnase::*)(void) const) & Rnase::speed)
      .def_property_readonly("footprint",
                             (int (Rnase::*)(void) const) & Rnase::footprint)
      .def_property("reading_frame",
                    (int (Rnase::*)(void) const) & Rnase::reading_frame,
                    (void (Rnase::*)(int)) & Rnase::reading_frame);

  py::class_<SpeciesReaction, std::shared_ptr<SpeciesReaction>>(
      m, "SpeciesReaction",
      R"doc(
            
            Defines reactions between two or fewer species (with stoichiometries
             of 1). For internal use only.
            
            )doc")
      .def(py::init<double, double, std::vector<std::string>,
                    std::vector<std::string>>())
      .def("caculate_propensity", &SpeciesReaction::CalculatePropensity)
      .def("execute", &SpeciesReaction::Execute)
      .def_property_readonly(
          "reactants",
          (std::vector<std::string>(SpeciesReaction::*)(void) const) &
              SpeciesReaction::reactants)
      .def_property_readonly(
          "products",
          (std::vector<std::string>(SpeciesReaction::*)(void) const) &
              SpeciesReaction::products);

  py::class_<MobileElementManager, std::shared_ptr<MobileElementManager>>(
      m, "MobileElementManager",
      R"doc(

      Manages MobileElements (polymerases, ribosomes, RNases) on a Polymer. For
      internal use only.
  
      )doc")
      .def(py::init<std::vector<double>>())
      .def("insert", &MobileElementManager::Insert)
      .def("delete", &MobileElementManager::Delete)
      .def("choose", &MobileElementManager::Choose)
      .def("valid_index", &MobileElementManager::ValidIndex)
      .def("get_pol", &MobileElementManager::GetPol)
      .def("get_attached", &MobileElementManager::GetAttached)
      .def("update_propensity", &MobileElementManager::UpdatePropensity)
      .def_property_readonly("prop_sum",
                             (double (MobileElementManager::*)(void)) &
                                 MobileElementManager::prop_sum)
      .def_property_readonly("pol_count",
                             (int (MobileElementManager::*)(void)) &
                                 MobileElementManager::pol_count);

  py::class_<Simulation, std::shared_ptr<Simulation>>(m, "Simulation",
                                                      R"doc(
            
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
           "rbs_start"_a, "rbs_stop"_a, "rbs_strength"_a)
      .def("add_rnase_site", &Genome::AddRnaseSite, "start"_a, "stop"_a);
}