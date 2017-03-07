#! /usr/bin/env python3

import random
import argparse
import yaml

from pysinthe.feature import Polymerase, Terminator, Promoter, Mask
from pysinthe.polymer import Genome
from pysinthe.eventsignal import Signal
from pysinthe.simulation import *


def main(my_params_file=""):
    """
    TODO: REFACTOR AND VALIDATE INPUT PARAMETERS
    """

    if my_params_file == "":
        parser = argparse.ArgumentParser(description="Simulate transcription \
            and translation.")

        # Add parameter file argument
        parser.add_argument("params", metavar="p", type=str, nargs=1,
                            help="parameter file")

        args = parser.parse_args()
        my_params_file = args.params[0]

    with open(my_params_file, "r") as my_file:
        params = yaml.safe_load(my_file)

    # Set seed
    random.seed(params["simulation"]["seed"])

    # Build simulation
    simulation = Simulation()

    if "reactions" in params.keys():
        for reaction in params["reactions"]:
            new_reaction = SpeciesReaction(simulation,
                                           reaction["propensity"],
                                           reaction["reactants"],
                                           reaction["products"])
            simulation.register_reaction(new_reaction)

    # Construct list of DNA elements and a transcript template that includes
    # RBS's and stop sites
    dna_elements = []
    transcript_template = []
    position = 0
    for element in params["elements"]:
        if element["type"] == "promoter":
            length = element["length"]
            if element["length"] < 10:
                # if promoter is smaller than RNApol footprint, expand out
                # promoter length to avoid polymerases getting stuck
                length = 10
            new_element = Promoter(element["name"],
                                   position,
                                   position + length,
                                   element["interactions"].keys())
        elif element["type"] == "terminator":
            new_element = Terminator(element["name"],
                                     position + element["length"] - 1,
                                     position + element["length"],
                                     element["interactions"])
        elif element["type"] == "transcript":
            transcript_template.append(element)
            position -= element["rbs"]
            new_element = False
        element["start"] = position
        element["stop"] = position + element["length"]
        if new_element != False:
            dna_elements.append(new_element)
        position += element["length"]

    # Build genome
    if "entered" in params["genome"]:
        genome_mask = Mask("mask", params["genome"]["entered"], position, ["rnapol", "ecolipol"])
    else:
        genome_mask = Mask("mask", position, position, [])
    genome = Genome(params["genome"]["name"],
                    position, dna_elements, transcript_template, genome_mask)


    # Add species-level polymerase counts and construct list of partners
    # that this polymerase interacts with (promoters, terminators, etc.)
    for pol in params["polymerases"]:
        # add polymerase as a reactant
        simulation.increment_reactant(pol["name"], pol["copy_number"])

    # Add binding reaction for each promoter-polymerase interaction pair
    for element in params["elements"]:
        if element["type"] == "promoter":
            simulation.increment_reactant(element["name"], 0)
            for partner, constant in element["interactions"].items():
                binding_constant = constant["binding_constant"]
                for pol in params["polymerases"]:
                    if pol["name"] == partner:
                        pol_args = [partner,
                                    element["start"],
                                    10,                 # footprint
                                    pol["speed"]
                                   ]
                        reaction = Bind(simulation,
                                        float(binding_constant),
                                        element["name"],
                                        pol_args)
                        simulation.register_reaction(reaction)
    # Add species level reactants for elements that are not masked
    for element in dna_elements:
        if element.type == "promoter" and element.stop < genome.mask.start:
            simulation.increment_reactant(element.name, 1)
        elif element.type == "promoter":
            simulation.increment_reactant(element.name, 0)

    # Register genome
    genome.termination_signal.connect(simulation.terminate_transcription)
    genome.promoter_signal.connect(simulation.free_promoter)
    genome.block_signal.connect(simulation.block_promoter)
    genome.transcript_signal.connect(simulation.register_transcript)
    simulation.register_polymer(genome)

    simulation.increment_reactant("rbs", 0)
    ribo_args = ["ribosome", 0, 10, #footprint
                 30]
    # Transcript-ribosome binding reaction
    reaction = Bind(simulation, float(1e7),
                    "rbs",
                    ribo_args)
    simulation.register_reaction(reaction)

    # Add species-level ribosomes
    simulation.increment_reactant("ribosome", params["simulation"]["ribosomes"])

    simulation.initialize_propensity()
    time_step = params["simulation"]["time_step"]
    old_time = 0
    # Print initial conditions
    print(simulation)
    while simulation.time < params["simulation"]["runtime"]:
        # Execute simulatio
        simulation.execute()
        if abs(simulation.time - old_time) > time_step:
            # Output data every ~time_step
            print(simulation)
            old_time += time_step
        elif params["simulation"]["debug"] is True:
            print(simulation)
            for pol in simulation.reactions:
                print(pol)

    # print("There are ", str(len(simulation.reactions)), "reactions.")
    # for pol in simulation.reactions:
    #     print(pol)

if __name__ == "__main__":
    main()
