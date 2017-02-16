#! /usr/bin/env python3

import random
import math
import argparse
import yaml

from feature import Polymerase, Terminator, Promoter, Mask
from polymer import Genome

class Reaction():
    """
    Generic class for species-level reaction. (Not currently used).
    """
    def __init__(self):
        pass

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction
        """
        pass

    def execute(self):
        """
        Execute the reaction.
        """
        pass

class SpeciesReaction(Reaction):
    """
    Generic class for species-level reaction.
    """
    def __init__(self, sim, rate_constant, reactants, products):
        super().__init__()
        self.sim = sim
        self.rate_constant = rate_constant
        self.reactants = reactants
        self.products = products

        for reactant in self.reactants:
            self.sim.increment_reactant(reactant, 0)

        for product in self.products:
            self.sim.increment_reactant(product, 0)


    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction
        """
        propensity = self.rate_constant
        for reactant in self.reactants:
            propensity *= self.sim.reactants[reactant]

        return propensity

    def execute(self):
        """
        Execute the reaction.
        """
        for reactant in self.reactants:
            self.sim.increment_reactant(reactant, -1)

        for product in self.products:
            self.sim.increment_reactant(product, 1)

class Bind(Reaction):
    """
    Bind a polymerase to a polymer.
    """

    def __init__(self, sim, polymer, rate_constant, promoter, product_args):
        """
        :param sim: reference to simulation object in which this reaction occurs
        :param polymer: polymer object involved in this reaction
        :param rate_constant: rate constant of reaction
        :param promoter: name of promoter involved in this reaction
        :param product_args: list of arguments to pass to polymerase
            constructor upon execution of this reaction.
        """
        super().__init__()
        self.sim = sim
        self.polymer = polymer
        self.polymerase_args = product_args
        self.rate_constant = rate_constant
        self.promoter = promoter # name of promoter

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction.

        TODO: Propensities could be cached to increase performance, because
        there will be many reactions in the simulation that have the exact same
        reactants.

        :returns: propensity of this reaction.
        """

        propensity = self.polymer.count_uncovered(self.promoter)

        if propensity == 0:
            return 0

        propensity = propensity * self.rate_constant * \
            self.sim.reactants[self.polymerase_args[0]]

        return propensity

    def execute(self):
        """
        Execute this reaction.
        """
        self.sim.increment_reactant(self.promoter, -1)
        self.sim.increment_reactant(self.polymerase_args[0], -1)
        # Construct and bind a new polymerase
        pol = Polymerase(*self.polymerase_args)
        self.polymer.bind_polymerase(pol, self.promoter)

class Bridge(Reaction):
    """
    Encapsulate polymer so it can participate in species-level reaction
    processing.
    """
    def __init__(self, polymer):
        """
        :param polymer: polymer object that this reaction is encapsulating
        """
        super().__init__()
        self.polymer = polymer

    def calculate_propensity(self):
        """
        Retrieve the total propensity of all reactions that may occur within
        this polymer.

        :returns: total propensity of reactions within polymer
        """
        return self.polymer.calculate_propensity()

    def execute(self):
        """
        Execute reaction within polymer (e.g. typically moving a polymerase).
        """
        self.polymer.execute()

    def __str__(self):
        """
        Retrieve string representation of polymer.

        :returns: string representation of polymer encapsulated by this object
        """
        return str(self.polymer)


class Simulation:
    """
    Coordinate polymers and species-level reactions.
    """
    def __init__(self):
        self.time = 0 # simulation time
        self.terminations = {} # track when polymerases terminate
        self.reactants = {} # species-level reactant counts
        self.reactions = [] # all reactions
        self.iteration = 0 # iteration counter

    def increment_reactant(self, name, copy_number):
        """
        Increment (or decrement) the copy number of a species-level reactant by
        copy_number.

        TODO: replace with some counter class?

        :param name: name of reactant
        :param copy_number: change in copy number (can be negative)
        """
        if name in self.reactants.keys():
            self.reactants[name] += copy_number
        else:
            self.reactants[name] = copy_number

    def register_reaction(self, reaction):
        """
        Add a SpeciesReaction object to the list of reactions.

        :param reaction: SpeciesReaction object
        """
        self.reactions.append(reaction)

    def register_polymer(self, polymer):
        """
        Add a polymer to the simulation.

        :param polymer: polymer object
        """
        # Encapsulate polymer in Bridge reaction and add to reaction list
        self.reactions.append(Bridge(polymer))

    def execute(self):
        """
        Execute one cycle of reaction.
        """
        # Generate random number
        random_num = random.random()

        # Sum propensities
        alpha_list = []
        for reaction in self.reactions:
            prop = reaction.calculate_propensity()
            alpha_list.append(prop)
        alpha = sum(alpha_list)

        # Calculate tau, i.e. time until next reaction
        tau = (1/alpha)*math.log(1/random_num)
        self.time += tau

        # Randomly select next reaction to execute, weighted by propensities
        next_reaction = random.choices(self.reactions, weights=alpha_list)[0]
        next_reaction.execute()

        self.iteration += 1

    def free_promoter(self, species):
        """
        Increment promoter count at species level.
        """
        self.increment_reactant(species, 1)

    def block_promoter(self, species):
        """
        Decrement promoter count at species level.
        """
        self.increment_reactant(species, -1)

    def register_transcript(self, polymer, reactants):
        """
        Register a new transcript with the simulation.
        """
        self.register_polymer(polymer)
        # Connect signals
        polymer.promoter_signal.connect(self.free_promoter)
        polymer.termination_signal.connect(self.terminate_translation)
        # Add species level reactants (i.e. proteins to be produced)
        for reactant in reactants:
            self.increment_reactant(reactant, 0)
        # Construct binding reactions
        for element in polymer.elements:
            if element.name == "rbs":
                # Template for ribosome to be constructed on transcript upon
                # binding.
                ribo_args = ["ribosome", element.start, 10, #footprint
                             10,
                             ["ribosome", "tstop", element.name]]
                # Transcript-ribosome binding reaction
                reaction = Bind(self, polymer, 0.05,
                                element.name,
                                ribo_args)
                self.register_reaction(reaction)

    def terminate_transcription(self, species):
        """
        Terminate transcription.

        :param polymer: the newly constructed transcript
        :param species: name of the polymerase completing transcription
        :param reactants: list of reactants to add to species-level pool
            (usually RBSs)
        """
        self.increment_reactant(species, 1)
        self.count_termination("transcript")

    def terminate_translation(self, protein, species):
        """
        Terminate translation.

        :param name: name of the protein being produced
        :param species: name of species that just translated this protein
            (usually 'ribosome')
        """
        self.increment_reactant(species, 1)
        self.increment_reactant(protein, 1)
        self.count_termination(protein)

    def count_termination(self, name):
        """
        Record the time at which a polymerase reaches a terminator.
        """
        name = name + "_total"
        if name not in self.terminations.keys():
            self.terminations[name] = 1
        else:
            self.terminations[name] += 1

    def __str__(self):
        """
        Convert `Simulation` to string representation showing names and
        termination times of `Polymerase` objects.
        """
        out_string = ""
        for name, count in self.terminations.items():
            out_string += str(self.iteration) + ", " + str(self.time) + ", " + \
                name + ", " + str(count) + "\n"
        for name, count in self.reactants.items():
            out_string += str(self.iteration) + ", " + str(self.time) + ", " + \
                name + ", " + str(count) + "\n"
        return out_string.strip()

def main():
    """
    TODO: REFACTOR AND VALIDATE INPUT PARAMETERS
    """
    parser = argparse.ArgumentParser(description="Simulate transcription \
        and translation.")

    # Add parameter file argument
    parser.add_argument("params", metavar="p", type=str, nargs=1,
                        help="parameter file")

    args = parser.parse_args()

    with open(args.params[0], "r") as my_file:
        params = yaml.safe_load(my_file)

    # Set seed
    random.seed(params["simulation"]["seed"])

    # Build simulation
    simulation = Simulation()

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
            new_element = Promoter(element["name"],
                                   position,
                                   position + element["length"],
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
    genome_mask = Mask("mask", 30, position, ["rnapol"])
    genome = Genome(params["genome"]["name"],
                    position, dna_elements, transcript_template, genome_mask)


    # Add species-level polymerase counts and construct list of partners
    # that this polymerase interacts with (promoters, terminators, etc.)
    for pol in params["polymerases"]:
        # add polymerase as a reactant
        simulation.increment_reactant(pol["name"], pol["copy_number"])
        # build interactions
        pol["interactions"] = [pol["name"]] # all pols interact with themselves
        for element in params["elements"]:
            if "interactions" in element.keys():
                if pol["name"] in element["interactions"].keys():
                    pol["interactions"].append(element["name"])

    # Add binding reaction for each promoter-polymerase interaction pair
    for element in params["elements"]:
        if element["type"] == "promoter":
            # simulation.increment_reactant(element["name"], 1)
            for partner, constant in element["interactions"].items():
                binding_constant = constant["binding_constant"]
                for pol in params["polymerases"]:
                    if pol["name"] == partner:
                        pol_args = [partner,
                                    element["start"],
                                    10,                 # footprint
                                    pol["speed"],
                                    pol["interactions"]
                                   ]
                        reaction = Bind(simulation,
                                        genome,
                                        binding_constant,
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

    # Add species-level ribosomes
    simulation.increment_reactant("ribosome", params["simulation"]["ribosomes"])

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
            # print(genome)
            for pol in simulation.reactions:
                print(pol)


if __name__ == "__main__":
    main()
