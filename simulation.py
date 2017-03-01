#! /usr/bin/env python3

import random
import math
import argparse
import yaml

from pysinthe.feature import Polymerase, Terminator, Promoter, Mask
from pysinthe.polymer import Genome
from pysinthe.eventsignal import Signal

AVAGADRO = float(6.0221409e+23)
# CELL_VOLUME = float(8e-16)
CELL_VOLUME = float(8e-16)

class Reaction():
    """
    Generic class for a reaction. (Not currently used).
    """
    def __init__(self):
        self.index = -1 # Index inside propensity list, sim.alpha_list

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
        if len(reactants) == 2:
            self.rate_constant = float(rate_constant)/(AVAGADRO*CELL_VOLUME)
        else:
            self.rate_constant = float(rate_constant)
        self.reactants = reactants
        self.products = products

        for reactant in self.reactants:
            self.sim.increment_reactant(reactant, 0)
            if reactant not in self.sim.reactant_bind_map:
                self.sim.reactant_bind_map[reactant] = [self]
            else:
                self.sim.reactant_bind_map[reactant].append(self)

        for product in self.products:
            self.sim.increment_reactant(product, 0)
            if product not in self.sim.reactant_bind_map:
                self.sim.reactant_bind_map[product] = [self]
            else:
                self.sim.reactant_bind_map[product].append(self)


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
    def __init__(self, sim, rate_constant, promoter, pol_args):
        """
        :param sim: reference to simulation object in which this reaction occurs
        :param rate_constant: rate constant of reaction
        :param promoter: name of promoter involved in this reaction
        :param pol_args: list of arguments to pass to polymerase
            constructor upon execution of this reaction.
        """
        super().__init__()
        self.sim = sim
        self.polymerase = pol_args[0]
        self.pol_args = pol_args
        self.rate_constant = float(rate_constant)/(AVAGADRO*CELL_VOLUME)
        self.promoter = promoter
        if self.polymerase not in self.sim.reactant_bind_map:
            self.sim.reactant_bind_map[self.polymerase] = [self]
        else:
            self.sim.reactant_bind_map[self.polymerase].append(self)
        if self.promoter not in self.sim.reactant_bind_map:
            self.sim.reactant_bind_map[self.promoter] = [self]
        else:
            self.sim.reactant_bind_map[self.promoter].append(self)

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction.

        TODO: Propensities could be cached to increase performance, because
        there will be many reactions in the simulation that have the exact same
        reactants.

        :returns: propensity of this reaction.
        """

        return self.rate_constant * self.sim.reactants[self.polymerase]\
            * self.sim.reactants[self.promoter]

    def execute(self):
        """
        Decrement reactants, choose a polymer, construct new polymerase, and
        bind the polymerase to the polymer.
        """
        # Find which polymer we should bind to
        weights = []
        for polymer in self.sim.promoter_polymer_map[self.promoter]:
            # Weight by the number of unbound promoters in each polymer
            weights.append(polymer.count_uncovered(self.promoter))
        polymer = random.choices(self.sim.promoter_polymer_map[self.promoter],
                                 weights=weights)[0]
        # Construct new polymerase
        new_pol = Polymerase(*self.pol_args)
        polymer.bind_polymerase(new_pol, self.promoter)

        self.sim.increment_reactant(self.promoter, -1)
        self.sim.increment_reactant(self.polymerase, -1)

        # self.sim.update_propensity(self.index)

    def __str__(self):
        return self.promoter + "-" + self.polymerase

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
        self.propensity_signal = Signal() # Signal to fire when propensity needs
                                          # to be updated
        polymer.propensity_signal.connect(self.update)

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

    def update(self):
        """
        Fires when polymer removes a polymerase from polymer and propensity
        needs to be updated.
        """
        self.propensity_signal.fire(self.index)

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
        self.promoter_polymer_map = {} # Map of which polymers contain a given
                                       # promoter
        self.reactant_bind_map = {} # Map of which binding reaction involves a
                                    # given reactant
        self.reactions = [] # all reactions
        self.iteration = 0 # iteration counter
        self.alpha_list = []

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

        if name in self.reactant_bind_map:
            for reaction in self.reactant_bind_map[name]:
                self.update_propensity(reaction.index)

    def register_reaction(self, reaction):
        """
        Add a SpeciesReaction object to the list of reactions.

        :param reaction: SpeciesReaction object
        """
        if reaction not in self.reactions:
            reaction.index = len(self.reactions)
            try:
                self.alpha_list.append(reaction.calculate_propensity())
            except KeyError:
                self.alpha_list.append(0.0)
            self.reactions.append(reaction)

    def register_polymer(self, polymer):
        """
        Add a polymer to the simulation.

        :param polymer: polymer object
        """
        # Add polymer to promoter-polymerase map
        for element in polymer.elements:
            if element.type == "promoter":
                try:
                    self.promoter_polymer_map[element.name].append(polymer)
                except KeyError:
                    self.promoter_polymer_map[element.name] = [polymer]

        # Encapsulate polymer in Bridge reaction and add to reaction list
        bridge = Bridge(polymer)
        bridge.propensity_signal.connect(self.update_propensity)
        self.register_reaction(bridge)

    def initialize_propensity(self):
        """
        Initialize all propensities before the start of the simulation.
        """
        for index, reaction in enumerate(self.reactions):
            self.alpha_list[index] = reaction.calculate_propensity()

    def update_propensity(self, index):
        """
        Update a propensity of a reaction at a given index.
        """
        self.alpha_list[index] = float(self.reactions[index].calculate_propensity())

    def execute(self):
        """
        Execute one cycle of reaction.
        """
        # Generate random number
        random_num = random.random()

        # Sum propensities
        alpha = sum(self.alpha_list)
        # Calculate tau, i.e. time until next reaction
        tau = (1/alpha)*math.log(1/random_num)
        self.time += tau

        # Randomly select next reaction to execute, weighted by propensities
        # print(self.reactions)
        next_reaction = random.choices(self.reactions,
                                       weights=self.alpha_list)[0]
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

    def register_transcript(self, polymer):
        """
        Register a new transcript with the simulation.

        :param polymer: transcript object to be added to simulation
        :param reactants: proteins to be produced by this transcript
        """
        self.register_polymer(polymer)
        # Connect signals
        polymer.promoter_signal.connect(self.free_promoter)
        polymer.termination_signal.connect(self.terminate_translation)

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
        Record when a polymerase reaches a terminator.
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
            out_string += str(self.iteration) + ", " + str(float(self.time)) + ", " + \
                name + ", " + str(count) + "\n"
        for name, count in self.reactants.items():
            out_string += str(self.iteration) + ", " + str(float(self.time)) + ", " + \
                name + ", " + str(count) + "\n"
        return out_string.strip()

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
    if "mask" in params["genome"]:
        genome_mask = Mask("mask", 30, position, ["rnapol", "ecolipol"])
    else:
        genome_mask = Mask("mask", position, position, [])
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
            simulation.increment_reactant(element["name"], 0)
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
                 30,
                 ["ribosome", "tstop", "rbs"]]
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
