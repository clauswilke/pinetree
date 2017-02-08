#! /usr/bin/env python3

import heapq
import random
import argparse
import yaml

from feature import *
from polymer import *

class SpeciesReaction:
    """
    Generic class for species-level reaction. (Not currently used).
    """
    def __init__(self):
        pass

    def calculate_time(self, current_time):
        """
        Calculate the time at which this reaction will occur next.

        :param current_time: current time of simulation
        """
        pass

    def execute(self, current_time):
        """
        Execute the reaction.

        :param current_time: current time of simulation
        """
        pass

class Bind(SpeciesReaction):
    """
    Bind a polymerase to a polymer.
    """
    def __init__(self, sim, polymer, rate_constant, reactants, product_args):
        """
        :param sim: reference to simulation object in which this reaction occurs
        :param polymer: polymer object involved in this reaction
        :param rate_constant: rate constant of reaction
        :param reactants: *names* of species-level reactants involved in this
            reaction
        :param product_args: list of arguments to pass to polymerase
            constructor upon execution of this reaction.
        """
        super().__init__()
        self.sim = sim
        self.polymer = polymer
        self.polymerase_args = product_args
        self.rate_constant = rate_constant
        self.reactants = reactants

    def calculate_time(self, current_time):
        """
        Calculate the time at which this reaction will occur next.

        :param current_time: the current time of the simulation
        :returns: time at which this reaction will next occur.
        """
        propensity = self.rate_constant
        for reactant in self.reactants:
            propensity *= self.sim.reactants[reactant]
        try:
            # calculate delta time according to exponential distribution with
            # the propensity as the mean
            delta_time = random.expovariate(propensity)
        except ZeroDivisionError:
            # propensity could be 0, causing divide-by-zero error
            delta_time = float('Inf')

        return current_time + delta_time

    def execute(self, current_time):
        """
        Execute this reaction.

        :param current_time: the current time of the simulation
        """
        for reactant in self.reactants:
            # Decrement each reactant by 1
            self.sim.increment_reactant(reactant, -1)
        # Construct and bind a new polymerase
        pol = Polymerase(*self.polymerase_args)
        self.polymer.bind_polymerase(pol, current_time)

class Bridge(SpeciesReaction):
    """
    Encapsulate polymer so it can participate in species-level reaction queue.
    """
    def __init__(self, polymer):
        """
        :param polymer: polymer object that this reaction is encapsulating
        """
        self.polymer = polymer

    def calculate_time(self, current_time):
        """
        Retrieve time of next reaction within the polymer (e.g. a polymerase
        moving).

        :param current_time: current time of simulation
        :returns: time at which next reaction will occur within polymer
        """
        return self.polymer.get_next_time()

    def execute(self, current_time):
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
        self.heap = [] # heap with next reaction on top
        self.iteration = 0 # iteration counter

    def increment_reactant(self, name, copy_number):
        """
        Increment (or decrement) the copy number of a species-level reactant by
        copy_number.

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
        if reaction not in self.reactions:
            self.reactions.append(reaction)

    def register_polymer(self, polymer):
        """
        Add a polymer to the simulation.

        :param polymer: polymer object
        """
        # Register this simulation as an observer of the polymer, so the polymer
        # can send messages to the simulation about its internal state.
        polymer.register_observer(self)
        # Encapsulate polymer in Bridge reaction and add to reaction list
        self.reactions.append(Bridge(polymer))

    def build_heap(self):
        """
        Construct reaction min-heap.
        """
        # Empty out heap just in case
        self.heap = []
        # Construct heap based on minimum next-reaction times for each reaction
        # in this simulation.
        for i in range(len(self.reactions)):
            self.heap.append((self.reactions[i].calculate_time(self.time), i))
        # Construct heap
        heapq.heapify(self.heap)

    def pop(self):
        """
        Remove next-occuring reaction from heap.

        :returns: reaction time, index of reaction in self.reactions
        """
        reaction = heapq.heappop(self.heap)
        return reaction[0], reaction[1]

    def push_reaction(self, reaction_index):
        """
        Add reaction back into min-heap.

        :param reaction_index: index of reaction within self.reactions
        """
        heapq.heappush(self.heap,
            (self.reactions[reaction_index].calculate_time(self.time),
             reaction_index))

    def execute(self):
        """
        Execute one cycle of reaction.

        TODO: avoid reconstructing heap on every iteration
        """
        time, index = self.pop()
        self.time = time
        self.reactions[index].execute(self.time)
        # Reconstruct heap on each iteration
        self.build_heap()
        self.iteration += 1

    def notify(self, observable, **kwargs):
        """
        Receive and respond to messages from polymers, and polymerases.

        Messages include:
        * terminate_transcript from polymerase: construct transcript, add
            polymerase back into species-level pool, register transcript with simulation, add binding reaction for RBS, count completed transcript
        * free_promoter from promoter: add promoter back to species-level pool
        * free_promoter from rbs: add an RBS back to species-level pool
        * terminate from ribosome: add ribosome and completed protein into
            species-level pool, count completed protein

        TODO: REFACTOR, DEAL WITH FOOTPRINT SIZES

        :param observable: object delivering messages, usually a polymer?
        :param kwargs: messages (e.g. "terminate_transcript", "free_promoter",
            etc.)
        """
        if kwargs["action"] == "terminate_transcript" and kwargs["type"] == "polymerase":
            self.increment_reactant(kwargs["species"], 1)
            self.register_polymer(kwargs["polymer"])
            for reactant in kwargs["reactants"]:
                self.increment_reactant(reactant, 1)
            # Construct binding reaction
            for element in kwargs["polymer"].elements:
                if element.name == "rbs":
                    # Template for ribosome to be constructed on transcript upon
                    # binding.
                    ribo_args = ["ribosome", element.start, 10, #footprint
                                  10,
                                 ["ribosome", "tstop", "rbs"]]
                    # Transcript-ribosome binding reaction
                    reaction = Bind(self, kwargs["polymer"], 0.05,
                                    ["rbs", "ribosome"],
                                    ribo_args)
                    self.register_reaction(reaction)
            # Count that a transcript has been constructed
            self.count_termination("full_transcript", self.time)
        elif kwargs["action"] == "free_promoter" and kwargs["type"] == "promoter":
            self.increment_reactant(kwargs["species"], 1)
        elif kwargs["action"] == "free_promoter" and kwargs["type"] == "rbs":
            # free ribosome binding site
            self.increment_reactant(kwargs["species"], 1)
        elif kwargs["action"] == "terminate" and kwargs["species"] == "ribosome":
            # complete protein synthesis
            self.increment_reactant(kwargs["species"], 1)
            self.increment_reactant(kwargs["name"], 1)
            self.count_termination(kwargs["name"], self.time)

    def count_termination(self, name, time):
        """
        Record the time at which a polymerase reaches a terminator.
        """
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
            out_string += str(self.iteration) + ", " + str(self.time) + ", " + name + ", " + str(count) + "\n"
        for name, count in self.reactants.items():
            out_string += str(self.iteration) + ", " + str(self.time) + ", " + name + ", " + str(count) + "\n"
        return out_string.strip()

def main():

    parser = argparse.ArgumentParser(description = "Simulate transcription \
        and translation.")

    # Add parameter file argument
    parser.add_argument("params", metavar = "p", type = str, nargs = 1,
                        help = "parameter file")

    args = parser.parse_args()

    with open(args.params[0], "r") as f:
        params = yaml.safe_load(f)

    random.seed(params["simulation"]["seed"])

    simulation = Simulation()

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
                                     position,
                                     position + element["length"],
                                     element["interactions"].keys())
        elif element["type"] == "transcript":
            transcript_template.append(element)
            new_element = False
        element["start"] = position
        element["stop"] = position + element["length"]
        if new_element != False:
            dna_elements.append(new_element)
        position += element["length"]

    genome = Genome(params["genome"]["name"], position, dna_elements, transcript_template)


    for pol in params["polymerases"]:
        # add polymerase as a reactant
        simulation.increment_reactant(pol["name"], pol["copy_number"])
        # build interactions
        pol["interactions"] = [pol["name"]] # all pols interact with themselves
        for element in params["elements"]:
            if "interactions" in element.keys():
                if pol["name"] in element["interactions"].keys():
                    pol["interactions"].append(element["name"])


    for element in params["elements"]:
        if element["type"] == "promoter":
            simulation.increment_reactant(element["name"], 1)
            for partner, constant in element["interactions"].items():
                binding_constant = constant["binding_constant"]
                interactions = [partner, element["name"]]
                for pol in params["polymerases"]:
                    if pol["name"] == partner:
                        pol_args = [partner,
                                    element["start"],
                                    10,                 # footpring
                                    pol["speed"],
                                    pol["interactions"]
                                    ]
                        reaction = Bind(simulation,
                                        genome,
                                        binding_constant,
                                        interactions,
                                        pol_args)
                        simulation.register_reaction(reaction)

    simulation.register_polymer(genome)

    # Add binding reactions
    # Add species

    #     # Construct polymerases
    #     rna_pol = Polymerase("rna_pol", 15, 10, args.speed[0], interactions)
    #     # Construct features
    #     promoter = Promoter("phi", 1, 10, ["rna_pol"])
    #     terminator = Terminator("T", 90, 100, ["rna_pol"])
    #     elements = [promoter, terminator]
    #     # Construct polymer
    #     tracker = Polymer("dna", 150, elements)
    #     # tracker.bind_polymerase(rna_pol)
    #
    #     simulation.increment_reactant("rna_pol", 10)
    #     simulation.increment_reactant("phi", 1)
    simulation.increment_reactant("ribosome", 100)
    #     pol_args = ["rna_pol", 1, 10, 4, ["rna_pol", "T", "phi"]]
    #     reaction = Bind(simulation, tracker, 0.1, ["rna_pol", "phi"], pol_args)
    #     simulation.register_reaction(reaction)
    #     simulation.register_polymer(tracker)
    simulation.build_heap()

    time_step = params["simulation"]["time_step"]
    old_time = 0
    while(simulation.time < params["simulation"]["runtime"]):
        simulation.execute()
        if abs(simulation.time - old_time) > time_step:
            print(simulation)
            old_time = simulation.time
        elif params["simulation"]["debug"] == True:
            print(simulation)
            print(genome)


if __name__ == "__main__":
    main()
