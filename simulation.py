#! /usr/bin/env python3

import heapq
import random
import argparse

from feature import *
from polymer import *

class SpeciesReaction:
    def __init__(self, sim):
        self.sim = sim

    def next_time(self, time):
        """
        Calculate the time at which this reaction will occur next.
        """
        pass

    def execute(self):
        pass

class Bind(SpeciesReaction):
    """
    Bind a polymerase to a polymer.
    """
    def __init__(self, sim, polymer, rate_constant, reactants, product_args):
        super().__init__(sim)
        self.polymer = polymer
        self.polymerase_args = product_args
        self.rate_constant = rate_constant
        self.reactants = reactants

    def next_time(self, time):
        propensity = self.rate_constant
        for reactant in self.reactants:
            propensity *= self.sim.reactants[reactant]

        try:
            delta_time = random.expovariate(propensity)
        except ZeroDivisionError:
            delta_time = float('Inf')

        return time + delta_time

    def execute(self):
        for reactant in self.reactants:
            self.sim.increment_reactant(reactant, -1)
        pol = Polymerase(*self.polymerase_args)
        self.polymer.bind_polymerase(pol)

class Bridge(SpeciesReaction):
    """
    Encapsulate polymer so it can participate in species-level reaction queue.
    """
    def __init__(self, sim, polymer):
        super().__init__(sim)
        self.polymer = polymer

    def next_time(self, time):
        if len(self.polymer.heap) > 0:
            return self.polymer.heap[0][0]
        else:
            return float('Inf')

    def execute(self):
        self.polymer.execute()

    def __str__(self):
        return str(self.polymer)


class Simulation:
    """
    Collect data from polymers.
    """
    def __init__(self):
        self.time = 0
        self.terminations = {}
        self.reactants = {}
        self.reactions = []
        self.heap = []
        self.rebuild_heap = False
        self.iteration = 0

    def increment_reactant(self, name, copy_number):
        if name in self.reactants.keys():
            self.reactants[name] += copy_number
        else:
            self.reactants[name] = copy_number

    def register_reaction(self, reaction):
        if reaction not in self.reactions:
            self.reactions.append(reaction)

    def register_polymer(self, polymer):
        polymer.register_observer(self)
        self.reactions.append(Bridge(self, polymer))

    def build_heap(self):
        self.heap = []
        for i in range(len(self.reactions)):
            self.heap.append((self.reactions[i].next_time(self.time), i))
        heapq.heapify(self.heap)

    def pop(self):
        reaction = heapq.heappop(self.heap)
        return reaction[0], reaction[1]

    def push_reaction(self, reaction_index):
        heapq.heappush(self.heap,
            (self.reactions[reaction_index].next_time(self.time),
             reaction_index))

    def execute(self):
        time, index = self.pop()
        self.time = time
        self.reactions[index].execute()
        self.build_heap()
        self.iteration += 1

    def notify(self, observable, **kwargs):
        if kwargs["action"] == "terminate" and kwargs["species"] == "rna_pol":
            self.increment_reactant(kwargs["species"], 1)
            promoter = Promoter("rbs", 1, 10, ["ribosome"])
            terminator = Terminator("tstop", 90, 100, ["ribosome"])
            elements = [promoter, terminator]
            polymer = Polymer("rna", 150, elements, self)
            self.register_polymer(polymer)
            self.increment_reactant("rbs", 1)
            self.register_reaction(Bind(self, polymer, 0.05, ["rbs", "ribosome"], ["ribosome", 1, 10, 4, ["ribosome", "tstop", "rbs"]]))
            self.count_termination("full_transcript", self.time)
        elif kwargs["action"] == "free_promoter" and kwargs["species"] == "phi":
            self.increment_reactant(kwargs["species"], 1)
        elif kwargs["action"] == "free_promoter" and kwargs["species"] == "rbs":
            self.increment_reactant(kwargs["species"], 1)
        elif kwargs["action"] == "terminate" and kwargs["species"] == "ribosome":
            self.increment_reactant(kwargs["species"], 1)
            self.increment_reactant("rna_pol", 1)
            self.count_termination("full_protein", self.time)

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

    parser = argparse.ArgumentParser(description = "Simulate a single \
        polymerase moving along a piece of DNA.")

    # Add polymerase speed argument
    parser.add_argument("speed", metavar = "s", type = int, nargs = 1,
                        help = "speed in nucleotides per second of polymerase.")

    args = parser.parse_args()

    random.seed(34)

    # Run simulation 50 times
    for i in range(0, 1):
        simulation = Simulation()
        # Construct interactions
        interactions = ["rna_pol", "T", "phi"]
        # Construct polymerases
        rna_pol = Polymerase("rna_pol", 15, 10, args.speed[0], interactions)
        # Construct features
        promoter = Promoter("phi", 1, 10, ["rna_pol"])
        terminator = Terminator("T", 90, 100, ["rna_pol"])
        elements = [promoter, terminator]
        # Construct polymer
        tracker = Polymer("dna", 150, elements, simulation)
        # tracker.bind_polymerase(rna_pol)

        simulation.increment_reactant("rna_pol", 10)
        simulation.increment_reactant("phi", 1)
        simulation.increment_reactant("ribosome", 100)
        pol_args = ["rna_pol", 1, 10, 4, ["rna_pol", "T", "phi"]]
        reaction = Bind(simulation, tracker, 0.1, ["rna_pol", "phi"], pol_args)
        simulation.register_reaction(reaction)
        simulation.register_polymer(tracker)
        simulation.build_heap()

        time_step = 5
        old_time = 0
        while(simulation.time < 100):
            simulation.execute()
            # print(abs(simulation.time - old_time))
            if abs(simulation.time - old_time) > time_step:
                print(simulation)
                old_time = simulation.time


if __name__ == "__main__":
    main()
