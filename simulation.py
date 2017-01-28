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
    def __init__(self, sim, polymer, rate_constant):
        super().__init__(sim)
        self.polymer = polymer
        self.rate_constant = rate_constant

    def next_time(self, time):
        propensity = self.sim.reactants["rna_pol"] * \
            self.sim.reactants["phi"] * self.rate_constant

        try:
            delta_time = random.expovariate(propensity)
        except ZeroDivisionError:
            delta_time = float('Inf')

        return time + delta_time

    def execute(self):
        self.sim.register_reactant("rna_pol", -1)
        self.sim.register_reactant("phi", -1)
        pol = Polymerase("rna_pol", 1, 10, 4, ["rna_pol", "T"])
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

    def register_reactant(self, name, copy_number):
        if name in self.reactants.keys():
            self.reactants[name] += copy_number
        else:
            self.reactants[name] = copy_number

    def register_reaction(self, reaction):
        self.reactions.append(reaction)

    def build_heap(self):
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
        print(self)
        self.time = time
        print("Executing reaction " + str(index))
        self.reactions[index].execute()
        self.push_reaction(index)


    def count_termination(self, name, time):
        """
        Record the time at which a polymerase reaches a terminator.
        """
        if name not in self.terminations.keys():
            self.terminations[name] = time

    def __str__(self):
        """
        Convert `Simulation` to string representation showing names and
        termination times of `Polymerase` objects.
        """
        out_string = ""
        for name, time in self.terminations.items():
            out_string += name + ", " + str(time)
        out_string += str(self.heap) + str(self.reactants)
        return out_string

def main():

    parser = argparse.ArgumentParser(description = "Simulate a single \
        polymerase moving along a piece of DNA.")

    # Add polymerase speed argument
    parser.add_argument("speed", metavar = "s", type = int, nargs = 1,
                        help = "speed in nucleotides per second of polymerase.")

    args = parser.parse_args()

    # Run simulation 50 times
    for i in range(0, 1):
        simulation = Simulation()
        # Construct interactions
        interactions = ["rna_pol", "T"]
        # Construct polymerases
        rna_pol = Polymerase("rna_pol", 15, 10, args.speed[0], interactions)
        # Construct features
        # promoter = Feature("phi", 1, 10, [])
        terminator = Terminator("T", 90, 100, interactions)
        elements = [terminator]
        # Construct polymer
        tracker = Polymer("dna", 150, elements)
        tracker.register_sim(simulation)
        tracker.bind_polymerase(rna_pol)

        simulation.register_reactant("rna_pol", 10)
        simulation.register_reactant("phi", 1)
        reaction = Bind(simulation, tracker, 0.1)
        reaction2 = Bridge(simulation, tracker)
        simulation.register_reaction(reaction)
        simulation.register_reaction(reaction2)
        simulation.build_heap()


        while(len(tracker.heap) > 0 and simulation.time < 50):
            simulation.execute()
            print(simulation)
            print(tracker)


if __name__ == "__main__":
    main()
