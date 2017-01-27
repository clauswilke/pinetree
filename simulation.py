#! /usr/bin/env python3

import heapq
import random
import argparse

from feature import *
from polymer import *

class SpeciesReaction:
    def __init__(self):
        pass

    def next_time(self):
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
    def __init__(self):
        super().__init__()

    def execute(self):
        # decrement free polymerase, add pol to polymer
        pass

class Bridge(SpeciesReaction):
    """
    Encapsulate polymer so it can participate in species-level reaction queue.
    """
    def __init__(self, polymer):
        super().__init__()
        self.polymer = polymer

    def next_time(self):
        return polymer.heap[0][0]

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

    def register_reaction(self, reactant1, reactant2, rate_constant):
        self.reactions.append({"reactants": [reactant1, reactant2],
                               "rate_constant": rate_constant})

    def register_polymer(self, polymer):
        self.heap.append((polymer.heap[0][0], polymer))

    def build_heap(self):
        for i in range(len(self.reactions)):
            self.heap.append((self.calculate_time(self.reactions[i]), i))
        heapq.heapify(self.heap)

    def pop(self):
        reaction = heapq.heappop(self.heap)
        return reaction[0], reaction[1]

    def push_reaction(self, reaction_index):
        heapq.heappush(self.heap,
            (self.calculate_time(self.reactions[reaction_index]),
             reaction_index))

    def push_polymer(self, polymer):
        heapq.heappush(self.heap, (polymer.heap[0][0], polymer))

    def calculate_time(self, reaction):
        propensity = self.reactants[reaction["reactants"][0]] * \
            self.reactants[reaction["reactants"][1]] * reaction["rate_constant"]

        return self.time + random.expovariate(propensity)

    def execute(self):
        time, index = self.pop()
        print(self)
        self.time = time
        if type(index) is int:
            print("Executing reaction " + str(index))
            self.push_reaction(index)
        else:
            print("Executing polymer")
            index.execute()
            self.push_polymer(index)


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
        rna_pol = Polymerase("rna_pol", 1, 10, args.speed[0], interactions)
        # Construct features
        promoter = Feature("phi", 1, 10, [])
        terminator = Terminator("T", 90, 100, interactions)
        elements = [promoter, terminator]
        # Construct polymer
        tracker = Polymer("dna", 150, elements)
        tracker.register_sim(simulation)
        tracker.bind_polymerase(rna_pol)

        simulation.register_reactant("rna_pol", 10)
        simulation.register_reactant("phi", 1)
        simulation.register_reaction("rna_pol", "phi", 0.1)
        simulation.register_polymer(tracker)
        simulation.build_heap()

        """
        while(len(tracker.heap) != 0):
            tracker.execute()
        """

        print(simulation)

        simulation.execute()

        print(simulation)

        simulation.execute()

        print(simulation)

if __name__ == "__main__":
    main()
