#! /usr/bin/env python3

class Molecule:

    def __init__(self, name, copy_number):
        self.name = name
        self.copy_number = copy_number

class Polymer(Molecule):

    def __init__(self, name, length, features):
        # Polymers are always modeled individually, so copy_number of 1
        super().__init__(name, 1)
        self.length = length
        self.features = features


class Polymerase(Molecule):

    def __init__(self):
        pass

class Feature:

    def __init__(start, stop, interactions):
        self.start = start
        self.stop = stop
        pass

class Reaction:

    def __init__(self, reactants):
        self.reactants = reactant

    def execute(self):
        pass

class Tracker(Reaction):

    def __init__(self, polymer, polymerases):
        self.polymer = polymer
        self.polymerases = polymerases
