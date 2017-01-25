#! /usr/bin/env python3

class Polymer:

    def __init__(self, name, length, features):
        # Polymers are always modeled individually, so copy_number of 1
        self.length = length
        self.features = features


class Polymerase:

    def __init__(self, name, footprint):
        self.name = name
        self.footprint = footprint

class Feature:

    def __init__(self, start, stop, interactions):
        self.start = start
        self.stop = stop
        self.interactions = interactions

class Tracker:

    def __init__(self, polymer, polymerase):
        self.polymer = polymer
        self.polymerases = polymerase

class PriorityQueue:

    def __init__(self):
        pass

def main():

    # Construct features
    interactions = ["rna_pol"]
    promoter = Feature(1, 10, interactions)
    terminator = Feature(90, 100, interactions)

    # Construct polymer
    features = [promoter, terminator]
    polymer = Polymer("dna", 150, features)

    # Construct polymerases
    rna_pol = Polymerase("rna_pol", 10)

    # Construct Tracker
    tracker = Tracker(polymer, rna_pol)

    # Construct PriorityQueue
    queue = PriorityQueue()

if __name__ == "__main__":
    main()
