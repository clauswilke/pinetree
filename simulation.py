#! /usr/bin/env python3

class Polymer:

    def __init__(self, name, length, features):
        # Polymers are always modeled individually, so copy_number of 1
        self.length = length
        self.features = features


class Polymerase:

    def __init__(self, name, start, footprint):
        self.name = name
        self.footprint = footprint
        self.start = start

    def move(self):
        self.start += 1

    def terminate(self):
        pass

class Feature:

    def __init__(self, start, stop, interactions):
        self.start = start
        self.stop = stop
        self.interactions = interactions

class Tracker:

    def __init__(self, polymer, polymerase):
        self.polymer = polymer
        self.polymerases = [polymerase]

    def execute(self):
        next_polymerase = self.polymerases.pop()
        next_polymerase.move()
        self.polymerases.append(next_polymerase)

    def __str__(self):
        out_string = "Polymerase position: " + str(self.polymerases[0].start)
        return out_string


def main():

    # Construct features
    interactions = ["rna_pol"]
    promoter = Feature(1, 10, interactions)
    terminator = Feature(90, 100, interactions)

    # Construct polymer
    features = [promoter, terminator]
    polymer = Polymer("dna", 150, features)

    # Construct polymerases
    rna_pol = Polymerase("rna_pol", 1, 10)

    # Construct Tracker
    tracker = Tracker(polymer, rna_pol)

    print(tracker)

    tracker.execute()

    print(tracker)

if __name__ == "__main__":
    main()
