#! /usr/bin/env python3

import heapq
import random

class Polymer:

    def __init__(self, name, length, features):
        # Polymers are always modeled individually, so copy_number of 1
        self.length = length
        self.features = features


class Polymerase:

    def __init__(self, name, start, footprint, speed):
        self.name = name
        self.footprint = footprint
        self.start = start
        self.speed = speed

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
        self.time = 0
        self.heap = self.build_heap()

    def build_heap(self):
        heap = []
        for pol in self.polymerases:
            heapq.heappush(heap, (self.calculate_time(pol), pol))
        return heap

    def calculate_time(self, pol):
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        next_pol = heapq.heappop(self.heap)
        next_pol[1].move()
        self.time = next_pol[0]
        heapq.heappush(self.heap,
                       (self.calculate_time(next_pol[1]), next_pol[1]))

    def __str__(self):
        out_string = "Polymerase position: " + str(self.polymerases[0].start)
        out_string += ", time: " + str(self.time)
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
    rna_pol = Polymerase("rna_pol", 1, 10, 10)

    # Construct Tracker
    tracker = Tracker(polymer, rna_pol)

    while(tracker.time < 100):
        tracker.execute()
        print(tracker)

if __name__ == "__main__":
    main()
