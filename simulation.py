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
        self.polymerases = []
        self.time = 0
        self.heap = []
        self.occupancy = [0]*polymer.length
        self.add_polymerase(polymerase)

    def add_polymerase(self, pol):
        self.polymerases.append(pol)
        for i in range(pol.start, pol.start + pol.footprint + 1):
            self.occupancy[i] = 1
        heapq.heappush(self.heap, (self.calculate_time(pol), pol))

    def calculate_time(self, pol):
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        next_pol = heapq.heappop(self.heap)
        if self.check_collision(next_pol[1]) == False:
            self.occupancy[next_pol[1].start] = 0
            self.occupancy[next_pol[1].start + next_pol[1].footprint + 1] = 1
            next_pol[1].move()
        self.time = next_pol[0]
        heapq.heappush(self.heap,
                       (self.calculate_time(next_pol[1]), next_pol[1]))

    def check_collision(self, pol):
        if self.occupancy[pol.start + pol.footprint + 1] == 0:
            return False
        else:
            return True


    def __str__(self):
        out_string = "Polymerase position: "
        for pol in self.polymerases:
            out_string += str(pol.start) + ", "
        out_string += "time: " + str(self.time)
        out_string += ", occupancy: \n" + ''.join(map(str, self.occupancy))
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
    rna_pol2 = Polymerase("rna_pol2", 40, 10, 1)

    # Construct Tracker
    tracker = Tracker(polymer, rna_pol)
    tracker.add_polymerase(rna_pol2)

    while(tracker.time < 20):
        tracker.execute()
        print(tracker)

if __name__ == "__main__":
    main()
