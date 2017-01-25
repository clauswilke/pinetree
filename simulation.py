#! /usr/bin/env python3

import heapq
import random

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

class Polymer:

    def __init__(self, name, length, features, polymerase):
        self.name = name
        self.length = length
        self.features = features
        self.polymerases = []
        self.time = 0
        self.heap = []
        self.occupancy = [0]*self.length
        self.add_polymerase(polymerase)

    def add_polymerase(self, pol):
        self.polymerases.append(pol)
        for i in range(pol.start, pol.start + pol.footprint + 1):
            self.occupancy[i] = 1
        self.push_polymerase(pol)

    def push_polymerase(self, pol):
        heapq.heappush(self.heap, (self.calculate_time(pol), pol))

    def pop_polymerase(self):
        pol = heapq.heappop(self.heap)
        return pol[0], pol[1]

    def calculate_time(self, pol):
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        time, pol = self.pop_polymerase()
        if self.check_collision(pol) == False:
            self.occupancy[pol.start] = 0
            self.occupancy[pol.start + pol.footprint + 1] = 1
            pol.move()
        self.time = time
        self.push_polymerase(pol)

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
    features = [promoter, terminator]

    # Construct polymerases
    rna_pol = Polymerase("rna_pol", 1, 10, 10)
    rna_pol2 = Polymerase("rna_pol2", 40, 10, 1)

    # Construct polymer
    tracker = Polymer("dna", 150, features, rna_pol)
    tracker.add_polymerase(rna_pol2)

    while(tracker.time < 20):
        tracker.execute()
        print(tracker)

if __name__ == "__main__":
    main()
