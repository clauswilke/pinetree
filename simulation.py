#! /usr/bin/env python3

import heapq
import random

class Polymerase:

    def __init__(self, name, start, footprint, speed):
        self.name = name
        self.footprint = footprint
        self.start = start
        self.speed = speed
        self.reenter_queue = True

    def move(self):
        self.start += 1


class Feature:

    def __init__(self, start, stop, interactions):
        self.start = start
        self.stop = stop
        self.interactions = interactions

    def react(self, pol):
        pass

class Terminator(Feature):

    def __init__(self, start, stop, interactions):
        super().__init__(start, stop, interactions)

    def react(self, pol):
        if pol.name in self.interactions:
            pol.reenter_queue = False


class Polymer:

    def __init__(self, name, length, features, polymerase):
        self.name = name
        self.length = length
        self.features = features
        self.polymerases = []
        self.time = 0
        self.heap = []
        self.occupancy = [0]*self.length
        self.feature_locs = [0]*self.length
        self.add_polymerase(polymerase)
        self.build_features()

    def build_features(self):
        for feature in self.features:
            for i in range(feature.start, feature.stop):
                self.feature_locs[i] = 1

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

    def remove_occupancy(self, start, stop):
        for i in range(start, stop + 1):
            self.occupancy[i] = 0

    def calculate_time(self, pol):
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        time, pol = self.pop_polymerase()
        if self.check_collision(pol) == False:
            self.occupancy[pol.start] = 0
            self.occupancy[pol.start + pol.footprint + 1] = 1
            pol.move()
        self.time = time

        if self.check_features(pol) == True:
            for feature in self.features:
                pol_loc = set(range(pol.start, pol.start + pol.footprint))
                feature_loc = set(range(feature.start, feature.stop))
                if pol_loc & feature_loc:
                    feature.react(pol)
                    break

        if pol.reenter_queue == True:
            self.push_polymerase(pol)
        else:
            self.remove_occupancy(pol.start, pol.start + pol.footprint)
            self.polymerases.remove(pol)

    def check_collision(self, pol):
        if self.occupancy[pol.start + pol.footprint + 1] == 0:
            return False
        else:
            return True

    def check_features(self, pol):
        if self.feature_locs[pol.start + pol.footprint + 1] == 0:
            return False
        else:
            return True

    def __str__(self):
        out_string = "Polymerase position: "
        for pol in self.polymerases:
            out_string += str(pol.start) + ", "
        out_string += "time: " + str(self.time)
        out_string += ", occupancy: \n" + ''.join(map(str, self.occupancy))
        out_string += "\nfeatures: \n" + ''.join(map(str, self.feature_locs))
        return out_string


def main():

    # Construct polymerases
    rna_pol = Polymerase("rna_pol", 1, 10, 4)
    rna_pol2 = Polymerase("rna_pol2", 40, 10, 4)

    # Construct features
    interactions = ["rna_pol", "rna_pol2"]
    promoter = Feature(1, 10, [])
    terminator = Terminator(90, 100, interactions)
    features = [promoter, terminator]

    # Construct polymer
    tracker = Polymer("dna", 150, features, rna_pol)
    tracker.add_polymerase(rna_pol2)

    while(tracker.time < 100):
        tracker.execute()
        print(tracker)
        if len(tracker.heap) == 0:
            break

if __name__ == "__main__":
    main()
