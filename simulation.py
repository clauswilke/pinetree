#! /usr/bin/env python3

import heapq
import random

class Polymerase:

    def __init__(self, name, start, footprint, speed):
        self.name = name
        self.footprint = footprint
        self.start = start
        self.stop = start + footprint
        self.speed = speed
        self.attached = True

    def move(self):
        self.start += 1
        self.stop += 1


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
            pol.attached = False


class Polymer:

    def __init__(self, name, length, features, polymerase):
        self.name = name
        self.length = length
        self.features = features
        self.polymerases = []
        self.time = 0
        self.heap = []
        self.bind_polymerase(polymerase)

    def bind_polymerase(self, pol):
        self.polymerases.append(pol)
        self.push(pol)

    def push(self, pol):
        heapq.heappush(self.heap, (self.calculate_time(pol), pol))

    def pop(self):
        pol = heapq.heappop(self.heap)
        return pol[0], pol[1]

    def calculate_time(self, pol):
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        time, pol = self.pop()

        collision = self.find_collision(pol)
        if collision == None:
            pol.move()
        self.time = time

        feature = self.find_feature(pol)
        if feature != None:
            feature.react(pol)

        if pol.attached == True:
            self.push(pol)
        else:
            self.polymerases.remove(pol)

    def find_collision(self, pol):
        for pol2 in self.polymerases:
            if self.segments_intersect(pol.start, pol.stop,
                                       pol2.start, pol2.stop):
                if pol2 != pol:
                    return pol2
        return None

    def find_feature(self, pol):
        for feature in self.features:
            if self.segments_intersect(pol.start, pol.stop,
                                       feature.start, feature.stop):
                return feature
        return None

    def segments_intersect(self, x1, x2, y1, y2):
        return x2 >= y1 and y2 >= x1

    def __str__(self):
        out_string = "Polymerase position: "
        polymerase_locs = polymerase_locs = [0]*self.length
        for pol in self.polymerases:
            out_string += str(pol.start) + ", "
            for i in range(pol.start - 1, pol.stop - 1):
                polymerase_locs[i] = 1
        out_string += "time: " + str(self.time)

        out_string += ", occupancy: \n" + ''.join(map(str, polymerase_locs))

        feature_locs = [0]*self.length
        for feature in self.features:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = 1
        out_string += "\nfeatures: \n" + ''.join(map(str, feature_locs))
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
    tracker.bind_polymerase(rna_pol2)

    while(tracker.time < 100):
        tracker.execute()
        print(tracker)
        if len(tracker.heap) == 0:
            break

if __name__ == "__main__":
    main()
