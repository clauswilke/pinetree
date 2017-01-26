#! /usr/bin/env python3

import heapq
import random

class Feature:
    """
    A generic feature in or on `Polymer`. Designed to be extended
    by `Terminator`, `Promoter`, `Polymerase`, etc.
    """
    def __init__(self, name, start, stop, interactions):
        self.name = name
        self.start = start
        self.stop = stop
        self.interactions = interactions

    def react(self, pol):
        """
        React with `pol`.

        :param pol: `Polymerase` interacting with this `Feature`.
        """
        pass

class Polymerase(Feature):
    """
    A molecule that binds to `Polymer` and moves. Extends `Feature`.
    """
    def __init__(self, name, start, footprint, speed):
        super().__init__(name, start, start + footprint, [])
        self.speed = speed
        self.attached = True

    def move(self):
        """
        Move one unit forward.
        """
        self.start += 1
        self.stop += 1

class Terminator(Feature):
    """
    Stops movement of `Polymerase` along a `Polymer`. Extends `Feature`.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)

    def react(self, pol):
        """
        Check for interaction with `pol` and detach `pol`.

        :param pol: `Polymerase`.
        """
        if pol.name in self.interactions:
            pol.attached = False

class Polymer:
    """
    Track `Feature` objects, `Polymerase` objects, and collisions on a single
    polymer. Move `Polymerase` objects along the polymer, maintain a priority
    queue of `Polymerase` objects that are expected to move, and calculate
    time-until-next move from an exponential distribution.
    """
    def __init__(self, name, length, features, polymerase):
        self.name = name
        self.length = length
        self.features = features
        self.polymerases = []
        self.time = 0
        self.heap = []
        self.bind_polymerase(polymerase)

    def bind_polymerase(self, pol):
        """
        Bind a `Polymerase` object to the polymer and add it to priority queue.

        :param pol: `Polymerase` object.
        """
        self.polymerases.append(pol)
        self.push(pol)

    def push(self, pol):
        """
        Calculate time-until-next reaction and add `Polymerase` object to
        priority queue.

        :param pol: `Polymerase` object.
        """
        heapq.heappush(self.heap, (self.calculate_time(pol), pol))

    def pop(self):
        """
        Remove and return `Polymerase` object from top of priority queue, along
        with its reaction time.

        :returns: time, `Polymerase` object
        """
        pol = heapq.heappop(self.heap)
        return pol[0], pol[1]

    def calculate_time(self, pol):
        """
        Calculate time-until-next reaction from an exponential distribution
        centered at a `Polymerase` object's `speed` attribute. Adds time to
        current simulation time.

        :param pol: `Polymerase` object.
        :returns: time that `pol` will move next.
        """
        return self.time + random.expovariate(pol.speed)

    def execute(self):
        """
        Process `Polymerase` object at the top of the priority queue. Check for
        collisions and terminations.
        """
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
        """
        Detect collision between a given `Polymerase` object and all other
        `Polymerase` objects currently on the polymer.

        :param pol: `Polymerase` object with which to check for collisions.
        :returns: interacting `Polymerase` object, or None.
        """
        for pol2 in self.polymerases:
            if self.segments_intersect(pol.start, pol.stop,
                                       pol2.start, pol2.stop):
                if pol2 != pol:
                    return pol2
        return None

    def find_feature(self, pol):
        """
        Detect collision between a given `Polymerase` object and all `Feature`
        objects currently in the polymer.

        :param pol: `Polymerase` object with which to check for collisions.
        :returns: interacting `Feature` object, or None.
        """
        for feature in self.features:
            if self.segments_intersect(pol.start, pol.stop,
                                       feature.start, feature.stop):
                return feature
        return None

    def segments_intersect(self, x1, x2, y1, y2):
        """
        Do two line segments (e.g. `Polymerase` objects) overlap?
        """
        return x2 >= y1 and y2 >= x1

    def __str__(self):
        """
        Convert `Polymer` object to string representation showing features and
        polymerases.
        """
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
    promoter = Feature("phi", 1, 10, [])
    terminator = Terminator("T", 90, 100, interactions)
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
