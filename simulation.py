#! /usr/bin/env python3

import heapq
import random
import argparse

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

    def check_interaction(self, feature):
        return feature.name in self.interactions

    def react(self):
        pass

class Polymerase(Feature):
    """
    A molecule that binds to `Polymer` and moves. Extends `Feature`.
    """
    def __init__(self, name, start, footprint, speed, interactions):
        super().__init__(name, start, start + footprint, interactions)
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
        pol.attached = False

class Polymer:
    """
    Track `Feature` objects, `Polymerase` objects, and collisions on a single
    polymer. Move `Polymerase` objects along the polymer, maintain a priority
    queue of `Polymerase` objects that are expected to move, and calculate
    time-until-next move from an exponential distribution.
    """
    def __init__(self, name, length, elements, polymerase):
        self.name = name
        self.length = length
        self.features = elements
        self.time = 0
        self.heap = []
        self.bind_polymerase(polymerase)

    def register_sim(self, sim):
        self.sim = sim

    def bind_polymerase(self, pol):
        """
        Bind a `Polymerase` object to the polymer and add it to priority queue.

        :param pol: `Polymerase` object.
        """
        self.features.append(pol)
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

        collisions = self.find_collisions(pol)
        if len(collisions) == 0:
            pol.move()
        else:
            for feature in collisions:
                feature.react(pol)
        self.time = time

        if pol.attached == True:
            self.push(pol)
        else:
            self.sim.count_termination(pol.name, self.time)
            self.features.remove(pol)

    def find_collisions(self, feature):
        """
        Detect collision between a given `Feature` object and all other
        `Feature` objects currently in or on the polymer.

        :param pol: `Feature` object with which to check for collisions.
        :returns: list of `Feature` objects.
        """
        collisions = []
        for feature2 in self.features:
            if feature == feature2:
                continue
            if self.segments_intersect(feature.start, feature.stop,
                                       feature2.start, feature2.stop):
                if feature2.check_interaction(feature):
                    collisions.append(feature2)
        return collisions

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
        out_string = "Time: " + str(self.time)

        feature_locs = [0]*self.length
        for feature in self.features:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = 1
        out_string += "\nfeatures: \n" + ''.join(map(str, feature_locs))
        return out_string

class Simulation:
    """
    Collect data from polymers.
    """
    def __init__(self):
        self.time = 0
        self.terminations = {}

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
        return out_string

def main():

    parser = argparse.ArgumentParser(description = "Simulate a single \
        polymerase moving along a piece of DNA.")

    # Add polymerase speed argument
    parser.add_argument("speed", metavar = "s", type = int, nargs = 1,
                        help = "speed in nucleotides per second of polymerase.")

    args = parser.parse_args()

    # Run simulation 50 times
    for i in range(0, 1000):
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
        tracker = Polymer("dna", 150, elements, rna_pol)
        tracker.register_sim(simulation)

        while(len(tracker.heap) != 0):
            tracker.execute()

        print(simulation)

if __name__ == "__main__":
    main()
