#! /usr/bin/env python3

import heapq
import random

class Polymer:
    """
    Track `Feature` objects, `Polymerase` objects, and collisions on a single
    polymer. Move `Polymerase` objects along the polymer, maintain a priority
    queue of `Polymerase` objects that are expected to move, and calculate
    time-until-next move from an exponential distribution.
    """
    def __init__(self, name, length, elements):
        self.name = name
        self.length = length
        self.features = elements
        self.heap = []

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
        return self.sim.time + random.expovariate(pol.speed)

    def execute(self):
        """
        Process `Polymerase` object at the top of the priority queue. Check for
        collisions and terminations.
        """
        time, pol = self.pop()
        pol.move()
        collisions = self.find_collisions(pol)
        if len(collisions) != 0:
            for feature in collisions:
                feature.react(pol)
        self.sim.time = time

        if pol.attached == True:
            self.push(pol)
        else:
            self.sim.count_termination(pol.name, self.sim.time)
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
        out_string = "Time: " + str(self.sim.time)

        feature_locs = [0]*self.length
        for feature in self.features:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = 1
        out_string += "\nfeatures: \n" + ''.join(map(str, feature_locs))
        return out_string
