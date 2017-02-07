#! /usr/bin/env python3

import heapq
import random

from feature import *

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
        self.polymerases = []
        self.elements = elements
        self.heap = []
        self.__observers = []

    def register_observer(self, observer):
        self.__observers.append(observer)

    def notify_observers(self, **kwargs):
        for observer in self.__observers:
            observer.notify(self, **kwargs)

    def bind_polymerase(self, pol, current_time):
        """
        Bind a `Polymerase` object to the polymer and add it to priority queue.

        :param pol: `Polymerase` object.
        """
        self.polymerases.append(pol)
        self.push(pol, current_time)

    def push(self, pol, current_time):
        """
        Calculate time-until-next reaction and add `Polymerase` object to
        priority queue.

        :param pol: `Polymerase` object.
        """
        heapq.heappush(self.heap, (self.calculate_time(pol, current_time), pol))

    def pop(self):
        """
        Remove and return `Polymerase` object from top of priority queue, along
        with its reaction time.

        :returns: time, `Polymerase` object
        """
        pol = heapq.heappop(self.heap)
        return pol[0], pol[1]

    def get_next_time(self):
        if len(self.heap) > 0:
            return self.heap[0][0]
        else:
            return float('Inf')

    def calculate_time(self, pol, current_time):
        """
        Calculate time-until-next reaction from an exponential distribution
        centered at a `Polymerase` object's `speed` attribute. Adds time to
        current simulation time.

        :param pol: `Polymerase` object.
        :returns: time that `pol` will move next.
        """
        return current_time + random.expovariate(pol.speed)

    def move_polymerase(self, pol):
        # Record old covered elements
        old_covered_elements = self.find_intersections(pol, self.elements)

        # Move polymerase
        pol.move()

        # Record new collisions
        collisions = self.find_intersections(pol, self.polymerases)

        # First resolve any collisions between polymerases
        for other_pol in collisions:
            other_pol.react(pol)

        # Now cover elements
        new_covered_elements = self.find_intersections(pol, self.elements)
        for element in new_covered_elements:
            element.cover()
            element.react(pol)

        # Check for uncovered elements
        for element in old_covered_elements:
            if element not in new_covered_elements:
                self.notify_observers(species = element.name,
                                      type = element.type,
                                      action = "free_promoter")
                # print("free promoter!")

    def execute(self):
        """
        Process `Polymerase` object at the top of the priority queue. Check for
        collisions, uncovering of elements, and terminations.
        """
        time, pol = self.pop()

        self.move_polymerase(pol)

        if pol.attached == True:
            self.push(pol, time)
            # print("move!")
        else:
            self.notify_observers(species = pol.name,
                                  action = "terminate",
                                  type = pol.type)
            self.polymerases.remove(pol)
            # print("terminate!")


    def find_intersections(self, pol, elements):
        """
        Detect overlap between a given `Polymerase` object and all other
        `Feature` objects in `elements` currently in or on the polymer that it
        interacts with.

        :param pol: `Polymerase` object with which to check for collisions.
        :returns: list of `Feature` objects that overlap and interact with
            `pol`.
        """
        intersections = []
        for element in elements:
            if pol == element:
                continue
            if self.segments_intersect(pol.start, pol.stop,
                                       element.start, element.stop):
                if element.check_interaction(pol):
                    intersections.append(element)
        return intersections

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
        feature_locs = [0]*self.length
        for feature in self.polymerases:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = 1
        out_string += "\nfeatures: \n" + ''.join(map(str, feature_locs))
        return out_string

class Genome(Polymer):

    def __init__(self, name, length, elements, transcript_template):
        super().__init__(name, length, elements)
        self.transcript_template = transcript_template

    def execute(self):
        """
        Process `Polymerase` object at the top of the priority queue. Check for
        collisions, uncovering of elements, and terminations.
        """
        time, pol = self.pop()

        self.move_polymerase(pol)

        if pol.attached == True:
            self.push(pol, time)
            # print("move!")
        else:
            polymer, species = self.build_transcript(pol.bound, pol.stop)
            self.notify_observers(species = pol.name,
                                  action = "terminate_transcript",
                                  type = pol.type,
                                  polymer = polymer,
                                  reactants = species)
            self.polymerases.remove(pol)
            # print("terminate!")

    def build_transcript(self, start, stop):
        """
        Build a transcript object corresponding to this genome.
        """
        species = []
        elements = []
        for element in self.transcript_template:
            if element["start"] >= start and element["stop"] <= stop:
                rbs = Promoter("rbs",
                               element["start"]-element["rbs"],
                               element["start"],
                               ["ribosome"])
                stop_site = Terminator("tstop",
                                       element["stop"],
                                       element["stop"] + 1,
                                       ["ribosome"])
                elements.append(rbs)
                elements.append(stop_site)
                species.append("rbs")
            polymer = Polymer("rna", 150, elements)
        return polymer, species

class Transcript(Polymer):
    def __init__():
        pass
