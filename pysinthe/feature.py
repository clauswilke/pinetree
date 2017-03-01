#! /usr/bin/env python3

"""
Defines promoters, terminators, polymerases, and other objects that are either
fixed components of a polymer or move on a polymer.
"""

import random

from . import eventsignal


class Feature:
    """
    A generic feature in or on `Polymer`. Designed to be extended
    by `Terminator`, `Promoter`, `Polymerase`, etc.

    TODO: make abstract?
    """
    def __init__(self, name, start, stop, interactions):
        """
        :param name: a (unique?) name for this feature
        :param start: start position of this feature within polymer
        :param stop: stop position of this feature within polymer
        :param interactions: list of names of other features/polymerases that
            this feature interacts with
        """
        self.name = name
        self.start = start
        self.stop = stop
        self.interactions = interactions
        self.type = ""  # type of feature, i.e., polymerase, promoter, etc.

    def check_interaction(self, feature_name):
        """
        Check to see if some other feature interacts with this feature.

        :param feature: a feature object
        """
        return feature_name in self.interactions

    def react(self, pol):
        """
        React with a feature.

        :param pol: a feature object
        """
        pass


class Polymerase(Feature):
    """
    A molecule that binds to `Polymer` and moves. Extends `Feature`.
    """
    def __init__(self, name, start, footprint, speed, interactions):
        """
        :param name: name of polymerase (unique?)
        :param start: current start position of polymerase
        :param footprint: polymerase footprint
        :param speed: speed of polymerase
        :param interactions: list of other features that this polymerase
            interacts with
        """
        super().__init__(name, start, start + footprint, interactions)
        self.speed = speed
        self.attached = True  # Is this polymerase attached to a polymer?
        self.bound = start  # Record where polymerase bound to genome
        self.type = "polymerase"
        self.footprint = footprint
        self.move_signal = eventsignal.Signal()  # signal to fire when this
        # polymerase moves
        self.termination_signal = eventsignal.Signal()

    def move(self):
        """
        Move one unit forward.
        """
        self.start += 1
        self.stop += 1

    def move_back(self):
        """
        Move one unit backwards
        """
        self.start -= 1
        self.stop -= 1


class Mask(Feature):
    """
    A pseudo-feature that tracks which portion of a genome or polymer are not
    yet accessible. For example, as the genome is entering the cell, or as a
    transcript is being synthesized.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)

    def recede(self):
        """
        Shrink mask by one base pair.

        TODO: have a dynamic step size?
        """
        self.start += 1


class Element(Feature):
    """
    A fixed feature in the polymer that can be covered or uncovered.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)
        self.covered = 0  # is this element covered? (i.e. inaccessible)
        self.old_covered = 0
        self.type = ""

    def save_state(self):
        """
        Save covering state.
        """
        self.old_covered = self.covered

    def was_uncovered(self):
        """
        Was this element just uncovered?
        """
        return self.old_covered >= 1 and self.covered == 0

    def was_covered(self):
        """
        Was this element just covered?
        """
        return self.old_covered == 0 and self.covered > 0

    def cover(self):
        """
        Cover this element. Elements can be covered by multiple features.
        """
        self.covered += 1

    def uncover(self):
        """
        Uncover element
        """
        if self.covered > 0:
            self.covered -= 1


class Promoter(Element):
    """
    A promoter class.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)
        self.type = "promoter"


class Terminator(Element):
    """
    Stops movement of `Polymerase` along a `Polymer`.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions.keys())
        self.type = "terminator"
        self.gene = ""
        self.efficiency = interactions
        self.readthrough = False

    def react(self, pol):
        """
        Check for interaction with `pol` and detach `pol`.

        TODO: add "attached" signal to event-firing system?

        :param pol: `Polymerase`.
        """
        if self.readthrough:
            return
        random_num = random.random()
        if random_num <= self.efficiency[pol.name]["efficiency"]:
            pol.attached = False  # signal to polymer to destroy polymerase
            # tell polymerase the last gene that it transcribed so it can
            # construct the correct transcript
            pol.last_gene = self.gene
            # Uncover terminator, mostly for debugging purposes
            pol.termination_signal.fire(self.stop)
            self.uncover()
        else:
            self.readthrough = True
