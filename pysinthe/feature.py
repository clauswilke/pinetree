#! /usr/bin/env python3

"""
Defines promoters, terminators, polymerases, and other objects that are either
fixed components of a polymer or move on a polymer.
"""

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


class Polymerase:
    """
    A molecule that binds to `Polymer` and moves.
    """
    def __init__(self, name, footprint, speed):
        """
        :param name: name of polymerase (unique?)
        :param footprint: polymerase footprint
        :param speed: speed of polymerase
        """
        self.name = name
        self.start = 0
        self.stop = footprint - 1
        self.speed = speed
        self.left_most_element = 0
        self.bound = 0  # Record where polymerase bound to genome
        self.type = "polymerase"
        self.footprint = footprint
        self.reading_frame = 0
        self.move_signal = eventsignal.Signal()  # signal to fire when this
        # polymerase moves
        self.release_signal = eventsignal.Signal()

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
        self._covered = 0  # is this element covered? (i.e. inaccessible)
        self._old_covered = 0
        self.type = ""
        self.cover_signal = eventsignal.Signal()
        self.uncover_signal = eventsignal.Signal()

    def save_state(self):
        """
        Save covering state.
        """
        self._old_covered = self._covered

    def was_uncovered(self):
        """
        Was this element just uncovered?
        """
        return self._old_covered >= 1 and self._covered == 0

    def was_covered(self):
        """
        Was this element just covered?
        """
        return self._old_covered == 0 and self._covered > 0

    def cover(self):
        """
        Cover this element. Elements can be covered by multiple features.
        """
        self._covered += 1

    def uncover(self):
        """
        Uncover element
        """
        if self._covered > 0:
            self._covered -= 1

    def is_covered(self):
        """
        Is this element covered?
        """
        return self._covered > 0

    def check_state(self):
        """
        Check for a change in state and react appropriately.
        """
        return


class Promoter(Element):
    """
    A promoter class.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)
        self.type = "promoter"

    def check_state(self):
        if self.was_covered():
            self.cover_signal.fire(self.name)
        elif self.was_uncovered():
            self.uncover_signal.fire(self.name)


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
        self.reading_frame = 0

    def check_state(self):
        if self.was_uncovered():
            self.readthrough = False
            self.uncover_signal.fire(self.name)
        elif self.was_covered():
            self.cover_signal.fire(self.name)

    def check_interaction(self, feature_name, reading_frame):
        """
        Check to see if some other feature interacts with this feature.

        :param feature: a feature object
        """
        return reading_frame is self.reading_frame and \
            feature_name in self.interactions
