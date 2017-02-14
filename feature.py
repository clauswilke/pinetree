#! /usr/bin/env python3

from signal import *

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
        self.type = "" # type of feature, i.e., polymerase, promoter, etc.

    def check_interaction(self, feature):
        """
        Check to see if some other feature interacts with this feature.

        :param feature: a feature object
        """
        return feature.name in self.interactions

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
        self.attached = True # Is this polymerase attached to a polymer?
        self.bound = start # Record where polymerase bound to genome
        self.type = "polymerase"
        self.footprint = footprint
        self.move_signal = Signal()

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

    def react(self, pol):
        """
        Move polymerase back one position of it collides with another
        polymerase.
        """
        pol.move_back()

class Mask(Feature):
    """
    A pseudo-feature that tracks which portion of a genome or polymer are not
    yet accessible. For example, as the genome is entering the cell, or as a
    transcript is being synthesized.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)

    def react(self, pol):
        """
        Shrink mask by one base pair.

        TODO: have a dynamic step size?
        """
        self.start += 1

class TranscriptMask(Feature):
    """
    A pseudo-feature that tracks which portion of a genome or polymer are not
    yet accessible. For example, as the genome is entering the cell, or as a
    transcript is being synthesized.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)

    def react(self, pol):
        """
        Shrink mask by one base pair.

        TODO: have a dynamic step size?
        """
        pol.move_back()

class Element(Feature):
    """
    A fixed feature in the polymer that can be covered or uncovered.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)
        self.covered = 0 # is this element covered? (i.e. inaccessible)
        self.old_covered = 0
        self.type = ""

    def save_state(self):
        self.old_covered = self.covered

    def was_uncovered(self):
        return self.old_covered >= 1 and self.covered == 0

    def was_covered(self):
        return self.old_covered == 0 and self.covered > 0

    def cover(self):
        """
        Cover this element.
        """
        self.covered += 1

    def uncover(self):
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
        super().__init__(name, start, stop, interactions)
        self.type = "terminator"
        self.gene = ""

    def react(self, pol):
        """
        Check for interaction with `pol` and detach `pol`.

        :param pol: `Polymerase`.
        """
        pol.attached = False # signal to polymer to destroy polymerase
        # tell polymerase the last gene that it transcribed so it can construct
        # the correct transcript
        pol.last_gene = self.gene
