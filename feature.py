#! /usr/bin/env python3

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
