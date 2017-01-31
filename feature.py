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
        self.__observers = []

    def register_observer(self, observer):
        self.__observers.append(observer)

    def notify_observers(self, **kwargs):
        for observer in self.__observers:
            observer.notify(self, **kwargs)

    def check_interaction(self, feature):
        return feature.name in self.interactions

    def react(self, pol):
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

    def move_back(self):
        """
        Move one unit backwards
        """
        self.start -= 1
        self.stop -= 1

    def react(self, pol):
        pol.move_back()

class Element(Feature):
    """
    A fixed element in the polymer that can be covered or uncovered.
    """
    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)
        self.covered = False

    def cover(self):
        self.covered = True


class Promoter(Element):

    def __init__(self, name, start, stop, interactions):
        super().__init__(name, start, stop, interactions)

    # def react(self, pol):
    #     self.cover()

class Terminator(Element):
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
