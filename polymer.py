#! /usr/bin/env python3

import random

from feature import *
from signal import *

class Polymer:
    """
    Track `Feature` objects, `Polymerase` objects, and collisions on a single
    polymer. Move `Polymerase` objects along the polymer. Handle logic for
    covering and uncovering of elements.

    TODO: make class abstract?
    """
    def __init__(self, name, length, elements):
        """
        :param name: name of this polymer (should it be unique?)
        :param length: length of polymer (used purely for debugging, may not be
            needed)
        :param elements: all elements on this polymer, including promoters,
            terminators, etc.
        """
        self.name = name
        self.length = length
        self.polymerases = []
        self.elements = elements
        self.termination_signal = Signal()
        self.promoter_signal = Signal()
        self.block_signal = Signal()
        self.mask = False

    def bind_polymerase(self, pol):
        """
        Bind a `Polymerase` object to the polymer and add it to min-heap.

        :param pol: `Polymerase` object.
        """
        self.polymerases.append(pol)

    def calculate_propensity(self):
        """
        Calculate the total propensity of all polymerase movement in this
        polymer.

        :returns: total propensity
        """
        prop = 0
        for pol in self.polymerases:
            prop += pol.speed
        return prop

    def move_polymerase(self, pol):
        """
        Move polymerase and deal with collisions and covering/uncovering of
        elements.

        :param pol: polymerase to move
        """
        # Find which elements this polymerase is covering and temporarily
        # uncover them
        for element in self.elements:
            if self.segments_intersect(pol.start, pol.stop,
                                       element.start, element.stop):
                element.uncover()
            if self.mask != False:
                if self.segments_intersect(self.mask.start, self.mask.stop,
                                           element.start, element.stop):
                    element.uncover()

        # Move polymerase
        pol.move()

        # First resolve any collisions between polymerases
        for other_pol in self.polymerases:
            if pol == other_pol:
                continue
            if self.segments_intersect(pol.start, pol.stop,
                                       other_pol.start, other_pol.stop):
                if other_pol.check_interaction(pol):
                    other_pol.react(pol)

        if self.mask != False:
            if self.segments_intersect(pol.start, pol.stop,
                                       self.mask.start, self.mask.stop):
                if self.mask.check_interaction(pol):
                    self.mask.react(pol)

        # Now recover elements
        for element in self.elements:
            if self.segments_intersect(pol.start, pol.stop,
                                       element.start, element.stop):
                element.cover()
                if element.check_interaction(pol):
                    # Resolve reactions between pol and element (e.g.,
                    # terminators)
                    element.react(pol)
            if self.mask != False:
                if self.segments_intersect(self.mask.start, self.mask.stop,
                                           element.start, element.stop):
                    element.cover()
            if element.old_covered == 0 and element.covered > 0 and element.type != "terminator":
                self.block_signal.fire(element.name)
                print(element.name, "covered!")
                element.old_covered = 1
            # Check for just-uncovered elements
            if element.old_covered >= 1 and element.covered == 0 and element.type != "terminator":
                self.promoter_signal.fire(element.name)
                print(element.name, "uncovered!")
                # Uncover element again in order to reset covering history
                # and avoid re-triggering an uncovering event.
                element.uncover()

    def execute(self):
        """
        Select a polymerase to move next and deal with terminations.
        """
        alpha_list = []

        for pol in self.polymerases:
            alpha_list.append(pol.speed)

        # Randomly select next polymerase to move, weighted by propensity
        pol = random.choices(self.polymerases, weights = alpha_list)[0]

        self.move_polymerase(pol)

        # Is the polymerase still attached?
        if pol.attached == False:
            self.termination_signal.fire(pol.last_gene, pol.name)
            self.polymerases.remove(pol)


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
        if self.mask != False:
            for i in range(self.mask.start - 1, self.mask.stop - 1):
                feature_locs[i] = "X"
        for feature in self.polymerases:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = "P"
        for feature in self.elements:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = feature.covered
        out_string = "\nfeatures: \n" + ''.join(map(str, feature_locs))
        return out_string

class Genome(Polymer):
    """
    Track polymerases on DNA, deal with collisions, promoters, terminators, and
    constructing transcripts.
    """
    def __init__(self, name, length, elements, transcript_template, mask):
        """
        :param name: name of this genome
        :param length: length of genome (do we still need this?)
        :param elements: DNA elements (promoters, terminators)
        :param transcript_template: list of parameters for all possible
            transcripts produced by this genome (i.e. the largest possible
            polycistronic transcript)
        """
        super().__init__(name, length, elements)
        self.transcript_template = transcript_template
        self.mask = mask
        for element in self.elements:
            if self.segments_intersect(element.start, element.stop,
                                       self.mask.start, self.mask.stop):
                element.cover()

    def execute(self):
        """
        Select next polymerase to move and deal with terminations.
        """
        alpha_list = []

        for pol in self.polymerases:
            alpha_list.append(pol.speed)

        # Randomly select next polymerase to move, weighted by propensity
        pol = random.choices(self.polymerases, weights = alpha_list)[0]

        self.move_polymerase(pol)

        if pol.attached == False:
            # Handle termination
            polymer, species = self.build_transcript(pol.bound, pol.stop)
            self.termination_signal.fire(polymer, pol.name, species)
            self.polymerases.remove(pol)

    def build_transcript(self, start, stop):
        """
        Build a transcript object corresponding to start and stop positions
        within this genome.

        :param start: start position of transcript within genome
        :param stop: stop position of transcript within genome

        :returns: polymer object, list of species that need to be added to
            species-level pool
        """
        species = []
        elements = []
        for element in self.transcript_template:
            if element["start"] >= start and element["stop"] <= stop:
                # Is this element within the start and stop sites?
                rbs = Promoter("rbs",
                               element["start"]-element["rbs"],
                               element["start"],
                               ["ribosome"])
                stop_site = Terminator("tstop",
                                       element["stop"],
                                       element["stop"],
                                       ["ribosome"])
                stop_site.gene = element["name"]
                elements.append(rbs)
                elements.append(stop_site)
                species.append("rbs")
            # build transcript
            polymer = Polymer("rna", 150, elements)
        return polymer, species

class Transcript(Polymer):
    """
    An mRNA transcript. Tracks ribosomes and protein production.
    """
    def __init__():
        pass
