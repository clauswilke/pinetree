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
    def __init__(self, name, length, elements, mask):
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
        self.mask = mask

        # Cover masked elements
        for element in self.elements:
            if self.segments_intersect(element.start, element.stop,
                                       self.mask.start, self.mask.stop):
                element.cover()

    def bind_polymerase(self, pol, promoter):
        """
        Bind a `Polymerase` object to the polymer.

        :param pol: `Polymerase` object.
        """

        found = False

        element_choices = []

        for element in self.elements:
            if element.name == promoter and element.covered == 0:
                element_choices.append(element)
                found = True

        element = random.choices(element_choices)[0]

        pol.start = element.start
        pol.stop = element.start + pol.footprint
        element.cover()
        element.save_state()
        self.polymerases.append(pol)
        assert(found == True)

    def count_uncovered(self, species):
        count = 0
        for element in self.elements:
            if element.name == species and element.covered == 0:
                count += 1
        return count

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

    def choose_polymerase(self):
        alpha_list = []

        for pol in self.polymerases:
            alpha_list.append(pol.speed)

        # Randomly select next polymerase to move, weighted by propensity
        pol = random.choices(self.polymerases, weights = alpha_list)[0]

        return pol

    def resolve_collisions(self, pol):
        collision = False
        for other_pol in self.polymerases:
            if pol == other_pol:
                continue
            if self.segments_intersect(pol.start, pol.stop,
                                       other_pol.start, other_pol.stop):
                if other_pol.check_interaction(pol):
                    other_pol.react(pol)
                    collision = True
        return collision

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
                element.save_state()
                element.uncover()
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.save_state()
                element.uncover()

        # Move polymerase
        pol.move()

        # First resolve any collisions between polymerases
        collision = self.resolve_collisions(pol)

        if collision == False:
            pol.move_signal.fire()

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
            # Re-cover masked elements
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.cover()
            # Check for newly-covered elements
            if element.was_covered() and element.type != "terminator":
                self.block_signal.fire(element.name)
                element.save_state()
            # Check for just-uncovered elements
            if element.was_uncovered() and element.type != "terminator":
                self.promoter_signal.fire(element.name)
                # Uncover element again in order to reset covering history
                # and avoid re-triggering an uncovering event.
                element.save_state()

    def uncover_elements(self):
        pass

    def execute(self):
        """
        Select a polymerase to move next and deal with terminations.
        """

        pol = self.choose_polymerase()

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
        feature_locs = ["o"]*self.length
        for i in range(self.mask.start - 1, self.mask.stop - 1):
            feature_locs[i] = "x"
        for feature in self.polymerases:
            for i in range(feature.start - 1, min(self.length, feature.stop - 1)):
                feature_locs[i] = "P"
        for feature in self.elements:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = feature.covered
        out_string = "\n"+self.name+": \n" + ''.join(map(str, feature_locs)) + "\n"
        print(self.elements)
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
        super().__init__(name, length, elements, mask)
        self.transcript_template = transcript_template
        self.transcript_signal = Signal()

    def bind_polymerase(self, pol, promoter):
        found = False

        for element in self.elements:
            if element.name == promoter and element.covered == 0:
                pol.start = element.start
                pol.stop = element.start + pol.footprint
                element.cover()
                element.save_state()
                self.polymerases.append(pol)
                found = True
                break

        assert(found == True)

        transcript, species = self.build_transcript(pol.start, self.length)
        pol.move_signal.connect(transcript.uncover_elements)
        self.transcript_signal.fire(transcript, species)

    def execute(self):
        """
        Select next polymerase to move and deal with terminations.
        """
        pol = self.choose_polymerase()

        self.move_polymerase(pol)

        if pol.attached == False:
            # Handle termination
            # polymer, species = self.build_transcript(pol.bound, pol.stop)
            self.termination_signal.fire(pol, pol.name, [])
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
                               element["start"]+element["rbs"],
                               element["start"],
                               ["ribosome"])
                stop_site = Terminator("tstop",
                                       element["stop"]-1,
                                       element["stop"],
                                       ["ribosome"])
                stop_site.gene = element["name"]
                elements.append(rbs)
                elements.append(stop_site)
                species.append("rbs")
        # build transcript
        polymer = Transcript("rna",
                          self.length,
                          elements,
                          TranscriptMask("mask", 0, self.length, ["ribosome"]))
        return polymer, species

class Transcript(Polymer):
    """
    An mRNA transcript. Tracks ribosomes and protein production.
    """
    def __init__(self, name, length, elements, mask):
        super().__init__(name, length, elements, mask)

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
                element.save_state()
                element.uncover()
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.save_state()
                element.uncover()

        # Move polymerase
        pol.move()

        # First resolve any collisions between polymerases
        collision = self.resolve_collisions(pol)

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
            # Re-cover masked elements
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.cover()
            # Check for newly-covered elements
            if element.was_covered() and element.type != "terminator":
                self.block_signal.fire(element.name)
                element.save_state()
            # Check for just-uncovered elements
            if element.was_uncovered() and element.type != "terminator":
                self.promoter_signal.fire(element.name)
                # Uncover element again in order to reset covering history
                # and avoid re-triggering an uncovering event.
                element.save_state()

    def uncover_elements(self):
        for element in self.elements:
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.save_state()
                element.uncover()

        self.mask.start += 1

        for element in self.elements:
            # Re-cover masked elements
            if self.segments_intersect(self.mask.start, self.mask.stop,
                                       element.start, element.stop):
                element.cover()
            # Check for just-uncovered elements
            if element.was_uncovered() and element.type != "terminator":
                self.promoter_signal.fire(element.name)
                # Uncover element again in order to reset covering history
                # and avoid re-triggering an uncovering event.
                element.save_state()
