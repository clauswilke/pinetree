#! /usr/bin/env python3

import random

from eventsignal import Signal
from feature import Promoter, Terminator, TranscriptMask

class Polymer:
    """
    Track element objects, polymerase objects, and collisions on a single
    polymer. Move polymerase objects along the polymer. Handle logic for
    covering and uncovering of elements. This class contains the core of the
    single-moleculre tracking used for both genomes (transcription) and
    transcripts (translation).

    The terms polymer, polymerase, promoter, and terminator are all used
    generically in this class. Each term could refer to a different biological
    definition in the context of transcription and translation.

    * polymer: genome, transcript
    * polymerase: RNA polymerase, ribosome, any object that binds to polymer
    * promoter: promoter, ribosome binding site, any site on a polymer in which
        a protein (polymerase) can bind
    * terminator: terminator, stop codon, any site on polymer that ends
        polymerization

    """
    def __init__(self, name, length, elements, mask):
        """
        :param name: name of this polymer (should it be unique?)
        :param length: length of polymer (used purely for debugging, may not be
            needed)
        :param elements: all elements on this polymer, including promoters,
            terminators, etc.
        :param mask: mask object which determines which portions of the polymer
            are currently inaccessible
        """
        self.index = 0
        self.name = name
        self.length = length
        self.polymerases = []
        self.elements = elements
        self.termination_signal = Signal() # Fires on termination
        self.promoter_signal = Signal() # Fires when promoter is freed
        self.block_signal = Signal() # Fires when promoter is blocked
        self.propensity_signal = Signal() # Fires when propensity changes
        self.mask = mask
        self.prop_sum = 0
        self.uncovered = {} # Running count of free promoters

        # Cover masked elements
        for element in self.elements:
            if self.elements_intersect(element, self.mask):
                element.cover()
                if element.name not in self.uncovered:
                    self.uncovered[element.name] = 0
            else:
                if element.name not in self.uncovered:
                    self.uncovered[element.name] = 1
                else:
                    self.uncovered[element.name] += 1

    def bind_polymerase(self, pol, promoter):
        """
        Bind a polymerase object to the polymer. Randomly select an open
        promoter with which to bind and update the polymerases position to the
        position of that promoter.

        :param pol: polymerase object.
        :param promoter: the name of a promoter that pol will bind

        """
        found = False
        element_choices = []
        # Make list of free promoters that pol can bind
        for element in self.elements:
            if element.name == promoter and element.covered == 0:
                element_choices.append(element)
                found = True

        # Randomly select promoter
        element = random.choices(element_choices)[0]
        # Update polymerase coordinates
        pol.start = element.start
        pol.stop = element.start + pol.footprint
        # Cover promoter
        element.cover()
        element.save_state()
        self.uncovered[element.name] -= 1
        # Add polymerase to tracked-polymerases list
        # Polymerases are maintained in order, such that higher-index
        # polymerases have moved further along the DNA
        # This make collision detection very efficient
        self.insert_polymerase(pol)
        # self.polymerases.append(pol)
        self.prop_sum += pol.speed
        # Sanity check; this function should never be called if there are no
        # free promoters with which to bind
        self.propensity_signal.fire()
        assert found

    def insert_polymerase(self, pol):
        found_position = False
        insert_position = 0
        for index, old_pol in enumerate(self.polymerases):
            # Find the first polymerase that is
            insert_position = index
            if pol.start < old_pol.start:
                found_position = True
                break
        if found_position == False:
            # Check to see if we're actually just adding to the end of the list
            insert_position += 1
        self.polymerases.insert(insert_position, pol)

    def count_uncovered(self, species):
        """
        Count the number of free promoters that match name `species`.

        TODO: cache uncovered values

        :param species: name of promoter to count
        """
        return self.uncovered[species]

    def calculate_propensity(self):
        """
        Calculate the total propensity of all polymerase movement in this
        polymer.

        :returns: total propensity
        """
        return self.prop_sum

    def execute(self):
        """
        Select a polymerase to move next and deal with terminations.
        """
        # Randomly choose polymerase to move
        pol = self.choose_polymerase()
        self.move_polymerase(pol)

        # Is the polymerase still attached?
        if not pol.attached:
            self.terminate(pol)

    def terminate(self, pol):
        self.prop_sum -= pol.speed
        self.polymerases.remove(pol)

    def choose_polymerase(self):
        """
        Randomly select next polymerase to move, weighted by move propensity
        (i.e. speed)

        :returns: selected polymerase
        """
        # Construct list of movement propensities
        prop_list = [pol.speed for pol in self.polymerases]
        # Randomly select next polymerase to move, weighted by propensity
        if len(prop_list) == 0:
            print("There are no active polymerases.", str(self.prop_sum))
            print(self.polymerases)
            raise
        pol = random.choices(self.polymerases, weights=prop_list)[0]

        return pol

    def move_polymerase(self, pol):
        """
        Move polymerase and deal with collisions and covering/uncovering of
        elements.

        :param pol: polymerase to move
        """
        # Find which elements this polymerase (or mask) is covering and
        # temporarily uncover them
        self.uncover_elements(pol)

        # Move polymerase
        pol.move()

        # First resolve any collisions between polymerases
        collision = self.resolve_collisions(pol)

        # If no collisions occurred, it's safe to broadcast that polymerase
        # has moved
        if not collision:
            pol.move_signal.fire()

        # If this polymerase interacts with the mask, push the mask back and
        # uncover more elements on the polymer
        if self.elements_intersect(pol, self.mask):
            if self.mask.check_interaction(pol):
                self.mask.react(pol)

        # Now recover elements and check for changes in covered elements
        self.recover_elements(pol)

    def resolve_collisions(self, pol):
        """
        Resolve collisions between polymerases.

        :param pol: polymerase with which to check for collisions
        :returns: True if there was at least 1 collision, False otherwise
        """
        collision = False
        index = self.polymerases.index(pol)
        # We only need to check the polymerase ahead of this polymerase
        if index + 1 > len(self.polymerases) - 1:
            return collision
        if self.elements_intersect(pol,
                                   self.polymerases[index + 1]):
            if self.polymerases[index + 1].check_interaction(pol):
                self.polymerases[index + 1].react(pol)
                collision = True
        return collision

    def uncover_elements(self, pol):
        """
        Uncover all elements currently covered by `pol` or the polymer's mask.
        Save the state of each element to check for uncoverings later.

        :param pol: polymerase
        """
        for element in self.elements:
            if self.elements_intersect(pol, element):
                element.save_state()
                element.uncover()
            if self.elements_intersect(self.mask, element):
                element.save_state()
                element.uncover()

    def recover_elements(self, pol):
        """
        Recover elements covered by `pol` and the polymer mask, and fire
        appropriate signals for elements that have changed state.
        """
        for element in self.elements:
            if self.elements_intersect(pol, element):
                element.cover()
                if element.check_interaction(pol):
                    # Resolve reactions between pol and element (e.g.,
                    # terminators)
                    element.react(pol)
            # Re-cover masked elements
            if self.elements_intersect(self.mask, element):
                element.cover()
            # Check for newly-covered elements
            if element.was_covered() and element.type != "terminator":
                self.block_signal.fire(element.name)
                element.save_state()
                self.uncovered[element.name] -= 1
            # Check for just-uncovered elements
            if element.was_uncovered() and element.type != "terminator":
                self.promoter_signal.fire(element.name)
                # Save current state to avoid re-triggering an uncovering event.
                element.save_state()
                self.uncovered[element.name] += 1
            if element.was_uncovered() and element.type == "terminator":
                # Reset readthrough state of terminator
                self.uncovered[element.name] += 1
                element.readthrough = False


    def elements_intersect(self, element1, element2):
        """
        Do two line segments (e.g. `Polymerase` objects) overlap?

        :param element1: first element
        :param element2: second element
        :returns: True if elements intersect
        """
        return element1.stop >= element2.start and \
            element2.stop >= element1.start

    def __str__(self):
        """
        Convert `Polymer` object to string representation showing features and
        polymerases.
        """
        feature_locs = ["o"]*self.length
        for i in range(self.mask.start - 1, self.mask.stop - 1):
            feature_locs[i] = "x"
        for index, feature in enumerate(self.polymerases):
            for i in range(feature.start - 1,
                           min(self.length, feature.stop - 1)):
                feature_locs[i] = "P" + str(index)
        for feature in self.elements:
            for i in range(feature.start - 1, feature.stop - 1):
                feature_locs[i] = feature.covered
        out_string = "\n"+self.name+": \n" + ''.join(map(str, feature_locs)) + \
            "\n"
        return out_string

class Genome(Polymer):
    """
    Track polymerases on DNA, deal with collisions, promoters, terminators, and
    constructing transcripts. Inherits from Polymer. Unlike Polymer, Genome
    must construct a transcript upon promoter binding.
    """
    def __init__(self, name, length, elements, transcript_template, mask):
        """
        :param name: name of this genome
        :param length: length of genome (do we still need this?)
        :param elements: DNA elements (promoters, terminators)
        :param transcript_template: list of parameters for all possible
            transcripts produced by this genome (i.e. the largest possible
            polycistronic transcript)
        :param mask: polymer mask (i.e. portion of genome that has not yet
            entered the cell and remains inaccessible)
        """
        super().__init__(name, length, elements, mask)
        self.transcript_template = transcript_template
        self.transcript_signal = Signal() # fires upon transcript construction

    def bind_polymerase(self, pol, promoter):
        """
        Bind a polymerase to genome and construct new transcript.

        :param pol: polymerase to bind
        :param promoter: name of promoter to which this polymerase binds
        """
        # Bind polymerase just like in parent Polymer
        super().bind_polymerase(pol, promoter)
        # Construct transcript
        transcript = self.build_transcript(pol.start, self.length)
        # Connect polymerase movement signal to transcript, so that the
        # transcript knows when to expose new elements
        pol.move_signal.connect(transcript.shift_mask)
        pol.termination_signal.connect(transcript.release)
        # Fire new transcript signal
        self.transcript_signal.fire(transcript)

    def terminate(self, pol):
        self.prop_sum -= pol.speed
        self.termination_signal.fire(pol.name)
        self.propensity_signal.fire() # Update propensities
        self.polymerases.remove(pol)

    def build_transcript(self, start, stop):
        """
        Build a transcript object corresponding to start and stop positions
        within this genome.

        TODO: Find less janky way of constructing a transcript

        :param start: start position of transcript within genome
        :param stop: stop position of transcript within genome
        :returns: polymer object, list of species that need to be added to
            species-level pool
        """
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
                                       {"ribosome": {"efficiency": 1.0}})
                stop_site.gene = element["name"]
                elements.append(rbs)
                elements.append(stop_site)
        # build transcript
        polymer = Transcript("rna",
                             self.length,
                             elements,
                             TranscriptMask("mask", 0, self.length,
                                            ["ribosome"]))
        return polymer

class Transcript(Polymer):
    """
    An mRNA transcript. Tracks ribosomes and protein production. Only differs
    from Polymer in capability to receive signals from a moving polymerase and
    uncover the appropriate, "newly-synthesized" elements.
    """
    def __init__(self, name, length, elements, mask):
        super().__init__(name, length, elements, mask)

    def terminate(self, pol):
        self.prop_sum -= pol.speed
        self.termination_signal.fire(pol.last_gene, pol.name)
        self.propensity_signal.fire() # Update propensities
        self.polymerases.remove(pol)

    def release(self, stop):
        jump = stop - self.mask.start
        self.mask.start += jump

    def shift_mask(self):
        """
        Shift start of mask by 1 base-pair and check for uncovered elements.
        """
        index = -1
        for index, element in enumerate(self.elements):
            if self.elements_intersect(self.mask, element):
                element.save_state()
                element.uncover()
                break

        self.mask.start += 1

        if index == -1:
            return
        # Re-cover masked elements
        if self.elements_intersect(self.mask, self.elements[index]):
            self.elements[index].cover()
        # Check for just-uncovered elements
        if self.elements[index].was_uncovered() and self.elements[index].type != "terminator":
            self.promoter_signal.fire(self.elements[index].name)
            # Uncover element again in order to reset covering history
            # and avoid re-triggering an uncovering event.
            self.elements[index].save_state()
            self.uncovered[self.elements[index].name] += 1
