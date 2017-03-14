#! /usr/bin/env python3

import random
import math

from .feature import Polymerase
from .eventsignal import Signal
from .choices import weighted_choice

CELL_VOLUME = float(8e-16)


class SpeciesTracker():
    """
    Tracks species' copy numbers and maintains promoter-to-polymer and
    species-to-reaction maps to easily look up which polymers contain and given
    promoter and which reactions involve a given species. These maps are
    necessarry to cache propensities in the simulation.

    TODO: Move propensity cache into this class?
    """
    def __init__(self):
        self.species = {}  # species-level reactant counts
        self.promoter_polymer_map = {}  # Map of which polymers contain a given
        # promoter
        self.species_reaction_map = {}  # Map of which binding reaction
        # involves a given reactant
        self.propensity_signal = Signal()

    def increment_species(self, species_name, copy_number):
        """
        Change a species count by a given value (positive or negative).

        :param species_name: name of species to change count
        :param copy_number: increase current copy number by this number
        """
        if species_name in self.species:
            self.species[species_name] += copy_number
        else:
            self.species[species_name] = copy_number

        if copy_number == 0:
            return

        # Tell simulation to update propensity cache
        if species_name in self.species_reaction_map:
            for reaction in self.species_reaction_map[species_name]:
                self.propensity_signal.fire(reaction.index)

    def add_reaction(self, species_name, reaction):
        """
        Add a species-reaction pair to species-reaction map.

        :param species_name: name of species
        :param reaction: reaction object taht involves species
        """
        if species_name not in self.species_reaction_map:
            self.species_reaction_map[species_name] = [reaction]
        else:
            self.species_reaction_map[species_name].append(reaction)

    def add_polymer(self, promoter_name, polymer):
        """
        Add promoter-polymer pair to promoter-polymer map.

        :param promoter_name: name of promoter
        :param polymer: polymer that contains the named promoter
        """
        if promoter_name not in self.promoter_polymer_map:
            self.promoter_polymer_map[promoter_name] = [polymer]
        else:
            self.promoter_polymer_map[promoter_name].append(polymer)

    def polymers_with_promoter(self, promoter_name):
        """
        Get polymers that contain a given promoter.

        :param promoter_name: name of promoter

        :returns: list of polymer objects
        """
        return self.promoter_polymer_map[promoter_name]

    def reactions_with_species(self, species_name):
        """
        Get reactions that involve a given species.

        :param species_name: name of species

        :returns: list of reaction objects
        """
        return self.species_reaction_map[species_name]


class Reaction():
    """
    Generic class for a reaction. (Not currently used).
    """
    _AVAGADRO = float(6.0221409e+23)

    def __init__(self):
        self.index = -1  # Index inside propensity list, sim.alpha_list

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction
        """
        pass

    def execute(self):
        """
        Execute the reaction.
        """
        pass


class SpeciesReaction(Reaction):
    """
    Generic class for species-level reaction.
    """
    def __init__(self, tracker, rate_constant, reactants, products):
        super().__init__()
        self.tracker = tracker
        if len(reactants) > 2:
            raise RuntimeError("Simulation does not support reactions with "
                               "more than two reactant species.")
        if len(reactants) == 2:
            self.rate_constant = float(rate_constant)/(self._AVAGADRO*CELL_VOLUME)
        else:
            self.rate_constant = float(rate_constant)
        self.reactants = reactants
        self.products = products

        for reactant in self.reactants:
            self.tracker.add_reaction(reactant, self)
            self.tracker.increment_species(reactant, 0)

        for product in self.products:
            self.tracker.add_reaction(product, self)
            self.tracker.increment_species(product, 0)

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction
        """
        propensity = self.rate_constant
        for reactant in self.reactants:
            propensity *= self.tracker.species[reactant]

        return propensity

    def execute(self):
        """
        Execute the reaction.
        """
        for reactant in self.reactants:
            self.tracker.increment_species(reactant, -1)

        for product in self.products:
            self.tracker.increment_species(product, 1)


class Bind(Reaction):
    """
    Bind a polymerase to a polymer.
    """
    def __init__(self, tracker, rate_constant, promoter_name, pol_args):
        """

        TODO: Refactor so Bind inherits from SpeciesReaction

        :param sim: reference to simulation object in which this reaction occurs
        :param rate_constant: rate constant of reaction
        :param promoter_name: name of promoter involved in this reaction
        :param pol_args: list of arguments to pass to polymerase
            constructor upon execution of this reaction.
        """
        super().__init__()
        self.tracker = tracker
        self.polymerase = pol_args[0]
        self.pol_args = pol_args
        self.rate_constant = float(rate_constant)/(self._AVAGADRO*CELL_VOLUME)
        self.promoter_name = promoter_name

        self.tracker.add_reaction(self.promoter_name, self)
        self.tracker.add_reaction(self.polymerase, self)

    def calculate_propensity(self):
        """
        Calculate the propensity of this reaction.

        :returns: propensity of this reaction.
        """

        return self.rate_constant * self.tracker.species[self.polymerase]\
            * self.tracker.species[self.promoter_name]

    def execute(self):
        """
        Decrement reactants, choose a polymer, construct new polymerase, and
        bind the polymerase to the polymer.
        """
        # Find which polymer we should bind to
        weights = []
        for polymer in self.tracker.promoter_polymer_map[self.promoter_name]:
            # Weight by the number of unbound promoters in each polymer
            weights.append(polymer.count_uncovered(self.promoter_name))
        polymer = weighted_choice(
            self.tracker.promoter_polymer_map[self.promoter_name],
            weights=weights
        )
        # Construct new polymerase
        new_pol = Polymerase(*self.pol_args)
        polymer.bind_polymerase(new_pol, self.promoter_name)

        self.tracker.increment_species(self.promoter_name, -1)
        self.tracker.increment_species(self.polymerase, -1)

    def __str__(self):
        return self.promoter_name + "-" + self.polymerase


class Bridge(Reaction):
    """
    Encapsulate polymer so it can participate in species-level reaction
    processing.
    """
    def __init__(self, polymer):
        """
        :param polymer: polymer object that this reaction is encapsulating
        """
        super().__init__()
        self.polymer = polymer
        self.propensity_signal = Signal()  # Signal to fire when propensity
        # needs to be updated
        polymer.propensity_signal.connect(self.update)

    def calculate_propensity(self):
        """
        Retrieve the total propensity of all reactions that may occur within
        this polymer.

        :returns: total propensity of reactions within polymer
        """
        return self.polymer.calculate_propensity()

    def execute(self):
        """
        Execute reaction within polymer (e.g. typically moving a polymerase).
        """
        self.polymer.execute()

    def update(self):
        """
        Fires when polymer removes a polymerase from polymer and propensity
        needs to be updated.
        """
        self.propensity_signal.fire(self.index)

    def __str__(self):
        """
        Retrieve string representation of polymer.

        :returns: string representation of polymer encapsulated by this object
        """
        return str(self.polymer)


class Simulation:
    """
    Coordinate polymers and species-level reactions.
    """
    def __init__(self):
        self.time = 0  # simulation time
        self.runtime = 0
        self.time_step = 0
        self.debug = False
        self.terminations = {}  # track when polymerases terminate
        self.tracker = SpeciesTracker()
        self.reactions = []  # all reactions
        self.iteration = 0  # iteration counter
        self.alpha_list = []

        self.tracker.propensity_signal.connect(self.update_propensity)

    def run(self):
        print(self)
        old_time = 0
        while self.time < self.runtime:
            # Execute simulatio
            self.execute()
            if abs(self.time - old_time) > self.time_step:
                # Output data every ~time_step
                print(self)
                old_time += self.time_step
            elif self.debug is True:
                print(self)
                for pol in self.reactions:
                    print(pol)

    def register_reaction(self, reaction):
        """
        Add a SpeciesReaction object to the list of reactions.

        :param reaction: SpeciesReaction object
        """
        if reaction not in self.reactions:
            reaction.index = len(self.reactions)
            try:
                self.alpha_list.append(reaction.calculate_propensity())
            except KeyError:
                self.alpha_list.append(0.0)
            self.reactions.append(reaction)

    def register_polymer(self, polymer):
        """
        Add a polymer to the simulation.

        :param polymer: polymer object
        """
        # Add polymer to promoter-polymerase map
        for element in polymer.elements:
            if element.type == "promoter":
                self.tracker.add_polymer(element.name, polymer)

        # Encapsulate polymer in Bridge reaction and add to reaction list
        bridge = Bridge(polymer)
        bridge.propensity_signal.connect(self.update_propensity)
        self.register_reaction(bridge)

    def initialize_propensity(self):
        """
        Initialize all propensities before the start of the simulation.
        """
        for index, reaction in enumerate(self.reactions):
            self.alpha_list[index] = reaction.calculate_propensity()

        prop_sum = sum(self.alpha_list)

        if prop_sum == 0:
            raise ValueError("Initial global reaction propensity is 0. Make "
                             "sure that at least one promoter is exposed in "
                             "the genome.")

    def update_propensity(self, index):
        """
        Update a propensity of a reaction at a given index.
        """
        self.alpha_list[index] = float(
            self.reactions[index].calculate_propensity()
            )

    def execute(self):
        """
        Execute one cycle of reaction.
        """
        # Generate random number
        random_num = random.random()

        # Sum propensities
        alpha = sum(self.alpha_list)

        # Calculate tau, i.e. time until next reaction
        tau = (1/alpha)*math.log(1/random_num)
        self.time += tau

        # Randomly select next reaction to execute, weighted by propensities
        # print(self.reactions)
        next_reaction = weighted_choice(self.reactions,
                                        weights=self.alpha_list)
        next_reaction.execute()

        self.iteration += 1

    def free_promoter(self, species):
        """
        Increment promoter count at species level.
        """
        self.tracker.increment_species(species, 1)

    def block_promoter(self, species):
        """
        Decrement promoter count at species level.
        """
        self.tracker.increment_species(species, -1)

    def register_transcript(self, polymer):
        """
        Register a new transcript with the simulation.

        :param polymer: transcript object to be added to simulation
        :param reactants: proteins to be produced by this transcript
        """
        self.register_polymer(polymer)
        # Connect signals
        polymer.promoter_signal.connect(self.free_promoter)
        polymer.termination_signal.connect(self.terminate_translation)

    def terminate_transcription(self, species):
        """
        Terminate transcription.

        :param polymer: the newly constructed transcript
        :param species: name of the polymerase completing transcription
        :param reactants: list of reactants to add to species-level pool
            (usually RBSs)
        """
        self.tracker.increment_species(species, 1)
        self.count_termination("transcript")

    def terminate_translation(self, protein, species):
        """
        Terminate translation.

        :param name: name of the protein being produced
        :param species: name of species that just translated this protein
            (usually 'ribosome')
        """
        self.tracker.increment_species(species, 1)
        self.tracker.increment_species(protein, 1)
        self.count_termination(protein)

    def count_termination(self, name):
        """
        Record when a polymerase reaches a terminator.
        """
        name = name + "_total"
        if name not in self.terminations.keys():
            self.terminations[name] = 1
        else:
            self.terminations[name] += 1

    def __str__(self):
        """
        Convert `Simulation` to string representation showing names and
        termination times of `Polymerase` objects.
        """
        out_string = ""
        for name, count in self.terminations.items():
            out_string += str(self.iteration) + ", " + str(float(self.time)) + \
                ", " + name + ", " + str(count) + "\n"
        for name, count in self.tracker.species.items():
            out_string += str(self.iteration) + ", " + str(float(self.time)) + \
                ", " + name + ", " + str(count) + "\n"
        return out_string.strip()
