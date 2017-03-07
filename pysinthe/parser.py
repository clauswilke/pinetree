#! /usr/bin/env python3

import random

import yaml
from voluptuous import Schema, Optional, Any, All, Range, Length

from .feature import Promoter, Terminator
from .polymer import Genome, Mask
from .simulation import Simulation, SpeciesReaction, Bind


class Parser:

    def __init__(self, myfile):
        self.params = yaml.safe_load(myfile)
        self.simulation = Simulation()
        self.construct_simulation()

    def construct_simulation(self):

        self._validate_schema(self.params)

        self.simulation.runtime = self.params["simulation"]["runtime"]
        self.simulation.time_step = self.params["simulation"]["time_step"]

        # Set seed
        random.seed(self.params["simulation"]["seed"])

        if "reactions" in self.params.keys():
            for reaction in self.params["reactions"]:
                new_reaction = SpeciesReaction(self.simulation,
                                               reaction["propensity"],
                                               reaction["reactants"],
                                               reaction["products"])
                self.simulation.register_reaction(new_reaction)

        # Construct list of DNA elements and a transcript template that includes
        # RBS's and stop sites
        dna_elements = []
        transcript_template = []
        position = 0
        for element in self.params["elements"]:
            if element["type"] == "promoter":
                length = element["length"]
                if element["length"] < 10:
                    # if promoter is smaller than RNApol footprint, expand out
                    # promoter length to avoid polymerases getting stuck
                    length = 10
                new_element = Promoter(element["name"],
                                       position,
                                       position + length,
                                       element["interactions"].keys())
            elif element["type"] == "terminator":
                new_element = Terminator(element["name"],
                                         position + element["length"] - 1,
                                         position + element["length"],
                                         element["interactions"])
            elif element["type"] == "transcript":
                transcript_template.append(element)
                position -= element["rbs"]
                new_element = False
            element["start"] = position
            element["stop"] = position + element["length"]
            if new_element != False:
                dna_elements.append(new_element)
            position += element["length"]

        # Build genome
        if "entered" in self.params["genome"]:
            genome_mask = Mask("mask", self.params["genome"]["entered"], position, ["rnapol", "ecolipol"])
        else:
            genome_mask = Mask("mask", position, position, [])
        genome = Genome(self.params["genome"]["name"],
                        position, dna_elements, transcript_template, genome_mask)


        # Add species-level polymerase counts and construct list of partners
        # that this polymerase interacts with (promoters, terminators, etc.)
        for pol in self.params["polymerases"]:
            # add polymerase as a reactant
            self.simulation.increment_reactant(pol["name"], pol["copy_number"])

        # Add binding reaction for each promoter-polymerase interaction pair
        for element in self.params["elements"]:
            if element["type"] == "promoter":
                self.simulation.increment_reactant(element["name"], 0)
                for partner, constant in element["interactions"].items():
                    binding_constant = constant["binding_constant"]
                    for pol in self.params["polymerases"]:
                        if pol["name"] == partner:
                            pol_args = [partner,
                                        element["start"],
                                        10,                 # footprint
                                        pol["speed"]
                                       ]
                            reaction = Bind(self.simulation,
                                            float(binding_constant),
                                            element["name"],
                                            pol_args)
                            self.simulation.register_reaction(reaction)
        # Add species level reactants for elements that are not masked
        for element in dna_elements:
            if element.type == "promoter" and element.stop < genome.mask.start:
                self.simulation.increment_reactant(element.name, 1)
            elif element.type == "promoter":
                self.simulation.increment_reactant(element.name, 0)

        # Register genome
        genome.termination_signal.connect(self.simulation.terminate_transcription)
        genome.promoter_signal.connect(self.simulation.free_promoter)
        genome.block_signal.connect(self.simulation.block_promoter)
        genome.transcript_signal.connect(self.simulation.register_transcript)
        self.simulation.register_polymer(genome)

        self.simulation.increment_reactant("rbs", 0)
        ribo_args = ["ribosome", 0, 10, #footprint
                     30]
        # Transcript-ribosome binding reaction
        reaction = Bind(self.simulation, float(1e7),
                        "rbs",
                        ribo_args)
        self.simulation.register_reaction(reaction)

        # Add species-level ribosomes
        self.simulation.increment_reactant("ribosome", self.params["simulation"]["ribosomes"])

        self.simulation.initialize_propensity()

    def _validate_schema(self, params):
        schema = Schema({
            'simulation': {Optional('seed'): All(int, Range(min=0)),
                           'runtime': All(int, Range(min=0)),
                           'time_step': All(int, Range(min=0)),
                           'ribosomes': All(int, Range(min=0)),
                           Optional('debug', default=False): bool},
            'genome': {'name': All(str, Length(min=1)),
                       'copy_number': All(int, Range(min=1)),
                       Optional('entered'): All(int, Range(min=1))},
            'polymerases': All([{'name': All(str, Length(min=1)),
                                 'copy_number': All(int, Range(min=0)),
                                 'speed': All(int, Range(min=0))}],
                               Length(min=1)),
            'elements': All([{'name': All(str, Length(min=1)),
                              'length': All(int, Range(min=0)),
                              'type': Any('promoter',
                                          'terminator',
                                          'transcript'),
                              'interactions': dict,
                              'rbs': int}],
                            Length(min=1)),
            Optional('reactions'): list
        }, required=True)
        schema(params)
