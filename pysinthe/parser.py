#! /usr/bin/env python3

import random

import yaml
from voluptuous import Schema, Optional, Any, All, Range, Length, Coerce

from .core import Promoter, Terminator, Polymerase
from .core import Genome, Mask
from .core import Reaction, Simulation, SpeciesReaction, Bind, SpeciesTracker
from .core import seed

class Parser:
    """
    A class to parse and validate simulation input parameters.
    """
    def __init__(self, myfile):
        """
        :param params: dictionary produced from parsing YAML file
        :param simulation: simulation object to populate
        """
        self.params = yaml.safe_load(myfile)
        self.simulation = Simulation()
        self.tracker = SpeciesTracker.get_instance()
        self.construct_simulation()

    def construct_simulation(self):
        """
        Construct a simulation according to set of YAML input parameters.
        """
        # validate parameters
        self._validate_schema(self.params)

        self.simulation.stop_time = self.params["simulation"]["runtime"]
        self.simulation.time_step = self.params["simulation"]["time_step"]
        # self.simulation.debug = self.params["simulation"]["debug"]
        self.tracker.cell_volume = float(
            self.params["simulation"]["cell_volume"]
        )

        # Set seed
        if "seed" in self.params["simulation"]:
            seed(self.params["simulation"]["seed"])

        # Parse reactions if present
        if "reactions" in self.params:
            self._parse_reactions(self.params["reactions"])

        genome, elements, mask = self._construct_genome(self.params["genome"],
                                        self.params["elements"])

        # Add species-level polymerase counts and construct list of partners
        # that this polymerase interacts with (promoters, terminators, etc.)
        self._parse_polymerases(self.params["polymerases"])

        # Add binding reaction for each promoter-polymerase interaction pair
        self._parse_binding_reactions(self.params["elements"],
                                      self.params["polymerases"])

        # Add species level reactants for elements that are not masked
        for element in elements:
            if element.type == "promoter" and element.stop < mask.start:
                self.tracker.increment(element.name, 1)
            elif element.type == "promoter":
                self.tracker.increment(element.name, 0)

        # Register genome
        self.simulation.register_genome(genome)

        self.tracker.increment("rbs", 0)
        self.tracker.increment(
            self.params["ribosomes"][0]["name"],
            int(self.params["ribosomes"][0]["copy_number"])
        )
        new_ribo = Polymerase(self.params["ribosomes"][0]["name"],
                     self.params["ribosomes"][0]["footprint"], # (10) footprint
                     self.params["ribosomes"][0]["speed"]) # (30) speed
        # Transcript-ribosome binding reaction
        reaction = Bind(float(self.params["ribosomes"][0]["binding_constant"]),
                        "rbs",
                        new_ribo)
        self.simulation.register_reaction(reaction)

        if "species" in self.params:
            self._parse_species(self.params["species"])


    def _validate_schema(self, params):
        """
        Validate parameter dictionary according to schema.
        """
        schema = Schema({
            'simulation': {Optional('seed'): All(int, Range(min=0)),
                           'runtime': All(Coerce(float), Range(min=0)),
                           'time_step': All(int, Range(min=0)),
                           'cell_volume': All(Coerce(float), Range(min=0)),
                           Optional('debug', default=False): bool},
            'genome': {'name': All(str, Length(min=1)),
                       'copy_number': All(int, Range(min=1)),
                       Optional('length'): All(int, Range(min=1)),
                       Optional('entered'): All(int, Range(min=1)),
                       Optional('mask_interactions'): list},
            'polymerases': All([{'name': All(str, Length(min=1)),
                                 'copy_number': All(int, Range(min=0)),
                                 'speed': All(int, Range(min=0)),
                                 'footprint': All(int, Range(min=1))}],
                               Length(min=1)),
            'ribosomes': All([{'name': All(str, Length(min=1)),
                               'copy_number': All(int, Range(min=0)),
                               'speed': All(int, Range(min=0)),
                               'footprint': All(int, Range(min=1)),
                               'binding_constant': All(Coerce(float),
                                                       Range(min=0))
                               }],
                             Length(min=1)),
            'elements': All([{'name': All(str, Length(min=1)),
                              'start': All(int, Range(min=0)),
                              'stop': All(int, Range(min=0)),
                              'type': Any('promoter',
                                          'terminator',
                                          'transcript'),
                              'interactions': dict,
                              'rbs': int}],
                            Length(min=1)),
            Optional('reactions'): list,
            Optional('species'): list
        }, required=True)
        schema(params)

    def _parse_species(self, species_params):
        for species in species_params:
            self.tracker.increment(
                species["name"],
                int(species["copy_number"])
            )

    def _parse_reactions(self, reaction_params):
        for reaction in reaction_params:
            new_reaction = SpeciesReaction(float(reaction["propensity"]),
                                           reaction["reactants"],
                                           reaction["products"])
            self.simulation.register_reaction(new_reaction)

    def _construct_genome(self, genome_params, element_params):
        """
        Build a genome.

        :param genome_params: genome parameters from YAML
        :param element_params: element parameters from YAML

        :return: genome object
        """
        # Construct list of DNA elements and a transcript template that
        # includes RBS's and stop sites
        dna_elements = []
        transcript_template = []
        last_position = 0
        old_transcript_stop = 0
        for element in element_params:
            if element["type"] == "promoter":
                # Add 1 because genome coordinates are inclusive
                new_element = Promoter(element["name"],
                                       element["start"],
                                       element["stop"],
                                       list(element["interactions"].keys()))
            elif element["type"] == "terminator":
                interactions = dict()
                for key, value in element["interactions"].items():
                    interactions[key] = value["efficiency"]
                new_element = Terminator(element["name"],
                                         element["start"],
                                         element["stop"],
                                         list(element["interactions"].keys()), 
                                         interactions)
            elif element["type"] == "transcript":
                old_transcript_stop = element["stop"]
                element["reading_frame"] = element["start"] % 3
                transcript_template.append(element)
                new_element = False
            if new_element is not False:
                dna_elements.append(new_element)
            last_position = element["stop"]

        # sort elements based on start position
        dna_elements.sort(key=lambda x: x.start)

        if "length" in genome_params:
            genome_length = genome_params["length"]
        else:
            genome_length = last_position

        # Build genome
        if "entered" in genome_params:
            if "mask_interactions" not in genome_params:
                raise ValueError(
                    "'mask_interactions' required if only a portion of the "
                    "genome has entered the cell."
                )
            genome_mask = Mask("mask",
                               genome_params["entered"],
                               genome_length,
                               genome_params["mask_interactions"])
        else:
            genome_mask = Mask("mask", genome_length + 1, genome_length, [])
        
        transcript_elems = self._build_transcript_template(transcript_template)

        genome = Genome(self.params["genome"]["name"],
                        genome_length,
                        dna_elements,
                        transcript_elems,
                        genome_mask)

        return genome, dna_elements, genome_mask
    
    def _build_transcript_template(self, template):
        elements = []
        for element in template:
            # Is this element within the start and stop sites?
            rbs = Promoter("rbs",
                            element["start"]+element["rbs"],
                            element["start"],
                            ["ribosome"])
            rbs.gene = element["name"]
            elements.append(rbs)
            stop_site = Terminator("tstop",
                                    element["stop"]-1,
                                    element["stop"],
                                    ["ribosome"],
                                    {"ribosome": 1.0})
            stop_site.reading_frame = element["start"] % 3
            stop_site.gene = element["name"]
            elements.append(stop_site)
        # Sort elements according to start position
        elements.sort(key=lambda x: x.start)
        return elements

    def _parse_polymerases(self, pol_params):
        """
        Add polymerases to species-level reactant count.

        :param pol_params: polymerase parameters
        """
        for pol in pol_params:
            # add polymerase as a reactant
            self.tracker.increment(pol["name"],
                                                      pol["copy_number"])

    def _parse_binding_reactions(self, element_params, pol_params):
        """
        Add binding reaction for each promoter-polymerase interaction pair

        :param element_params: element parameters
        :param pol_params: polymerase parameters
        """
        # Add binding reaction for each promoter-polymerase interaction pair
        seen = list()
        for element in element_params:
            if element["type"] == "promoter":
                self.tracker.increment(element["name"], 0)
                for partner, constant in element["interactions"].items():
                    binding_constant = constant["binding_constant"]
                    for pol in pol_params:
                        if pol["name"] == partner:
                            new_pol = Polymerase(partner,
                                        pol["footprint"],  # (10) footprint
                                        pol["speed"])
                            reaction = Bind(float(binding_constant),
                                            element["name"],
                                            new_pol)
                            # if (partner, element["name"]) not in seen:
                            self.simulation.register_reaction(reaction)
                            # seen.append((partner, element["name"]))
