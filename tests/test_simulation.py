#! /usr/bin/env python3

import unittest
from pysinthe.simulation import Simulation, Reaction, SpeciesReaction, Bridge


class TestReactionMethods(unittest.TestCase):

    def test_init(self):
        reaction = Reaction()
        self.assertEqual(reaction.index, -1)
        self.assertIsNone(reaction.calculate_propensity())
        self.assertIsNone(reaction.execute())


class TestSpeciesReaction(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation()
        self.reaction = SpeciesReaction(self.sim,
                                        1000,
                                        ["reactant1", "reactant2"],
                                        ["product1", "product2"])
        self.sim.register_reaction(self.reaction)

    def test_init(self):
        with self.assertRaises(RuntimeError):
            SpeciesReaction(self.sim,
                            1000,
                            ["reactant1", "reactant2", "reactant3"],
                            ["product1", "product2"])

        self.assertTrue(
            self.reaction in self.sim.reactant_bind_map["reactant1"]
        )
        self.assertTrue(
            self.reaction in self.sim.reactant_bind_map["reactant2"]
        )
        self.assertTrue(self.reaction in self.sim.reactant_bind_map["product1"])
        self.assertTrue(self.reaction in self.sim.reactant_bind_map["product2"])

    def test_calculate_propensity(self):
        self.sim.increment_reactant("reactant1", 2)
        self.sim.increment_reactant("reactant2", 3)
        prop = self.reaction.rate_constant * 2 * 3
        self.assertEqual(prop, self.reaction.calculate_propensity())

    def test_execute(self):
        self.setUp()
        self.sim.increment_reactant("reactant1", 2)
        self.sim.increment_reactant("reactant2", 3)
        self.reaction.execute()
        self.assertEqual(self.sim.reactants["reactant1"], 1)
        self.assertEqual(self.sim.reactants["reactant2"], 2)
        self.assertEqual(self.sim.reactants["product1"], 1)
        self.assertEqual(self.sim.reactants["product2"], 1)


class TestBridge(unittest.TestCase):

    def setUp(self):
        pass
