#! /usr/bin/env python3

import unittest
from unittest.mock import patch, MagicMock, Mock

import pysinthe
from pysinthe.simulation import (Simulation, Reaction, SpeciesReaction, Bind,
                                 Bridge)


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


class TestBind(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation()
        self.reaction = Bind(self.sim,
                             1000,
                             "promoter1",
                             ["pol1", 10, 20])
        self.sim.register_reaction(self.reaction)

    def test_init(self):
        self.assertTrue(
            self.reaction in self.sim.reactant_bind_map["promoter1"]
        )
        self.assertTrue(
            self.reaction in self.sim.reactant_bind_map["pol1"]
        )

    def test_calculate_propensity(self):
        self.sim.increment_reactant("promoter1", 2)
        self.sim.increment_reactant("pol1", 3)
        prop = self.reaction.rate_constant * 2 * 3
        self.assertEqual(prop, self.reaction.calculate_propensity())

    @patch("pysinthe.polymer.Polymer")
    @patch("pysinthe.simulation.Polymerase")
    def test_execute(self, mock_polymerase, mock_polymer):
        self.setUp()
        self.sim.increment_reactant("promoter1", 2)
        self.sim.increment_reactant("pol1", 3)
        mock_polymer.count_uncovered = Mock(return_value=3)
        self.sim.promoter_polymer_map["promoter1"] = [mock_polymer]
        self.reaction.execute()
        self.assertTrue(mock_polymerase.called)
        self.assertTrue(mock_polymer.bind_polymerase.called)
        self.assertEqual(self.sim.reactants["promoter1"], 1)
        self.assertEqual(self.sim.reactants["pol1"], 2)


class TestBridge(unittest.TestCase):

    @patch("pysinthe.polymer.Polymer")
    def setUp(self, mock_polymer):
        self.mock_polymer = mock_polymer
        self.reaction = Bridge(self.mock_polymer)
        self.reaction.propensity_signal.connect(self.fired)
        self.has_fired = False

    def fired(self, index):
        self.has_fired = True

    def test_calculate_propensity(self):
        self.reaction.calculate_propensity()
        self.assertTrue(self.mock_polymer.calculate_propensity.called)

    def test_execute(self):
        self.reaction.execute()
        self.assertTrue(self.mock_polymer.execute.called)

    def test_update(self):
        self.reaction.update()
        self.assertTrue(self.has_fired)
