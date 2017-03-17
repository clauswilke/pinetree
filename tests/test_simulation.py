#! /usr/bin/env python3

import unittest
import random
from unittest.mock import patch, MagicMock

from pysinthe.simulation import (Simulation, Reaction, SpeciesReaction, Bind,
                                 Bridge, SpeciesTracker)


class TestSpeciesTrackerMethods(unittest.TestCase):

    @patch("pysinthe.simulation.Signal")
    def test_increment_species(self, mock_signal):
        tracker = SpeciesTracker()
        tracker.increment_species("reactant1", 1)
        self.assertEqual(tracker.species["reactant1"], 1)
        tracker.increment_species("reactant2", 1)
        self.assertEqual(tracker.species["reactant2"], 1)
        tracker.increment_species("reactant1", 2)
        self.assertEqual(tracker.species["reactant1"], 3)
        self.assertFalse(mock_signal.return_value.fire.called)

        mock_reaction = MagicMock()
        mock_reaction.index = 3
        tracker.add_reaction("reactant2", mock_reaction)
        tracker.increment_species("reactant2", 1)
        mock_signal.return_value.fire.assert_called_with(3)

    def test_add_reaction(self):
        tracker = SpeciesTracker()
        mock_reaction = MagicMock()
        tracker.add_reaction("reaction3", mock_reaction)
        self.assertTrue(
            mock_reaction in tracker.species_reaction_map["reaction3"]
        )
        tracker.add_reaction("reaction4", mock_reaction)
        self.assertTrue(
            mock_reaction in tracker.species_reaction_map["reaction3"]
        )
        self.assertTrue(
            mock_reaction in tracker.species_reaction_map["reaction4"]
        )

    def test_add_polymer(self):
        tracker = SpeciesTracker()
        mock_polymer = MagicMock()
        tracker.add_polymer("promoter1", mock_polymer)
        self.assertTrue(
            mock_polymer in tracker.promoter_polymer_map["promoter1"]
        )
        tracker.add_polymer("promoter2", mock_polymer)
        self.assertTrue(
            mock_polymer in tracker.promoter_polymer_map["promoter2"]
        )
        self.assertTrue(
            mock_polymer in tracker.promoter_polymer_map["promoter1"]
        )


class TestReactionMethods(unittest.TestCase):

    def test_init(self):
        reaction = Reaction()
        self.assertEqual(reaction.index, -1)
        self.assertIsNone(reaction.calculate_propensity())
        self.assertIsNone(reaction.execute())


class TestSpeciesReactionMethods(unittest.TestCase):

    def setUp(self):
        self.tracker = SpeciesTracker()
        self.reaction = SpeciesReaction(self.tracker,
                                        1000,
                                        ["reactant1", "reactant2"],
                                        ["product1", "product2"])

    def test_init(self):
        with self.assertRaises(RuntimeError):
            SpeciesReaction(self.tracker,
                            1000,
                            ["reactant1", "reactant2", "reactant3"],
                            ["product1", "product2"])

        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["reactant1"]
        )
        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["reactant2"]
        )
        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["product1"]
        )
        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["product2"]
        )

    def test_calculate_propensity(self):
        self.tracker.increment_species("reactant1", 2)
        self.tracker.increment_species("reactant2", 3)
        prop = self.reaction.rate_constant * 2 * 3
        self.assertEqual(prop, self.reaction.calculate_propensity())

    def test_execute(self):
        self.setUp()
        self.tracker.increment_species("reactant1", 2)
        self.tracker.increment_species("reactant2", 3)
        self.reaction.execute()
        self.assertEqual(self.tracker.species["reactant1"], 1)
        self.assertEqual(self.tracker.species["reactant2"], 2)
        self.assertEqual(self.tracker.species["product1"], 1)
        self.assertEqual(self.tracker.species["product2"], 1)


class TestBindMethods(unittest.TestCase):

    def setUp(self):
        self.tracker = SpeciesTracker()
        Reaction._CELL_VOLUME = float(8e-15)
        self.reaction = Bind(self.tracker,
                             1000,
                             "promoter1",
                             ["pol1", 10, 20])

    def test_init(self):
        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["promoter1"]
        )
        self.assertTrue(
            self.reaction in self.tracker.species_reaction_map["pol1"]
        )

    def test_calculate_propensity(self):
        self.tracker.increment_species("promoter1", 2)
        self.tracker.increment_species("pol1", 3)
        prop = self.reaction.rate_constant * 2 * 3
        self.assertEqual(prop, self.reaction.calculate_propensity())

    @patch("pysinthe.simulation.Polymerase")
    def test_execute(self, mock_polymerase):
        self.setUp()
        self.tracker.increment_species("promoter1", 2)
        self.tracker.increment_species("pol1", 3)
        mock_polymer = MagicMock()
        mock_polymer.count_uncovered.return_value = 3
        self.tracker.promoter_polymer_map["promoter1"] = [mock_polymer]
        self.reaction.execute()
        self.assertTrue(mock_polymerase.called)
        self.assertTrue(mock_polymer.bind_polymerase.called)
        self.assertEqual(self.tracker.species["promoter1"], 1)
        self.assertEqual(self.tracker.species["pol1"], 2)


class TestBridgeMethods(unittest.TestCase):

    @patch("pysinthe.simulation.Signal")
    def setUp(self, mock_signal):
        self.mock_signal = mock_signal
        self.mock_polymer = MagicMock()
        self.reaction = Bridge(self.mock_polymer)

    def test_calculate_propensity(self):
        self.reaction.calculate_propensity()
        self.assertTrue(self.mock_polymer.calculate_propensity.called)

    def test_execute(self):
        self.reaction.execute()
        self.assertTrue(self.mock_polymer.execute.called)

    def test_update(self):
        self.reaction.update()
        # index will not yet have been set
        self.mock_signal.return_value.fire.assert_called_with(-1)


class TestSimulationMethods(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation()

    def test_register_reaction(self):
        # Add one reaction
        reaction1 = MagicMock()
        reaction1.calculate_propensity.return_value = 0.5
        self.sim.register_reaction(reaction1)
        self.assertEqual(self.sim.alpha_list[0], 0.5)
        self.assertEqual(self.sim.reactions.index(reaction1), 0)
        # Add a second reaction
        reaction2 = MagicMock()
        reaction2.calculate_propensity.return_value = 0.8
        self.sim.register_reaction(reaction2)
        self.assertEqual(self.sim.alpha_list[1], 0.8)
        self.assertEqual(self.sim.reactions.index(reaction2), 1)
        # Make sure that alpha_list is really acting as a cache
        reaction1.calculate_propensity.return_value = 0.1
        self.assertNotEqual(self.sim.alpha_list[0], 0.1)

    def test_register_polymer(self):
        self.setUp()
        polymer = MagicMock()
        mock_promoter = MagicMock()
        mock_promoter.type = "promoter"
        mock_promoter.name = "mypromoter"
        polymer.elements = [mock_promoter]
        self.sim.register_polymer(polymer)
        # Make sure promoter/polymer was added to map
        self.assertEqual(self.sim.tracker.polymers_with_promoter("mypromoter"),
                         [polymer])
        # Make sure bridge reaction was added
        self.assertEqual(self.sim.reactions[0].polymer, polymer)

    def test_initialize_propensity(self):
        # add reactions to simulation as before
        self.setUp()
        self.test_register_reaction()
        # Make sure that cache has not been updated
        self.assertNotEqual(self.sim.alpha_list[0], 0.1)
        # Update cache
        self.sim.initialize_propensity()
        self.assertEqual(self.sim.alpha_list[0], 0.1)
        # reaction2 propensity should not have changed
        self.assertEqual(self.sim.alpha_list[1], 0.8)

    def test_update_propensity(self):
        # add reactions to simulation as before
        self.setUp()
        self.test_register_reaction()
        # Make sure that cache has not been updated
        self.assertNotEqual(self.sim.alpha_list[0], 0.1)
        # Update cache
        self.sim.update_propensity(0)
        self.assertEqual(self.sim.alpha_list[0], 0.1)
        # Update propensity in reaction object
        self.sim.reactions[1].calculate_propensity.return_value = 7
        self.assertNotEqual(self.sim.alpha_list[1], 7)
        self.sim.update_propensity(1)
        self.assertEqual(self.sim.alpha_list[1], 7)

    def test_execute(self):
        self.setUp()
        reaction1 = MagicMock()
        reaction1.calculate_propensity.return_value = 0.5
        self.sim.register_reaction(reaction1)
        reaction2 = MagicMock()
        reaction2.calculate_propensity.return_value = 0.8
        self.sim.register_reaction(reaction2)

        self.sim.initialize_propensity()
        random.seed(17)
        old_iteration = self.sim.iteration
        old_time = self.sim.time
        new_time = old_time + 0.5000911660616504
        self.sim.execute()
        self.assertEqual(self.sim.time, new_time)
        self.assertTrue(reaction2.execute.called)
        self.assertEqual(self.sim.iteration, old_iteration + 1)
        self.assertFalse(reaction1.execute.called)
        self.sim.execute()
        self.assertTrue(reaction1.execute.called)
