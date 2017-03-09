#! /usr/bin/env python3

import unittest
from unittest.mock import patch, MagicMock, Mock

import pysinthe
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


class TestSimulation(unittest.TestCase):

    def setUp(self):
        self.sim = Simulation()


    def test_register_reaction(self):
        pass

    def test_initialize_propensity(self):
        pass

    def test_update_propensity(self):
        pass

    def test_execute(self):
        pass

    def test_register_transcript(self):
        pass

    def test_terminate_transcription(self):
        pass

    def test_terminate_translation(self):
        pass

    def test_count_termination(self):
        pass