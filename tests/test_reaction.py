import unittest
import pinetree as pt


class TestSpeciesReactionMethods(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(ValueError):
            rxn = pt.SpeciesReaction(-1.0, 1.0, ["x"], ["y"])
        with self.assertRaises(ValueError):
            rxn = pt.SpeciesReaction(1.0, -1.0, ["x"], ["y"])
        with self.assertRaises(ValueError):
            rxn = pt.SpeciesReaction(1.0, 1.0, ["x", "z", "m"], ["y"])
        with self.assertRaises(ValueError):
            rxn = pt.SpeciesReaction(1.0, 1.0, [], [])

