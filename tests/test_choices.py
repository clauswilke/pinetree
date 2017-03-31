import unittest
import random
# from pysinthe.choices import weighted_choice
import cppimport
module = cppimport.imp("pysinthe.c_choices")
weighted_choice = module.c_choices.weighted_choice


class TestFeatureMethods(unittest.TestCase):

    def test_weighted_choice_series(self):
        population = list(range(0, 50))
        weights = list(range(0, 50))

        random.seed(5)
        output = weighted_choice(population, weights)
        self.assertEqual(output, 39)

        output = weighted_choice(population, weights)
        self.assertEqual(output, 43)

        output = 0
        for i in range(100):
            output += weighted_choice(population, weights)
        self.assertEqual(output, 3228)

        # Shuffle weights
        random.shuffle(weights)
        output = 0
        for i in range(100):
            output += weighted_choice(population, weights)

        self.assertEqual(output, 2432)
