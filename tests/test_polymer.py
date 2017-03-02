import unittest
import random
from pysinthe import feature, polymer


class TestPolymerMethods(unittest.TestCase):

    def setUp(self):
        promoter = feature.Promoter("promoter1",
                                    5,
                                    15,
                                    ["ecolipol", "rnapol"]
                                    )
        terminator = feature.Terminator("myterm",
                                        50,
                                        55,
                                        {"rnapol": {"efficiency": 1.0},
                                         "ecolipol": {"efficiency": 0.6}
                                         })
        mask = feature.Mask("mask",
                            10,
                            100,
                            ["ecolipol"])

        self.polymer = polymer.Polymer("mygenome",
                                       100,
                                       [promoter, terminator],
                                       mask)

        self.assertTrue(terminator.is_covered())
        self.assertTrue(promoter.is_covered())
        self.assertEqual(self.polymer.uncovered["promoter1"], 0)
        self.assertEqual(self.polymer.uncovered["myterm"], 0)
        self.fired = False
        self.polymer.propensity_signal.connect(self.fire)

    def fire(self):
        self.fired = True

    def test_bind_polymerase(self):
        random.seed(22)
        pol = feature.Polymerase("ecolipol",
                                 20,
                                 10,
                                 30,
                                 ["ecolipol", "terminator", "rnapol"]
                                 )
        # Promoter should be covered and inaccessible
        self.assertRaises(RuntimeError,
                          self.polymer.bind_polymerase, pol, "promoter1")
        # Shift mask back 10 positions
        for i in range(10):
            self.polymer.shift_mask()
        # Bind polymerase
        self.polymer.bind_polymerase(pol, "promoter1")
        # Check changes in coverings and positions
        self.assertEqual(pol.start, 5)
        self.assertEqual(pol.stop, 15)
        self.assertEqual(self.polymer.uncovered["promoter1"], 0)
        self.assertEqual(self.polymer.prop_sum, 30)
        self.assertTrue(self.fired)

    def test_insert_polymerase(self):
        pol = feature.Polymerase("ecolipol",
                                 20,
                                 10,
                                 30,
                                 []
                                 )
        pol2 = feature.Polymerase("ecolipol",
                                  60,
                                  10,
                                  30,
                                  []
                                  )
        pol3 = feature.Polymerase("ecolipol",
                                  40,
                                  10,
                                  30,
                                  []
                                  )
        # Clear out any polymerases that may be on polymer
        self.polymer.polymerases = []

        self.polymer._insert_polymerase(pol2)
        self.assertEqual(self.polymer.polymerases.index(pol2), 0)

        self.polymer._insert_polymerase(pol)
        self.assertEqual(self.polymer.polymerases.index(pol), 0)
        self.assertEqual(self.polymer.polymerases.index(pol2), 1)

        self.polymer._insert_polymerase(pol3)
        self.assertEqual(self.polymer.polymerases.index(pol), 0)
        self.assertEqual(self.polymer.polymerases.index(pol2), 2)
        self.assertEqual(self.polymer.polymerases.index(pol3), 1)

        # Trying to insert the same pol object twice should throw an error
        self.assertRaises(RuntimeError, self.polymer._insert_polymerase, pol2)

    def test_count_uncovered(self):
        # Check that cached count matches true count
        count = 0
        for element in self.polymer.elements:
            if element.name == "promoter1" and not element.is_covered():
                count += 1
        self.assertEqual(self.polymer.count_uncovered("promoter1"), count)

        count = 0
        for element in self.polymer.elements:
            if element.name == "myterm" and not element.is_covered():
                count += 1
        self.assertEqual(self.polymer.count_uncovered("myterm"), count)

    def test_calculate_propensity(self):
        prop_sum = 0
        for pol in self.polymer.polymerases:
            prop_sum += pol.speed
        self.assertEqual(prop_sum, self.polymer.calculate_propensity())
