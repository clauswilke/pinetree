import unittest
import random
from pysinthe import feature


class TestFeatureMethods(unittest.TestCase):

    def test_interactions(self):
        my_feature = feature.Feature("test",
                                     0,
                                     10,
                                     ["test", "test2"])
        self.assertTrue(my_feature.check_interaction("test2"))
        self.assertFalse(my_feature.check_interaction("test3"))


class TestPolymeraseMethods(unittest.TestCase):

    def setUp(self):
        self.pol = feature.Polymerase("mypol",
                                      40,
                                      30
                                      )
        self.pol2 = feature.Polymerase("mypol2",
                                       40,
                                       30
                                       )

    def test_init(self):
        self.assertEqual(self.pol.footprint, self.pol.stop - self.pol.start)
        self.assertEqual(self.pol.type, "polymerase")

    def test_move(self):
        old_start = self.pol.start
        old_stop = self.pol.stop
        self.pol.move()
        self.assertEqual(self.pol.start, old_start + 1)
        self.assertEqual(self.pol.stop, old_stop + 1)

    def test_move_back(self):
        old_start = self.pol.start
        old_stop = self.pol.stop
        self.pol.move_back()
        self.assertEqual(self.pol.start, old_start - 1)
        self.assertEqual(self.pol.stop, old_stop - 1)


class TestMaskMethods(unittest.TestCase):
    def setUp(self):
        self.mask = feature.Mask("mask",
                                 10,
                                 100,
                                 ["pol"]
                                 )

    def test_recede(self):
        old_start = self.mask.start
        self.mask.recede()
        self.assertEqual(self.mask.start, old_start + 1)


class TestElementMethods(unittest.TestCase):
    def setUp(self):
        self.element = feature.Element("myelement",
                                       23,
                                       60,
                                       ["ribosome"]
                                       )

    def test_save_state(self):
        self.element._covered = 3
        self.element.save_state()
        self.assertEqual(self.element._old_covered, 3)
        self.assertEqual(self.element._covered, self.element._old_covered)

    def test_cover(self):
        old_covered = self.element._covered
        self.element.cover()
        self.assertEqual(self.element._covered, old_covered + 1)

    def test_uncover(self):
        self.element._covered = 2
        self.element.uncover()
        self.assertEqual(self.element._covered, 1)
        # Make sure covering count doesn't go below zero
        self.element._covered = 0
        self.element.uncover()
        self.assertTrue(self.element._covered >= 0)

    def test_was_covered(self):
        self.element._covered = 0
        self.element.save_state()
        self.element.cover()
        self.assertTrue(self.element.was_covered())
        self.element.save_state()
        self.element.cover()
        self.assertFalse(self.element.was_covered())

    def test_was_uncovered(self):
        self.element._covered = 3
        self.element.save_state()
        self.element.uncover()
        self.element.uncover()
        self.element.save_state()
        self.assertFalse(self.element.was_uncovered())
        self.element.save_state()
        self.element.uncover()
        self.assertTrue(self.element.was_uncovered())

    def test_is_covered(self):
        self.element._covered = 1
        self.assertTrue(self.element.is_covered())
        self.element._covered = 0
        self.assertFalse(self.element.is_covered())


class TestTerminatorMethods(unittest.TestCase):
    def setUp(self):
        random.seed(22)
        self.term = feature.Terminator("myterm",
                                       23,
                                       60,
                                       {"rnapol": {"efficiency": 1.0},
                                        "ecolipol": {"efficiency": 0.6}
                                        }
                                       )
        self.pol = feature.Polymerase("rnapol",
                                      40,
                                      30
                                      )
        self.pol2 = feature.Polymerase("ecolipol",
                                       40,
                                       30
                                       )
        self.term.gene = "mygene"

        # Create temp termination signal
        self.pol.termination_signal.connect(lambda x: self.assertEqual(x, 60))

    def test_init(self):
        self.assertFalse(self.term.readthrough)


if __name__ == '__main__':
    unittest.main()
