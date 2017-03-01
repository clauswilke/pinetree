import unittest
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
                                      20,
                                      40,
                                      30,
                                      ["mypol", "terminator", "mypol2"]
                                      )
        self.pol2 = feature.Polymerase("mypol2",
                                       20,
                                       40,
                                       30,
                                       ["mypol", "terminator"]
                                       )

    def test_init(self):
        self.assertEqual(self.pol.footprint, self.pol.stop - self.pol.start)

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


if __name__ == '__main__':
    unittest.main()
