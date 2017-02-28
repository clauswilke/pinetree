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

if __name__ == '__main__':
    unittest.main()
