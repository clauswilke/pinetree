import unittest
from pysinthe import feature


class TestFeatureMethods(unittest.TestCase):

    def test_interactions(self):
        my_feature = feature.Feature("test",
                                     0,
                                     10,
                                     ["test", "test2"])
        my_feature2 = feature.Feature("test2",
                                     0,
                                     10,
                                     ["test2"])
        self.assertTrue(my_feature.check_interaction(my_feature2))
        self.assertFalse(my_feature2.check_interaction(my_feature))

if __name__ == '__main__':
    unittest.main()
