import unittest
import pinetree as pt


class TestPromoterMethods(unittest.TestCase):
    def test_init(self):
        with self.assertRaises(ValueError):
            prom = pt.Promoter("mypromoter", -1, 10, dict())
        with self.assertRaises(ValueError):
            prom = pt.Promoter("mypromoter", -1, -10, dict())
        with self.assertRaises(ValueError):
            prom = pt.Promoter("mypromoter", 1, -10, dict())


if __name__ == '__main__':
    unittest.main()
