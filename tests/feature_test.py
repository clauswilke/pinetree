import unittest
import pinetree as pt


class TestBindingSiteMethods(unittest.TestCase):

    def test_init(self):
        with self.assertRaises(ValueError):
            prom = pt.BindingSite("promoter", -1, 10, dict())
        with self.assertRaises(ValueError):
            site = pt.BindingSite("promoter", -1, -10, dict())
        with self.assertRaises(ValueError):
            site = pt.BindingSite("promoter", 1, -10, dict())

    def test_coverings(self):
        site = pt.BindingSite("promoter", 1, 10, {"rnapol": 1.0})
        self.assertFalse(site.was_covered())
        self.assertFalse(site.is_covered())
        self.assertFalse(site.was_uncovered())

        site.cover()
        self.assertTrue(site.is_covered())
        self.assertTrue(site.was_covered())
        self.assertFalse(site.was_uncovered())

        site.reset_state()
        self.assertTrue(site.is_covered())
        self.assertFalse(site.was_covered())
        self.assertFalse(site.was_uncovered())

        site.uncover()
        self.assertFalse(site.is_covered())
        self.assertFalse(site.was_covered())
        self.assertTrue(site.was_uncovered())

        site.reset_state()
        self.assertFalse(site.is_covered())
        self.assertFalse(site.was_covered())
        self.assertFalse(site.was_uncovered())

    def test_interaction(self):
        site = pt.BindingSite("promoter", 1, 10, {"rnapol": 1.0})
        self.assertTrue(site.check_interaction("rnapol"))
        self.assertFalse(site.check_interaction("otherpol"))


if __name__ == '__main__':
    unittest.main()
