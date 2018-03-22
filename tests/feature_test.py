import unittest
import pinetree as pt


class TestFixedElementMethods(unittest.TestCase):
    def test_init(self):
        # FixedElement is an abstract class, we shouldn't be able to instantiate
        with self.assertRaises(TypeError):
            element = pt.FixedElement()


if __name__ == '__main__':
    unittest.main()
