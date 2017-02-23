# Test simulation
import unittest
import sys
import io

import simulation

class MainTest(unittest.TestCase):

    def test_three_genes(self):
        stdout = sys.stdout  #keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results #Choose a file-like object to write to
        simulation.main("examples/three_genes.yml")
        sys.stdout = stdout

        with open('tests/integration/three_genes.out') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)

if __name__ == '__main__':
    unittest.main()
