# Test simulation
import unittest
import sys
import io

import run_simulation as simulation

class MainTest(unittest.TestCase):

    def test_three_genes(self):
        stdout = sys.stdout  # keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results  # Choose a file-like object to write to
        simulation.main("tests/params/three_genes.yml")
        sys.stdout = stdout

        with open('tests/output/three_genes_out.csv') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)

    def test_dual_polymerases(self):
        stdout = sys.stdout  # keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results  # Choose a file-like object to write to
        simulation.main("tests/params/dual_polymerases.yml")
        sys.stdout = stdout

        with open('tests/output/dual_polymerases_out.csv') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)

    def test_dual_promoters(self):
        stdout = sys.stdout  # keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results  # Choose a file-like object to write to
        simulation.main("tests/params/dual_promoter.yml")
        sys.stdout = stdout

        with open('tests/output/dual_promoter_out.csv') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)

    def test_readthrough(self):
        stdout = sys.stdout  # keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results  # Choose a file-like object to write to
        simulation.main("tests/params/readthrough.yml")
        sys.stdout = stdout

        with open('tests/output/readthrough_out.csv') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)

    def test_genome_entry(self):
        stdout = sys.stdout  # keep a handle on the real standard output
        results = io.StringIO()
        sys.stdout = results  # Choose a file-like object to write to
        simulation.main("tests/params/genome_entry.yml")
        sys.stdout = stdout

        with open('tests/output/genome_entry_out.csv') as f:
            text = f.read()

        self.assertEqual(results.getvalue(), text)


if __name__ == '__main__':
    unittest.main()
