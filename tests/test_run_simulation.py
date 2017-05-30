# Test simulation
import unittest
import subprocess
import tempfile

class MainTest(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
    
    def tearDown(self):
        self.tempdir.cleanup()
    
    def run_test(self, prefix):
        out_prefix = self.tempdir.name + "/" + prefix
        subprocess.call(['bin/pysinthe-run', 
                         'tests/params/' + prefix + '.yml', 
                         '-o', out_prefix])

        with open('tests/output/' + prefix + '_counts.tsv') as f:
            text = f.read()
        with open(out_prefix + '_counts.tsv') as results_file:
            results = results_file.read()

        self.assertEqual(results, text)

    def test_three_genes(self):
        self.run_test('three_genes')

    def test_dual_polymerases(self):
        self.run_test('dual_polymerases')

    def test_dual_promoters(self):
        self.run_test('dual_promoter')

    def test_readthrough(self):
        self.run_test('readthrough')

    def test_genome_entry(self):
        self.run_test('genome_entry')

    def test_consecutive_promoters(self):
        self.run_test('consecutive_promoters')

    def test_lotka_voltera(self):
        self.run_test('lotka_voltera')

    # def test_promoter_gene_overlap(self):
    #     stdout = sys.stdout  # keep a handle on the real standard output
    #     results = io.StringIO()
    #     sys.stdout = results  # Choose a file-like object to write to
    #     simulation.main("tests/params/promoter_gene_overlap.yml")
    #     sys.stdout = stdout

    #     with open('tests/output/promoter_gene_overlap_out.csv') as f:
    #         text = f.read()

    #     self.assertEqual(results.getvalue(), text)

    # def test_three_genes_runoff(self):
    #     stdout = sys.stdout  # keep a handle on the real standard output
    #     results = io.StringIO()
    #     sys.stdout = results  # Choose a file-like object to write to
    #     simulation.main("tests/params/three_genes_runoff.yml")
    #     sys.stdout = stdout

    #     with open('tests/output/three_genes_runoff_out.csv') as f:
    #         text = f.read()

    #     self.assertEqual(results.getvalue(), text)

    # def test_overlapping_genes(self):
    #     stdout = sys.stdout  # keep a handle on the real standard output
    #     results = io.StringIO()
    #     sys.stdout = results  # Choose a file-like object to write to
    #     simulation.main("tests/params/overlapping_genes.yml")
    #     sys.stdout = stdout

    #     with open('tests/output/overlapping_genes_out.csv') as f:
    #         text = f.read()

    #     self.assertEqual(results.getvalue(), text)


if __name__ == '__main__':
    unittest.main()
