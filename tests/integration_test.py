# Test simulation
import unittest
import subprocess
import tempfile
import importlib


class MainTest(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def run_test(self, prefix):
        test_mod = importlib.import_module('.models.' + prefix, 'tests')
        out_prefix = self.tempdir.name + "/" + prefix
        test_mod.execute(out_prefix)
        with open('tests/output/' + prefix + '_counts.tsv') as f:
            text = f.read()
        with open(out_prefix + '_counts.tsv') as results_file:
            results = results_file.read()
        self.assertEqual(results, text)

    def test_single_gene(self):
        self.run_test('single_gene')
    
    def test_circular_genome(self):
        self.run_test("circular_genome")

    # def test_three_genes(self):
    #     self.run_test('three_genes')

    # def test_dual_polymerases(self):
    #     self.run_test('dual_polymerases')

    # def test_dual_promoters(self):
    #     self.run_test('dual_promoter')

    # def test_readthrough(self):
    #     self.run_test('readthrough')

    # def test_genome_entry(self):
    #     self.run_test('genome_entry')

    # def test_consecutive_promoters(self):
    #     self.run_test('consecutive_promoters')

    # def test_lotka_voltera(self):
    #     self.run_test('lotka_voltera')

    # def test_promoter_gene_overlap(self):
    #     self.run_test('promoter_gene_overlap')

    # def test_three_genes_runoff(self):
    #     self.run_test('three_genes_runoff')

    # def test_overlapping_genes(self):
    #     self.run_test('overlapping_genes')


if __name__ == '__main__':
    unittest.main()
