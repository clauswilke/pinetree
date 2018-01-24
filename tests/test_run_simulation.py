# Test simulation
import unittest
import subprocess
import tempfile
import importlib
import pinetree.pinetree as pt


class MainTest(unittest.TestCase):

    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tempdir.cleanup()

    def run_test(self, prefix):
        # test_mod = importlib.import_module('params.' + prefix)
        out_prefix = self.tempdir.name + "/" + prefix
        # test_mod.execute(out_prefix)
        with open('tests/output/' + prefix + '_counts.tsv') as f:
            text = f.read()

        sim = pt.Simulation(cell_volume=8e-16)
        sim.seed(34)

        sim.add_polymerase(name="rnapol", copy_number=1,
                           speed=40, footprint=10)
        sim.add_polymerase(name="ribosome", copy_number=1,
                           speed=30, footprint=10)

        plasmid = pt.Genome(name="T7", length=605)

        plasmid.add_promoter(name="phi1", start=1, stop=10,
                             interactions={"rnapol": 2e8})
        plasmid.add_terminator(name="t1", start=604, stop=605,
                               efficiency={"rnapol": 1.0})

        plasmid.add_gene(name="rnapol", start=26, stop=225,
                         rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)
        plasmid.add_gene(name="proteinX", start=241, stop=280,
                         rbs_start=(241 - 15), rbs_stop=241, rbs_strength=1e7)
        plasmid.add_gene(name="proteinY", start=296, stop=595,
                         rbs_start=(296 - 15), rbs_stop=296, rbs_strength=1e7)

        sim.register_genome(plasmid)
        sim.run(stop_time=40, time_step=1, output_name=out_prefix)

        with open(out_prefix + '_counts.tsv') as results_file:
            results = results_file.read()
        self.assertEqual(results, text)

    def test_single_gene(self):
        self.run_test('single_gene')

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
