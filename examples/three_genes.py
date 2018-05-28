import pinetree as pt


def execute(output):

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=10, speed=40, footprint=10)
    sim.add_ribosome(copy_number=100, speed=30, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=605)

    plasmid.add_promoter(name="p1", start=1, stop=10,
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

    sim.simulate(time_limit=60, time_step=1, output=output + "_counts.tsv")


if __name__ == "__main__":
    execute("three_genes")
