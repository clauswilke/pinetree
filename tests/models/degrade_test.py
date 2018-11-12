import pinetree as pt


def execute(output):

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=10, speed=30, footprint=10)
    sim.add_ribosome(copy_number=100,
                     speed=20, footprint=10)

    plasmid = pt.Genome(name="T7", length=305,
                        transcript_degradation_rate=1e-2,
                        transcript_degradation_rate_ext=1e-2,
                        rnase_footprint=9,
                        rnase_speed=20)

    # plasmid = pt.Genome(name="T7", length=305)

    plasmid.add_promoter(name="phi1", start=1, stop=10,
                         interactions={"rnapol": 2e8})
    plasmid.add_terminator(name="t1", start=304, stop=305,
                           efficiency={"rnapol": 1.0})
    plasmid.add_gene(name="proteinX", start=30, stop=99,
                     rbs_start=(30 - 10), rbs_stop=30, rbs_strength=1e7)
    # plasmid.add_promoter(name="phi2", start=100, stop=109,
    #                      interactions={"rnapol": 2e8})
    # plasmid.add_terminator(name="t1", start=98, stop=99,
    #                        efficiency={"rnapol": 1.0})
    # plasmid.add_rnase_site(100, 110)
    plasmid.add_gene(name="proteinY", start=120, stop=199,
                     rbs_start=(120 - 10), rbs_stop=120, rbs_strength=1e7)
    # plasmid.add_promoter(name="phi3", start=200, stop=209,
    #                     interactions={"rnapol": 2e8})
    plasmid.add_gene(name="proteinZ", start=220, stop=300,
                     rbs_start=(220 - 10), rbs_stop=220, rbs_strength=1e7)

    sim.register_genome(plasmid)

    sim.simulate(time_limit=500, time_step=1, output=output + "counts.tsv")


if __name__ == "__main__":
    execute("three_genes_rnase_")
