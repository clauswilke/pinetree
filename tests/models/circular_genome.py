import pinetree as pt


def execute(output):

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase_with_readthrough(name="rnapol", copy_number=5, speed=1, footprint=10)
    sim.add_ribosome(copy_number=5, speed=1, footprint=10)

    plasmid = pt.Genome(name="plasmid", length=450)

    plasmid.add_promoter(name="p1", start=1, stop=10,
                        interactions={"rnapol": 2e10})
    plasmid.add_gene(name="proteinX", start=26, stop=448,
                    rbs_start=(26 - 15), rbs_stop=26, rbs_strength=1e7)

    sim.register_genome(plasmid)

    sim.simulate(time_limit=500, time_step=1, output=output + "_counts.tsv")


if __name__ == "__main__":
    execute("circular_genome")
