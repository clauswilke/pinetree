import pinetree as pt


def execute(output):

    sim = pt.Simulation(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=1, speed=30, footprint=10)
    sim.add_polymerase(name="ribosome", copy_number=1,
                       speed=5, footprint=10)

    plasmid = pt.Genome(name="T7", length=205,
                        transcript_degradation_rate=1e9)

    plasmid.add_promoter(name="phi1", start=1, stop=10,
                         interactions={"rnapol": 2e8})
    plasmid.add_terminator(name="t1", start=204, stop=205,
                           efficiency={"rnapol": 1.0})
    plasmid.add_gene(name="proteinX", start=30, stop=100,
                     rbs_start=(30 - 10), rbs_stop=30, rbs_strength=1e7)
    plasmid.add_gene(name="proteinY", start=120, stop=200,
                     rbs_start=(120 - 10), rbs_stop=120, rbs_strength=1e7)

    sim.register_genome(plasmid)

    sim.run(stop_time=200, time_step=1, output_prefix=output)


if __name__ == "__main__":
    execute("")
