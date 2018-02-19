import pinetree as pt


def execute(output):

    sim = pt.Simulation(cell_volume=8e-16)
    sim.seed(34)
    sim.add_polymerase(name="rnapol", copy_number=20, speed=40, footprint=10)
    sim.add_polymerase(name="ribosome", copy_number=1,
                       speed=30, footprint=10)

    plasmid = pt.Genome(name="T7", length=110,
                        transcript_degradation_rate=1e8)

    plasmid.add_promoter(name="phi1", start=1, stop=10,
                         interactions={"rnapol": 2e8})
    plasmid.add_terminator(name="t1", start=109, stop=110,
                           efficiency={"rnapol": 1.0})
    plasmid.add_gene(name="proteinX", start=50, stop=100,
                     rbs_start=(50 - 15), rbs_stop=50, rbs_strength=1e7)

    sim.register_genome(plasmid)

    sim.run(stop_time=500, time_step=25, output_prefix=output)


if __name__ == "__main__":
    execute("")
