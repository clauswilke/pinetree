import pinetree as pt

def execute(output):
    ## parameters
    RB_COPY = 10 # number of ribosomes in simulation
    TS_COPY = 100 # number of transcripts (mRNAs) in simulation
    RBS_STRENGTH = 10000000.0 # strength of ribosome binding to an mRNA
    TRNA_CHRG_RATES = [50.0, 50.0] # strength of tRNA re-charging reaction [tRNA_a, tRNA_b]
    TRNA_PROPORTIONS = (0.1, 0.9) # tRNA proportions, i.e. 90% total tRNA is type A, other 10% is type B
    TOTAL_TRNA = 2500 # total tRNA
    SPEED = 0.5 # rate constant describing how efficient ribosome movement is
    
    # 300 nt coding region + 50 nt buffer (30 on the left, 20 on the right)
    # transcripts coding region is 10% codon "AAA" and 90% codon "TAT"
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATTATTATTATAAATATTATTATTATTATTATTATAAATATTATTATTATTATTATTATTATTATAAATATAAATATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATAAATATTATTATTATTATAAATATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATAAATATTATTATTATTATTATTATTATTATAAATATTATTATAAAAAATATTATTATTATTATTATTATTATAAAAAAAAAAAAAAAAAAAA"

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(1)
    sim.add_ribosome(copy_number=RB_COPY, speed=SPEED, footprint=15)
    
    i = 0
    while i < TS_COPY: # add 100 copies of transcript to the simulation (each copy is identical)
        transcript = pt.Transcript("transcript", 350)
        transcript.add_gene(name="proteinX", start=31, stop=330,
                     rbs_start=(31 - 15), rbs_stop=31, rbs_strength=RBS_STRENGTH)
        transcript.add_seq(seq=seq)
        sim.register_transcript(transcript)
        i += 1
    
    ## tRNA/codon mapping: 
    #
    # tRNA "TTT" corresponds to codon "AAA", tRNA "ATA" corresponds to codon "TAT"
    tRNA_map = {"AAA": ["TTT"], "TAT": ["ATA"]} 
    # initial tRNA species counts:
    #
    # TTT: 250 (total tRNA * 0.1) charged, 0 uncharged
    # ATA: 2250 (total tRNA * 0.9) charged, 0 uncharged
    tRNA_counts = {"TTT": [int(TOTAL_TRNA*TRNA_PROPORTIONS[0]), 0], "ATA": [int(TOTAL_TRNA*TRNA_PROPORTIONS[1]), 0]}
    # tRNA charging rates: TTT: 100.0, ATA: 100.0
    tRNA_rates = {"TTT": TRNA_CHRG_RATES[0], "ATA": TRNA_CHRG_RATES[1]}
    sim.add_trna(tRNA_map, tRNA_counts, tRNA_rates)
    
    sim.simulate(time_limit=15, time_step=5, output=output + "_counts.tsv")

    
if __name__ == "__main__":
    execute("dynamic_trnas")