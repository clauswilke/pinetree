from Bio import Entrez, SeqIO
from yaml import dump, safe_load

opt_codons_E_coli = { 'A':['GCT'], 'R':['CGT', 'CGC'], 'N':['AAC'], 'D':['GAC'], 'C':['TGC'], 'Q':['CAG'], 'E':['GAA'], 'G':['GGT','GGC'], 'H':['CAC'], 'I':['ATC'], 'L':['CTG'], 'F':['TTC'], 'P':['CCG'], 'S':['TCT','TCC'], 'T':['ACT','ACC'], 'Y':['TAC'], 'V':['GTT','GTA'] }

def set_up():
    data = """
    simulation:
        seed: 34
        runtime: 1800
        time_step: 5
        cell_volume: 1.1e-15
        debug: False
    # Genome parameters
    genome:
        name: T7
        copy_number: 1
        entered: 500
        mask_interactions:
            - rnapol-1
            - rnapol-3.5
            - ecolipol
            - ecolipol-p
            - ecolipol-2
            - ecolipol-2-p

    polymerases:
    - name: rnapol-1
      copy_number: 0
      speed: 230
      footprint: 35
    - name: rnapol-3.5
      copy_number: 0
      speed: 230
      footprint: 35
    - name: ecolipol
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-p
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-2
      copy_number: 0
      speed: 45
      footprint: 35
    - name: ecolipol-2-p
      copy_number: 0
      speed: 45
      footprint: 35

    ribosomes:
    - name: ribosome
      copy_number: 0
      speed: 30
      footprint: 30
      binding_constant: 1e7

    species:
    - name: bound_ribosome
      copy_number: 10000 # Assume 10000 total ribosomes bound
    - name: bound_ecolipol
      copy_number: 1800
    - name: bound_ecolipol_p
      copy_number: 0
    - name: ecoli_genome
      copy_number: 0
    - name: ecoli_transcript
      copy_number: 0

    reactions:
    - name: ecoli_transcripts # Assumes that each ribosome binds to one RBS
      propensity: 1e6
      reactants:
          - ecoli_transcript
          - ribosome
      products:
          - bound_ribosome
    - name: ecoli_transcripts_reverse # Assumes that ~1000nt transcript and 30bp/s
      propensity: 0.04
      reactants:
          - bound_ribosome
      products:
          - ribosome
          - ecoli_transcript
    - name: ecoli_transcript_degradation
      propensity: 0.001925  # 1st-order decay based on 6min half-life
      reactants:
          - ecoli_transcript
      products:
          - degraded_transcript
    - name: ecoli_promoters
      propensity: 1e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol
          - ecoli_genome
      products:
          - bound_ecolipol
    - name: ecoli_promoters_p
      propensity: 0.3e7 # Corresponds to weak E. coli promoter
      reactants:
          - ecolipol-p
          - ecoli_genome
      products:
          - bound_ecolipol_p
    - name: ecoli_promoters_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol
      products:
          - ecolipol
          - ecoli_genome
          - ecoli_transcript
    - name: ecoli_promoters_p_reverse # Assumes 1 new RBS per 1000nt gene
      propensity: 0.04
      reactants:
          - bound_ecolipol_p
      products:
          - ecolipol-p
          - ecoli_genome
          - ecoli_transcript
    - name: reaction1
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol
      products:
          - ecolipol-p
          - protein_kinase-0.7
    - name: reaction2
      propensity: 3.8e7
      reactants:
          - protein_kinase-0.7
          - ecolipol-2
      products:
          - ecolipol-2-p
          - protein_kinase-0.7
    - name: reaction3
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol
      products:
          - ecolipol-2
    - name: reaction4
      propensity: 3.8e7
      reactants:
          - gp-2
          - ecolipol-p
      products:
          - ecolipol-2-p
    - name: reaction5
      propensity: 1.1
      reactants:
          - ecolipol-2
      products:
          - gp-2
          - ecolipol
    - name: reaction6
      propensity: 1.1
      reactants:
          - ecolipol-2-p
      products:
          - gp-2
          - ecolipol-p
    - name: reaction7
      propensity: 3.8e7
      reactants:
          - lysozyme-3.5
          - rnapol-1
      products:
          - rnapol-3.5
    - name: reaction8
      propensity: 3.5
      reactants:
          - rnapol-3.5
      products:
          - rnapol-1
          - lysozyme-3.5
    """

    return safe_load(data)

def get_promoter_interactions(name):
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    phi10_bind = 1.82e7 # Binding constant for phi10

    if name in ecoli_strong:
        return {'ecolipol': {'binding_constant': 10e4},
                'ecolipol-p': {'binding_constant': 3e4}}
    elif name in ecoli_weak:
        return {'ecolipol': {'binding_constant': 1e4},
                'ecolipol-p': {'binding_constant': 0.3e4}}
    elif name in phi1_3:
        return {'rnapol-1': {'binding_constant': phi10_bind*0.01},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.01*0.5}}
    elif name in phi3_8:
        return {'rnapol-1': {'binding_constant': phi10_bind*0.01},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.01*0.5}}
    elif name in phi6_5:
        return {'rnapol-1': {'binding_constant': phi10_bind*0.05},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.05*0.5}}
    elif name in phi9:
        return {'rnapol-1': {'binding_constant': phi10_bind*0.2},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.2*0.5}}
    elif name in phi10:
        return {'rnapol-1': {'binding_constant': phi10_bind},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.5}}
    elif name in phi13:
        return {'rnapol-1': {'binding_constant': phi10_bind*0.1},
                'rnapol-3.5': {'binding_constant': phi10_bind*0.1*0.5}}
    else:
        raise ValueError("Promoter strength for {0} not assigned.".format(name))

def get_terminator_interactions(name):

    if name == "E. coli transcription terminator TE":
        return {'ecolipol': {'efficiency': 1.0},
                'ecolipol-p': {'efficiency': 1.0},
                'rnapol-1': {'efficiency': 0.0},
                'rnapol-3.5': {'efficiency': 0.0}}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': {'efficiency': 0.85},
                'rnapol-3.5': {'efficiency': 0.85}}
    else:
        return {'name': {'efficiency': 0.0}}

ignore_genes = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 4.1",
                "gene 4B",
                "gene 0.6A",
                "gene 0.6B",
                "possible gene 0.6B",
                "gene 0.5",
                "gene 0.4"]
ignore_regulatory = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

Entrez.email = "benjamin.r.jack@gmail.com"

# NOTE: Read genbank record notes carefully to interpret genomic coordinates
#
# "Promoters" are really transcription start sites!!
# Likewise, "terminators" are really the final position in a transcript

# Download T7 wild-type genbank records
handle = Entrez.efetch(db="nuccore",
                       id=["NC_001604"],
                       rettype="gb",
                       retmode="text")

records = SeqIO.parse(handle, "genbank")

output = set_up()

output["elements"] =[]
old_start = 0
old_stop = 0
start = 0
stop = 0



for record in records:

    translation_scale_factors = [1.0]*len(record.seq)

    for feature in record.features:
        # print(feature)
        start = feature.location.start.position
        stop = feature.location.end.position
        # Grab promoters and terminators
        if feature.type == "regulatory":
            name = feature.qualifiers["note"][0]
            if name in ignore_regulatory:
                continue
            # Construct promoter params
            if "promoter" in feature.qualifiers["regulatory_class"]:
                if stop - start < 35:
                    promoter_start = start - 35
                interactions = get_promoter_interactions(name)
                output["elements"].append({
                     "type": "promoter",
                     "name": name,
                     "start": promoter_start,
                     "stop": stop,
                     "interactions": interactions
                })
            # Construct terminator params
            if "terminator" in feature.qualifiers["regulatory_class"]:
                interactions = get_terminator_interactions(name)
                output["elements"].append({
                    "type": "terminator",
                    "name": name,
                    "start": start,
                    "stop": stop,
                    "interactions": interactions
                })
        # Grab genes/CDSes
        if feature.type == "gene":
            temp_name = feature.qualifiers["note"][0]
            if temp_name in ignore_genes:
                continue
            if temp_name == "gene 2":
                gene_name = "gp-2"
            elif temp_name == "gene 1":
                gene_name = "rnapol-1"
            elif temp_name == "gene 3.5":
                gene_name = "lysozyme-3.5"
            elif temp_name == "gene 0.7":
                gene_name = "protein_kinase-0.7"
            else:
                gene_name = temp_name
            # Construct CDS parameters for this gene
            transcript = {"type": "transcript",
                          "name": gene_name,
                          "start": start,
                          "stop": stop,
                          "rbs": -30}
            output["elements"].append(transcript)
        if feature.type == "CDS":
            # Grab the gene name
            nuc_seq = feature.location.extract(record).seq
            aa_seq = feature.qualifiers["translation"][0]
            for index, nuc in enumerate(nuc_seq):
                aa_index = int(index / 3)
                codon_start = aa_index * 3
                codon = nuc_seq[codon_start:codon_start+3]
                genome_index = feature.location.start + index
                if aa_index < len(aa_seq):
                    if aa_seq[aa_index] in opt_codons_E_coli:
                        if codon in opt_codons_E_coli[aa_seq[aa_index]]:
                            translation_scale_factors[genome_index] = 2.1
                        else:
                            translation_scale_factors[genome_index] = 0.15
                # if feature.qualifiers["protein_id"][0] == 'NP_041998.1':
                #     translation_scale_factors[genome_index] = 0.15

# output["genome"]["translation_weights"] = translation_scale_factors
output["genome"]["length"] = len(record.seq)

print(dump(output, default_flow_style=False))
