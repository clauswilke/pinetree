from Bio import Entrez, SeqIO
from yaml import dump


def get_promoter_interactions(name):
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phage_weak = ["T7 promoter phi1.1A",
                  "T7 promoter phi1.1B",
                  "T7 promoter phi1.5",
                  "T7 promoter phi1.6"]
    phage_strong = ["T7 promoter phi2.5",
                    "T7 promoter phi3.8",
                    "T7 promoter phi4c",
                    "T7 promoter phi4.3",
                    "T7 promoter phi4.7",
                    "T7 promoter phi6.5",
                    "T7 promoter phi9",
                    "T7 promoter phi10",
                    "T7 promoter phi13",
                    "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': {'binding_constant': 10e7},
                'ecolipol-p': {'binding_constant': 3e7}}
    elif name in ecoli_weak:
        return {'ecolipol': {'binding_constant': 1e7},
                'ecolipol-p': {'binding_constant': 0.3e7}}
    elif name in phage_weak:
        return {'rnapol-1': {'binding_constant': 3.3e7},
                'rnapol-3.5': {'binding_constant': 1.6e7}}
    elif name in phage_strong:
        return {'rnapol-1': {'binding_constant': 1.8e8},
                'rnapol-3.5': {'binding_constant': 0.9e8}}
    else:
        return {'name': {'binding_constant': 0.0}}

def get_terminator_interactions(name):

    if name == "E. coli transcription terminator TE":
        return {'ecolipol': {'efficiency': 1.0},
                'ecolipol-p': {'efficiency': 1.0},
                'rnapol-1': {'efficiency': 0.0},
                'rnapol-3.5': {'efficiency': 0.0}}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': {'efficiency': 0.8},
                'rnapol-3.5': {'efficiency': 0.8}}
    else:
        return {'name': {'efficiency': 0.0}}

ignore_genes = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 4.1",
                "gene 4B",
                "gene 1.1",
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

output = {"elements": []}
old_start = 0
old_stop = 0
start = 0
stop = 0

for record in records:

    for feature in record.features:
        print(feature)
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
            # Grab the gene name
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

print(dump(output, default_flow_style=False))
