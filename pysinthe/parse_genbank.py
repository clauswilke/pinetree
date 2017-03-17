from Bio import Entrez, SeqIO
from yaml import dump


Entrez.email = "benjamin.r.jack@gmail.com"

# NOTE: Read genbank record notes carefully to interpret genomic coordinates
#
# "Promoters" are really transcription start sites!!
# Likewise, "terminators" are really the final position in a transcript

# Download T7 wild-type and E. coli K12 genbank records
handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")

records = SeqIO.parse(handle, "genbank")

output = {"elements": []}

length_sum = 0

for record in records:

    for feature in record.features:
        # Grab coding sequences
        print(feature)
        start = feature.location.start.position
        stop = feature.location.end.position
        if feature.type == "regulatory":
            length = feature.location.end.position - feature.location.start.position
            if "promoter" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({"type": "promoter",
                                           "name": feature.qualifiers["note"][0],
                                           "start": start,
                                           "stop": stop,
                                           "interactions": {"name": {"binding_constant": " "}}})
            if "terminator" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({"type": "terminator",
                                           "name": feature.qualifiers["note"][0],
                                            "start": start,
                                            "stop": stop,
                                           "interactions": {"name": {"efficiency": " "}}})
        if feature.type == "gene":
            # Record some informationa about the sequence for the FASTA header
            gene_name = feature.qualifiers["note"][0]


            # Construct a string in FASTA format
            output["elements"].append({"type": "transcript",
                                       "name": gene_name,
                                        "start": start,
                                        "stop": stop,
                                       "rbs": -15})
print(dump(output, default_flow_style=False))
