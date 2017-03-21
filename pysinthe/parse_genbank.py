from Bio import Entrez, SeqIO
from yaml import dump


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

for record in records:

    for feature in record.features:
        # print(feature)
        start = feature.location.start.position
        stop = feature.location.end.position
        # Grab promoters and terminators
        if feature.type == "regulatory":
            # Construct promoter params
            if "promoter" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({
                     "type": "promoter",
                     "name": feature.qualifiers["note"][0],
                     "start": start,
                     "stop": stop,
                     "interactions": {"name": {"binding_constant": " "}}
                })
            # Construct terminator params
            if "terminator" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({
                    "type": "terminator",
                    "name": feature.qualifiers["note"][0],
                    "start": start,
                    "stop": stop,
                    "interactions": {"name": {"efficiency": " "}}
                })
        # Grab genes/CDSes
        if feature.type == "gene":
            # Grab the gene name
            gene_name = feature.qualifiers["note"][0]
            # Construct CDS parameters for this gene
            output["elements"].append({"type": "transcript",
                                       "name": gene_name,
                                       "start": start,
                                       "stop": stop,
                                       "rbs": -15})

print(dump(output, default_flow_style=False))
