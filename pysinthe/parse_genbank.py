from Bio import Entrez, SeqIO
from yaml import dump


Entrez.email = "benjamin.r.jack@gmail.com"

# Download T7 wild-type and E. coli K12 genbank records
handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")

records = SeqIO.parse(handle, "genbank")

output = {"elements": []}

length_sum = 0

for record in records:

    for feature in record.features:
        # Grab coding sequences
        if feature.type == "regulatory":
            print(feature)
            length = feature.location.end.position - feature.location.start.position
            if "promoter" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({"type": "promoter",
                                           "name": feature.qualifiers["note"][0],
                                           "length": length,
                                           "interactions": {"name": {"binding_constant": " "}}})
            if "terminator" in feature.qualifiers["regulatory_class"]:
                output["elements"].append({"type": "terminator",
                                           "name": feature.qualifiers["note"][0],
                                           "length": length,
                                           "interactions": {"name": {"efficiency": " "}}})
            length_sum += length
        if feature.type == "CDS":
            # Record some informationa about the sequence for the FASTA header
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
            elif "product" in feature.qualifiers:
                gene_name = feature.qualifiers["product"][0]
            if gene_name == "hypothetical protein":
                continue

            # Everything should have a locus tag
            if "locus_tag" in feature.qualifiers:
                id = feature.qualifiers["locus_tag"][0]

            # Grab protein ID if applicable
            prot_id = "NA"
            if "protein_id" in feature.qualifiers:
                prot_id = feature.qualifiers["protein_id"][0]

            length = feature.location.end.position - feature.location.start.position

            # Construct a string in FASTA format
            output["elements"].append({"type": "transcript",
                                       "name": gene_name,
                                       "length": length,
                                       "rbs": -15})
            length_sum += length
print(length_sum)
print(dump(output, default_flow_style=False))
