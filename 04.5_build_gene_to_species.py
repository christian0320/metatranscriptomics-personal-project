#!/usr/bin/env python3

import os
import pandas as pd
from Bio import SeqIO

# Paths
PROJECT = "/home/christianarturo/project1"
CDS_DIR = f"{PROJECT}/data/CDS/individual"
OUTFILE = f"{PROJECT}/data/processed/gene_to_species.csv"

entries = []

print("Scanning .ffn files to extract gene → species mapping...\n")

for f in os.listdir(CDS_DIR):
    if not f.endswith(".ffn"):
        continue
    
    species = f.replace(".ffn", "")  # e.g., Root100.ffn → Root100
    filepath = os.path.join(CDS_DIR, f)

    for record in SeqIO.parse(filepath, "fasta"):
        gene_id = record.id   # the full ORF name
        entries.append([gene_id, species])

df = pd.DataFrame(entries, columns=["gene_id", "species"])

df.to_csv(OUTFILE, index=False)

print("gene_to_species.csv successfully created!")
print(f"Output file: {OUTFILE}")
print(f"Total genes processed: {len(df)}")
