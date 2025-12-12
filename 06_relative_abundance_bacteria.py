#!/usr/bin/env python3

import os
import pandas as pd

# ================================
# PATHS
# ================================
PROJECT = "/home/christianarturo/project1"
TPM_FILE = f"{PROJECT}/data/processed/tpm_matrix_bacteria.csv"
MAP_FILE = f"{PROJECT}/data/processed/gene_to_species.csv"
OUTFILE = f"{PROJECT}/data/processed/relative_abundance_species.csv"

# ================================
# LOAD TPM AND MAPPING
# ================================
print("Loading TPM matrix...")
tpm = pd.read_csv(TPM_FILE)

print("Loading gene-to-species mapping...")
mapping = pd.read_csv(MAP_FILE)

# ================================
# MERGE TPM WITH SPECIES MAP
# ================================
print("Merging TPM matrix with species mapping...")
df = tpm.merge(mapping, left_on="Name", right_on="gene_id", how="left")

df = df.drop(columns=["gene_id"])

# ================================
# AGGREGATE TPM BY SPECIES
# ================================
print("Aggregating TPM by species...")

species_abundance = df.groupby("species").sum(numeric_only=True)

# ================================
# SAVE OUTPUT
# ================================
species_abundance.to_csv(OUTFILE)

print("-------------------------------------------------------")
print("Relative species abundance successfully generated.")
print(f"Output file: {OUTFILE}")
print(f"Dimensions: {species_abundance.shape}")
print("-------------------------------------------------------")
