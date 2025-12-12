#!/usr/bin/env python3

import os
import pandas as pd

# ==========================================================
# CONFIGURATION
# ==========================================================

PROJECT = "/home/christianarturo/project1"
SALMON_DIR = f"{PROJECT}/data/salmon_quant_bacteria"
OUTDIR = f"{PROJECT}/data/processed"
os.makedirs(OUTDIR, exist_ok=True)

print("Searching for Salmon quantification outputs...")

# ==========================================================
# 1. Locate all quant.sf files
# ==========================================================

samples = []
quant_files = []

for sample in os.listdir(SALMON_DIR):
    sample_dir = os.path.join(SALMON_DIR, sample)
    quant_file = os.path.join(sample_dir, "quant.sf")

    if os.path.isfile(quant_file):
        samples.append(sample)
        quant_files.append(quant_file)

samples = sorted(samples)

print(f"Detected {len(samples)} samples:")
for s in samples:
    print(" -", s)

# ==========================================================
# 2. Load all TPM columns into a single matrix
# ==========================================================

tpm_matrix = pd.DataFrame()

for sample, qf in zip(samples, quant_files):
    df = pd.read_csv(qf, sep="\t")
    df = df[["Name", "TPM"]].rename(columns={"TPM": sample})

    if tpm_matrix.empty:
        tpm_matrix = df.copy()
    else:
        tpm_matrix = tpm_matrix.merge(df, on="Name", how="outer")

tpm_matrix = tpm_matrix.fillna(0)

# ==========================================================
# 3. Save output matrix
# ==========================================================

output_file = f"{OUTDIR}/tpm_matrix_bacteria.csv"
tpm_matrix.to_csv(output_file, index=False)

print("-------------------------------------------------------")
print("TPM matrix successfully generated.")
print(f"Output file: {output_file}")
print("Dimensions:", tpm_matrix.shape)
print("-------------------------------------------------------")
