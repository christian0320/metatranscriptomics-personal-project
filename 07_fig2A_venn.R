#!/usr/bin/env Rscript

library(tidyverse)
library(VennDiagram)

# ======================================================
# PATHS
# ======================================================

PROJECT <- "/home/christianarturo/project1"
TPM_FILE <- file.path(PROJECT, "data/processed/tpm_matrix_bacteria.csv")
OUTDIR <- file.path(PROJECT, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ======================================================
# 1. Load TPM matrix
# ======================================================

tpm <- read_csv(TPM_FILE)

gene_ids <- tpm$Name
tpm_data <- tpm %>% select(-Name)

# ======================================================
# 2. Apply expression threshold (TPM > 1)
# ======================================================

THRESHOLD <- 1
binary_matrix <- tpm_data > THRESHOLD

# Remove genes not expressed in ANY sample
expressed_anywhere <- rowSums(binary_matrix) > 0

gene_ids <- gene_ids[expressed_anywhere]
binary_matrix <- binary_matrix[expressed_anywhere, ]

# ======================================================
# 3. Groups used in the paper (Matrix vs Root)
# ======================================================

matrix_samples <- c("SRR24611159", "SRR24611160", "SRR24611161")
root_samples   <- c("SRR24611162", "SRR24611163", "SRR24611164")

genes_matrix <- gene_ids[rowSums(binary_matrix[, matrix_samples]) > 0]
genes_root   <- gene_ids[rowSums(binary_matrix[, root_samples]) > 0]

# ======================================================
# 4. Plot 2-circle Venn Diagram
# ======================================================

venn.diagram(
  x = list(
    Matrix = genes_matrix,
    Root   = genes_root
  ),
  filename = file.path(OUTDIR, "Fig2A_Venn.png"),
  output = TRUE,
  height = 1800,
  width = 1800,
  resolution = 300,
  fill = c("#7EC8E3", "#ACE1AF"),  # pastel blue & pastel green
  alpha = 0.55,
  cex = 2.2,
  cat.cex = 2.2,
  main = "Fig 2A – Venn Diagram of expressed bacterial genes"
)

cat("Venn diagram saved to:", file.path(OUTDIR, "Fig2A_Venn.png"), "\n")
