
---

## Pipeline Description

### Step 1 – Data and Reference Download
**Script:** `01_download_data.sh`

- Downloads the GEO reference genome for dataset **GSE231841**
- Downloads raw RNA-seq reads from SRA using `fasterq-dump`
- Retrieves bacterial coding sequences (CDS) from **AT-SPHERE**
- Decompresses and concatenates all bacterial CDS into a single reference file

**Output:**
- Raw FASTQ files  
- Bacterial CDS reference (`allCDS_bacteria.ffn`)

---

### Step 2 – Read Trimming and Quality Control
**Script:** `02_trim_reads_fastp.sh`

- Adapter trimming and quality filtering using **fastp**
- Removal of low-quality bases and short reads
- Generation of QC reports

**Output:**
- Quality-filtered FASTQ files

---

### Step 3 – Host Read Removal
**Script:** `03_filter_host_bowtie2.sh`

- Alignment of reads to the plant host genome using **Bowtie2**
- Removal of host-derived reads
- Retention of microbial reads only

**Output:**
- Host-filtered FASTQ files

---

### Step 4 – Bacterial Gene Expression Quantification
**Script:** `04_salmon_bacteria.sh`

- Construction of a Salmon index using bacterial CDS
- Transcript-level quantification of gene expression
- TPM-based abundance estimates

**Output:**
- Salmon quantification directories per sample

---

### Step 4.5 – Gene-to-Species Mapping
**Script:** `04.5_build_gene_to_species.py`

- Mapping of bacterial genes to their corresponding species
- Generation of a gene-to-species lookup table

**Output:**
- Gene-to-species mapping file

---

### Step 5 – TPM Matrix Construction
**Script:** `05_build_tpm_matrix.py`

- Extraction of TPM values from Salmon outputs
- Merging of all samples into a unified TPM matrix
- Preparation of expression data for downstream analysis

**Output:**
- TPM expression matrix (genes × samples)

---

### Step 6 – Relative Abundance Calculation
**Script:** `06_relative_abundance_bacteria.py`

- Aggregation of gene-level TPM values by species
- Normalization to relative abundance
- Separation of samples by experimental condition (root vs matrix)

**Output:**
- Species-level relative abundance table

---

### Step 7 – Venn Diagram (Figure 2A)
**Script:** `07_fig2A_venn.R`

- Identification of shared and unique bacterial species between conditions
- Generation of a Venn diagram reproducing **Figure 2A**

**Output:**
- Venn diagram figure

---

### Step 8 – Relative Abundance Visualization (Figure 2C)
**Script:** `08_fig2C_relative_abundance.R`

- Visualization of species-level relative abundance
- Reproduction of **Figure 2C** from the reference study

**Output:**
- Relative abundance plots

---

## Requirements

### Software
- Bash
- Python ≥ 3.8
- R ≥ 4.2
- fastp
- Bowtie2
- Salmon
- SRA Toolkit

### R Packages
- tidyverse
- VennDiagram
- ggplot2

---

## Usage

Scripts should be executed **sequentially**, following their numbering:

```bash
sbatch 01_download_data.sh
sbatch 02_trim_reads_fastp.sh
sbatch 03_filter_host_bowtie2.sh
sbatch 04_salmon_bacteria.sh
python 04.5_build_gene_to_species.py
python 05_build_tpm_matrix.py
python 06_relative_abundance_bacteria.py
Rscript 07_fig2A_venn.R
Rscript 08_fig2C_relative_abundance.R
