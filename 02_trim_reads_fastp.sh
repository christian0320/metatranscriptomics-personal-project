#!/bin/bash
#SBATCH --job-name=step2_trim
#SBATCH --output=/home/christianarturo/project1/logs/step2_trim_%j.log
#SBATCH --error=/home/christianarturo/project1/logs/step2_trim_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -A introtogds
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=christianarturo@vt.edu

set -euo pipefail

echo "========================================="
echo " STEP 2 – Read trimming using fastp"
echo "========================================="

PROJECT_DIR="/home/christianarturo/project1"
RAW_DIR="${PROJECT_DIR}/data/raw_fastq"
TRIM_DIR="${PROJECT_DIR}/data/trimmed_fastq"
QC_DIR="${PROJECT_DIR}/data/qc_fastp"

mkdir -p "${TRIM_DIR}" "${QC_DIR}"

module load fastp || true

echo "Files detected in raw_fastq:"
ls ${RAW_DIR}

for R1 in ${RAW_DIR}/*_1.fastq.gz; do

    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    BASENAME=$(basename "${R1/_1.fastq.gz/}")

    echo "-----------------------------------------"
    echo " Processing sample: ${BASENAME}"
    echo "-----------------------------------------"

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${TRIM_DIR}/${BASENAME}_1.trim.fastq.gz" \
        -O "${TRIM_DIR}/${BASENAME}_2.trim.fastq.gz" \
        --thread ${SLURM_CPUS_PER_TASK} \
        --detect_adapter_for_pe \
        --length_required 50 \
        --cut_mean_quality 20 \
        --json "${QC_DIR}/${BASENAME}.fastp.json" \
        --html "${QC_DIR}/${BASENAME}.fastp.html"

    echo "Sample ${BASENAME} completed"
done

echo "========================================="
echo " STEP 2 COMPLETED"
echo " Trimmed FASTQ files: ${TRIM_DIR}"
echo " QC reports: ${QC_DIR}"
echo "========================================="
