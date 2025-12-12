#!/bin/bash
#SBATCH --job-name=step3_filter_host
#SBATCH --output=/home/christianarturo/project1/logs/step3_filter_host_%j.log
#SBATCH --error=/home/christianarturo/project1/logs/step3_filter_host_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -A introtogds
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=christianarturo@vt.edu

set -euo pipefail

# ======================================================
# PATHS
# ======================================================

PROJECT="/home/christianarturo/project1"

TRIMMED="${PROJECT}/data/trimmed_fastq"
PLANT_GENOME="${PROJECT}/data/plant_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

PLANT_INDEX="${PROJECT}/data/plant_index"
FILTERED_FASTQ="${PROJECT}/data/bacteria_clean_fastq"

mkdir -p "${PLANT_INDEX}" "${FILTERED_FASTQ}"

# ======================================================
# MODULES
# ======================================================

module load Bowtie2/2.5.4-GCC-13.2.0
module load SAMtools/1.18-GCC-12.3.0

# ======================================================
# BUILD BOWTIE2 INDEX (RUN ONCE)
# ======================================================

if [ ! -f "${PLANT_INDEX}/Arabidopsis.1.bt2" ]; then
  echo "Building Bowtie2 index for Arabidopsis..."
  bowtie2-build "${PLANT_GENOME}" "${PLANT_INDEX}/Arabidopsis"
fi

# ======================================================
# HOST REMOVAL (PLANT FILTERING)
# ======================================================

for R1 in ${TRIMMED}/*_1.trim.fastq.gz; do

  SAMPLE=$(basename ${R1} _1.trim.fastq.gz)
  R2="${TRIMMED}/${SAMPLE}_2.trim.fastq.gz"

  echo "-----------------------------------------"
  echo " Filtering host (Arabidopsis) for sample: ${SAMPLE}"
  echo "-----------------------------------------"

  bowtie2 \
    -x "${PLANT_INDEX}/Arabidopsis" \
    -1 "${R1}" \
    -2 "${R2}" \
    -p 12 \
    --very-sensitive \
    -S "${FILTERED_FASTQ}/${SAMPLE}.plant.sam"

  samtools view -b -f 12 \
    "${FILTERED_FASTQ}/${SAMPLE}.plant.sam" \
    > "${FILTERED_FASTQ}/${SAMPLE}.unmapped.bam"

  samtools sort -n \
    "${FILTERED_FASTQ}/${SAMPLE}.unmapped.bam" \
    -o "${FILTERED_FASTQ}/${SAMPLE}.unmapped.sorted.bam"

  samtools fastq \
    -1 "${FILTERED_FASTQ}/${SAMPLE}_bacteria_1.fastq.gz" \
    -2 "${FILTERED_FASTQ}/${SAMPLE}_bacteria_2.fastq.gz" \
    "${FILTERED_FASTQ}/${SAMPLE}.unmapped.sorted.bam"

  echo "Sample ${SAMPLE} filtered (bacterial reads only)"
done

echo "======================================"
echo " STEP 3 COMPLETED"
echo " Bacterial FASTQ files located at:"
echo " ${FILTERED_FASTQ}"
echo "======================================"
