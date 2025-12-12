#!/bin/bash
#SBATCH --job-name=step4_salmon_bacteria
#SBATCH --output=/home/christianarturo/project1/logs/step4_salmon_bacteria_%j.log
#SBATCH --error=/home/christianarturo/project1/logs/step4_salmon_bacteria_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -A introtogds
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=christianarturo@vt.edu

set -euo pipefail

# =====================================================
# PROJECT PATHS
# =====================================================

PROJECT="/home/christianarturo/project1"
DATA="${PROJECT}/data"

CDS_DIR="${DATA}/CDS"
BACT_FASTQ="${DATA}/bacteria_clean_fastq"

SALMON_INDEX="${DATA}/salmon_index_bacteria"
SALMON_QUANT="${DATA}/salmon_quant_bacteria"

REF_CDS="${CDS_DIR}/allCDS_bacteria.ffn"

mkdir -p "${SALMON_INDEX}" "${SALMON_QUANT}"

# =====================================================
# 1. LOAD SALMON MODULE
# =====================================================

module load Salmon/1.10.3-GCC-12.3.0

echo "Using Salmon executable: $(which salmon)"

# =====================================================
# 2. BUILD SALMON INDEX (IF NOT ALREADY PRESENT)
# =====================================================

if [ ! -f "${SALMON_INDEX}/hash.bin" ]; then
  echo "Building Salmon index for bacterial CDS..."

  salmon index \
    -t "${REF_CDS}" \
    -i "${SALMON_INDEX}" \
    -k 31 \
    --threads "${SLURM_CPUS_PER_TASK}"

  echo "Index generated at: ${SALMON_INDEX}"
else
  echo "Salmon index already exists at: ${SALMON_INDEX}"
fi

# =====================================================
# 3. QUANTIFY EACH BACTERIAL SAMPLE
# =====================================================

echo "Starting Salmon quantification for bacterial reads..."

for R1 in ${BACT_FASTQ}/*_bacteria_1.fastq.gz; do
  [ -e "$R1" ] || continue

  SAMPLE=$(basename "${R1}" _bacteria_1.fastq.gz)
  R2="${BACT_FASTQ}/${SAMPLE}_bacteria_2.fastq.gz"

  if [ ! -f "${R2}" ]; then
    echo "Missing R2 file for sample ${SAMPLE}. Skipping."
    continue
  fi

  OUTDIR="${SALMON_QUANT}/${SAMPLE}"
  mkdir -p "${OUTDIR}"

  echo "----------------------------------------------"
  echo " Salmon quantification for sample: ${SAMPLE}"
  echo " R1: ${R1}"
  echo " R2: ${R2}"
  echo " Output directory: ${OUTDIR}"
  echo "----------------------------------------------"

  salmon quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${R1}" \
    -2 "${R2}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    --validateMappings \
    -o "${OUTDIR}"

  echo "Finished processing sample: ${SAMPLE}"
done

echo "==================================================="
echo " STEP 4 COMPLETED (SALMON BACTERIAL QUANTIFICATION)"
echo " Results stored in: ${SALMON_QUANT}"
echo "==================================================="
