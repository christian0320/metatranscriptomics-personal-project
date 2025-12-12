#!/bin/bash
#SBATCH --job-name=step1_download
#SBATCH --output=/home/christianarturo/project1/logs/step1_download_%j.log
#SBATCH --error=/home/christianarturo/project1/logs/step1_download_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -A introtogds
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=christianarturo@vt.edu

set -euo pipefail

# ====================================================
# GENERAL CONFIGURATION
# ====================================================

WORKDIR="/home/christianarturo/project1"
cd "$WORKDIR"

mkdir -p logs data/reference data/raw_fastq data/CDS/individual tmp

# ========================================
# 1 DOWNLOAD ARABIDOPSIS GENOME FOR HOST FILTERING
# ========================================

PLANT_GENOME_DIR="${WORKDIR}/data/plant_genome"
PLANT_GENOME="${PLANT_GENOME_DIR}/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"

mkdir -p "${PLANT_GENOME_DIR}"

if [ ! -f "${PLANT_GENOME}" ]; then
  echo ">>> Downloading Arabidopsis TAIR10 genome (Ensembl Plants)..."

  cd "${PLANT_GENOME_DIR}"

  wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

  echo ">>> Unzipping genome..."
  gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

  echo ">>> Genome ready at: ${PLANT_GENOME}"
else
  echo ">>> Arabidopsis genome already exists, skipping download."
fi

# ====================================================
#  2 DOWNLOAD GEO REFERENCE
# ====================================================

echo ">>> Downloading reference GSE231841..."

REF_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/GSE231841/suppl/GSE231841_reference.fasta.gz"

wget -c -O data/reference/GSE231841_reference.fasta.gz "$REF_URL"
gunzip -f data/reference/GSE231841_reference.fasta.gz

echo ">>> Reference saved at: data/reference/GSE231841_reference.fasta"

# ====================================================
#  3 INSTALL SRA TOOLKIT (ONLY IF NOT PRESENT)
# ====================================================

if [ ! -d "$WORKDIR/sratoolkit.3.3.0-ubuntu64" ]; then
  echo ">>> Installing SRA Toolkit..."
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz -O sratoolkit.tar.gz
  tar -xvzf sratoolkit.tar.gz
  rm sratoolkit.tar.gz
fi

SRATOOLKIT_BIN="$WORKDIR/sratoolkit.3.3.0-ubuntu64/bin/fasterq-dump"

# ====================================================
#  4 SRR LIST (RESOLVED IDENTIFIERS)
# ====================================================

declare -a SAMPLES=(
  SRR24611159
  SRR24611160
  SRR24611161
  SRR24611162
  SRR24611163
  SRR24611164
)

# ====================================================
#  5 ROBUST FASTQ DOWNLOAD (WITHOUT PREFETCH)
# ====================================================

for SRR in "${SAMPLES[@]}"; do

  echo "======================================"
  echo " Downloading $SRR using fasterq-dump..."
  echo "======================================"

  $SRATOOLKIT_BIN "$SRR" \
    --outdir data/raw_fastq \
    --threads 8 \
    --temp tmp \
    --split-files \
    --skip-technical

  echo ">>> Compressing FASTQ files for $SRR..."
  gzip data/raw_fastq/${SRR}*.fastq

  echo "$SRR COMPLETED"
done

# ====================================================
#  6 DOWNLOAD BACTERIAL CDS FILES (AT-SPHERE)
# ====================================================

cd data/CDS/individual

echo ">>> Creating download link list for CDS files..."

cat > cds_links.txt <<EOF
http://www.at-sphere.com/download/ORFs_NT/Root100.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root101.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root102.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root105.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root107.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root11.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root112D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1203.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1204.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1212.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1217.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root122.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1220.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1221.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1237.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1238.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root123D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1240.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1252.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1257.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root127.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1272.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1277.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1279.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1280.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1290.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1293.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1294.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1295.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1298.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1304.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root131.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1310.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1312.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1319.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root133.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1334.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root135.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root136.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root137.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root140.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root142.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1423.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1433D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1444.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1455.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1462.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1464.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root147.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1471.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1472.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1480D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1485.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root149.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root1497.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root151.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root154.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root157.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root166.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root16D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root170.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root172.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root179.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root180.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root181.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root186.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root187.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root189.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root190.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root198D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root209.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root214.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root217.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root219.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root22.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root224.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root227.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root231.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root236.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root239.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root240.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root241.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root258.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root264.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root265.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root267.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root268.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root274.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root275.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root278.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root280D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root29.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root31.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root318D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root322.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root329.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root332.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root335.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root336D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root342.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root343.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root344.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root351.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root369.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root381.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root4.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root401.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root402.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root404.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root405.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root411.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root413D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root418.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root420.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root423.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root431.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root434.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root436.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root444D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root456.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root472D3.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root473.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root480.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root482.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root483D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root483D2.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root485.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root487D2Y.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root491.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root494.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root495.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root50.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root52.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root53.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root55.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root552.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root553.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root554.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root558.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root559.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root561.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root562.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root563.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root564.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root565.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root568.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root569.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root60.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root604.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root608.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root61.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root614.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root627.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root63.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root630.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root635.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root65.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root651.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root655.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root656.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root662.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root667.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root66D1.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root670.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root672.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root68.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root682.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root685.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root690.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root695.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root70.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root700.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root708.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root71.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root710.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root720.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root73.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root74.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root76.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root77.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root79.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root81.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root83.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root85.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root9.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root901.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root916.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root918.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root920.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root930.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root935.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root954.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root96.ffn.gz
http://www.at-sphere.com/download/ORFs_NT/Root983.ffn.gz
EOF

echo ">>> Downloading CDS files from AT-Sphere..."

while read -r URL; do
    echo ">>> Downloading $URL"
    curl -L -k -O "$URL"
done < cds_links.txt

gunzip -f *.ffn.gz

# ====================================================
#  7 CONCATENATE ALL CDS FILES
# ====================================================

cat *.ffn > "${WORKDIR}/data/CDS/allCDS_bacteria.ffn"

echo ">>> Total number of CDS files:"
ls *.ffn | wc -l

echo ">>> Final concatenated file:"
ls -lh "${WORKDIR}/data/CDS/allCDS_bacteria.ffn"

echo "======================================"
echo " STEP 1 COMPLETED SUCCESSFULLY"
echo "======================================"
