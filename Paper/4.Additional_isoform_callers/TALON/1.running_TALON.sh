#!/bin/bash
#SBATCH --job-name=talon_batches_MD_label_resumable
#SBATCH --account=stergachislab
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=120:00:00
#SBATCH --mem=400G
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --chdir=/mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker_testing/TALON/8.30.25_allow_resume
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yhhc@uw.edu
#SBATCH --export=ALL

###########################
# Env
###########################
source /mmfs1/gscratch/stergachislab/yhhc/tools/miniconda3/miniconda3/bin/activate
conda activate talon

set -euo pipefail
trap 'echo "FAILED at line $LINENO" >&2' ERR
log(){ echo "[$(date +'%F %T')] $*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }

###########################
# User settings
###########################
CONFIG_ALL="/mmfs1/gscratch/stergachislab/yhhc/projects/IsoRanker_testing/TALON/8.30.25_allow_resume/configV2.csv"

OUT_PREFIX="talon_run"
DB="${OUT_PREFIX}.db"
BUILD_NAME="hg38"
ANNOT_NAME="gencode_v47"
ANNOT_GTF="/mmfs1/gscratch/stergachislab/yhhc/tools/Gencode/gencode.v47.annotation.gtf"
GENOME_FA="/gscratch/stergachislab/assemblies/simple-names/hg38.fa"

THREADS=7
BATCH_SIZE=2

###########################
# Dirs
###########################
LOGDIR="logsV2"
OUTDIR="outputs"
CHUNK_DIR="${OUTDIR}/chunks"
MD_DIR="md_bams"
TMP_SAM_DIR="tmp_sams"
LABEL_SAM_DIR="labeled_sams"
LABEL_BAM_DIR="labeled_bams"

mkdir -p "$LOGDIR" "$OUTDIR" "$CHUNK_DIR" "$MD_DIR" "$TMP_SAM_DIR" "$LABEL_SAM_DIR" "$LABEL_BAM_DIR"

###########################
# Checks
###########################
[[ -s "$CONFIG_ALL" ]] || die "Missing config: $CONFIG_ALL"
command -v samtools >/dev/null 2>&1 || die "samtools not found"
command -v talon >/dev/null 2>&1 || die "talon not found"
command -v talon_label_reads >/dev/null 2>&1 || die "talon_label_reads not found"
[[ -s "$ANNOT_GTF" ]] || die "Missing GTF: $ANNOT_GTF"
[[ -s "$GENOME_FA" ]] || die "Missing FASTA: $GENOME_FA"
[[ -s "${GENOME_FA}.fai" ]] || { log "Indexing $GENOME_FA"; samtools faidx "$GENOME_FA"; }

###########################
# Helpers
###########################
make_md_bam(){
  local inpath="$1" dataset="$2"
  local outbam="${MD_DIR}/${dataset}.md.bam"
  local clog="${LOGDIR}/calmd.${dataset}.log"

  case "$inpath" in
    *.bam|*.BAM)
      log "[MD] ${dataset}: BAM -> calmd"
      samtools calmd -bAr "$inpath" "$GENOME_FA" 2> "$clog" > "$outbam"
      ;;
    *.sam|*.SAM)
      log "[MD] ${dataset}: SAM -> calmd"
      samtools calmd -bAr "$inpath" "$GENOME_FA" 2> "$clog" > "$outbam"
      ;;
    *.cram|*.CRAM)
      log "[MD] ${dataset}: CRAM -> BAM + calmd"
      samtools view -b -T "$GENOME_FA" -@ "$THREADS" "$inpath" \
        | samtools calmd -bAr - "$GENOME_FA" 2> "$clog" > "$outbam"
      ;;
    *) die "Unsupported alignment format for $dataset: $inpath" ;;
  esac

  [[ -s "$outbam" ]] || die "MD BAM not created: $outbam"
  samtools quickcheck -v "$outbam" || die "samtools quickcheck failed: $outbam"
  samtools index -@ "$THREADS" "$outbam"
  echo "$outbam"
}

label_to_bam(){
  local dataset="$1" md_bam="$2"
  local fifo_in="${TMP_SAM_DIR}/${dataset}.fifo.sam"
  local lbl_prefix="${LABEL_SAM_DIR}/${dataset}"
  local lbl_sam="${lbl_prefix}_labeled.sam"
  local lbl_log="${LOGDIR}/labelreads.${dataset}.log"
  local lbam="${LABEL_BAM_DIR}/${dataset}.labeled.bam"

  rm -f "$fifo_in" "$lbl_sam" "$lbam"
  mkfifo "$fifo_in"

  (
  set -euo pipefail
  talon_label_reads \
    --f "$fifo_in" \
    --g "$GENOME_FA" \
    --t "$THREADS" \
    --deleteTmp \
    --o "$lbl_prefix" >> "$lbl_log" 2>&1
) & pid_label=$!

  samtools view -@ "$THREADS" -h -F 0x904 "$md_bam" > "$fifo_in" || true
  wait "$pid_label" || { tail -n 80 "$lbl_log" || true; die "talon_label_reads failed for $dataset"; }
  rm -f "$fifo_in"

  [[ -s "$lbl_sam" ]] || { tail -n 80 "$lbl_log" || true; die "Missing labeled SAM: $lbl_sam"; }

  samtools view -@ "$THREADS" -b "$lbl_sam" | samtools sort -@ "$THREADS" -o "$lbam" -
  samtools index -@ "$THREADS" "$lbam"
  rm -f "$lbl_sam"

  echo "$lbam"
}

###########################
# Split config into chunks
###########################
log "[*] Splitting $CONFIG_ALL into chunks of $BATCH_SIZE lines"
rm -f "${CHUNK_DIR}/batch_"*.csv 2>/dev/null || true
split -d -l "$BATCH_SIZE" -a 3 --additional-suffix=.csv "$CONFIG_ALL" "${CHUNK_DIR}/batch_"
mapfile -t CHUNKS < <(ls -1 "${CHUNK_DIR}"/batch_*.csv | sort -V)
[[ ${#CHUNKS[@]} -gt 0 ]] || die "No chunks created"

# ---------- init DB ----------
if [[ -s "$DB" ]]; then
  log "[*] DB exists: $DB (skip talon_initialize_database)"
else
  log "[*] Initializing TALON DB: $DB"
  talon_initialize_database \
    --f "$ANNOT_GTF" \
    --a "$ANNOT_NAME" \
    --g "$BUILD_NAME" \
    --o "$OUT_PREFIX"
  [[ -s "$DB" ]] || die "DB initialization failed: $DB not created"
fi

###########################
# Completed batches tracker
###########################
DONE_FILE="${OUTDIR}/completed_batches.txt"
touch "$DONE_FILE"

###########################
# Process each chunk
###########################
for CHUNK in "${CHUNKS[@]}"; do
  base="$(basename "$CHUNK" .csv)"   # e.g., batch_000

  # Skip if already completed. Very important, because allows resuming of this script
  if grep -qx "$base" "$DONE_FILE"; then
    log "[*] ${base}: already completed, skipping"
    continue
  fi

  chunk_md_csv="${CHUNK_DIR}/${base}.withMD.csv"
  chunk_bam_cfg="${CHUNK_DIR}/${base}.labeledBAM.csv"
  batch_prefix="${OUT_PREFIX}.${base}"

  log "[*] ${base}: Making MD-tagged BAMs -> $(basename "$chunk_md_csv")"
  : > "$chunk_md_csv"
  while IFS=, read -r DATASET SAMPLE PLATFORM LOCATION || [[ -n "${DATASET:-}" ]]; do
    [[ -z "${DATASET:-}" || -z "${LOCATION:-}" ]] && continue
    [[ -s "$LOCATION" ]] || die "File not found for $DATASET: $LOCATION"
    md_bam="$(make_md_bam "$LOCATION" "$DATASET")"
    echo "${DATASET},${SAMPLE},${PLATFORM},${md_bam}" >> "$chunk_md_csv"
  done < <(sed 's/\r$//' "$CHUNK")
  [[ -s "$chunk_md_csv" ]] || die "${base}: empty MD CSV"

  log "[*] ${base}: Labeling (adds fA) -> $(basename "$chunk_bam_cfg")"
  : > "$chunk_bam_cfg"
  while IFS=, read -r DATASET SAMPLE PLATFORM MD_BAM || [[ -n "${DATASET:-}" ]]; do
    [[ -z "${DATASET:-}" || -z "${MD_BAM:-}" ]] && continue
    lbam="$(label_to_bam "$DATASET" "$MD_BAM")"
    echo "${DATASET},${SAMPLE},${PLATFORM},${lbam}" >> "$chunk_bam_cfg"
  done < "$chunk_md_csv"
  [[ -s "$chunk_bam_cfg" ]] || die "${base}: empty labeled BAM CSV"

  log "[*] ${base}: TALON annotate (threads=$THREADS)"
  talon \
    --f "$chunk_bam_cfg" \
    --db "$DB" \
    --build "$BUILD_NAME" \
    --o "$batch_prefix" \
    --threads "$THREADS" \
    |& tee "${LOGDIR}/talon.annotate.${base}.log"

  # Mark batch as completed
  echo "$base" >> "$DONE_FILE"
done

###########################
# Abundance and filtering
###########################
log "[*] talon_abundance (unfiltered, FINAL)..."
talon_abundance \
  --db "$DB" \
  -a "$ANNOT_NAME" \
  --build "$BUILD_NAME" \
  --o "${OUT_PREFIX}.unfiltered" \
  |& tee "$LOGDIR/talon.abundance.unfiltered.log"

MIN_DATASETS=1
MIN_COUNT=5
MAX_FRAC_A=0.5
PASS_LIST="${OUT_PREFIX}_filtered_transcripts.csv"

log "[*] talon_filter_transcripts ..."
talon_filter_transcripts \
  --db "$DB" \
  -a "$ANNOT_NAME" \
  --maxFracA "$MAX_FRAC_A" \
  --minCount "$MIN_COUNT" \
  --minDatasets "$MIN_DATASETS" \
  --o "$PASS_LIST" \
  |& tee "$LOGDIR/talon.filter.log"

log "[*] talon_abundance (filtered, whitelist=$PASS_LIST)..."
talon_abundance \
  --db "$DB" \
  --whitelist "$PASS_LIST" \
  -a "$ANNOT_NAME" \
  --build "$BUILD_NAME" \
  --o "${OUT_PREFIX}.filtered" \
  |& tee "$LOGDIR/talon.abundance.filtered.log"

log "[*] ALL DONE."
