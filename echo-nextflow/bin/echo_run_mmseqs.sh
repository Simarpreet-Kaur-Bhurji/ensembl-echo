#!/usr/bin/env bash
set -euo pipefail

# Wrapper for MMseqs2 easy-cluster via Singularity
# Produces: mmseqs_results_cluster.tsv

usage() {
  echo "Usage: $0 --input_fasta <file> --out_prefix <prefix> --tmp_dir <dir> --singularity_image <img> --min_seq_id <float> --coverage <float> --cov_mode <int> --threads <int>"
  exit 1
}

INPUT_FASTA=""
OUT_PREFIX=""
TMP_DIR=""
SIF=""
MIN_SEQ_ID=""
COVERAGE=""
COV_MODE=""
THREADS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_fasta) INPUT_FASTA="$2"; shift 2 ;;
    --out_prefix) OUT_PREFIX="$2"; shift 2 ;;
    --tmp_dir) TMP_DIR="$2"; shift 2 ;;
    --singularity_image) SIF="$2"; shift 2 ;;
    --min_seq_id) MIN_SEQ_ID="$2"; shift 2 ;;
    --coverage) COVERAGE="$2"; shift 2 ;;
    --cov_mode) COV_MODE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; usage ;;
  esac
done

[[ -z "$INPUT_FASTA" || -z "$OUT_PREFIX" || -z "$TMP_DIR" || -z "$SIF" || -z "$MIN_SEQ_ID" || -z "$COVERAGE" || -z "$COV_MODE" || -z "$THREADS" ]] && usage

mkdir -p "$TMP_DIR"

echo "[INFO] Running MMseqs2 easy-cluster"
echo "[INFO] input_fasta      = $INPUT_FASTA"
echo "[INFO] out_prefix       = $OUT_PREFIX"
echo "[INFO] tmp_dir          = $TMP_DIR"
echo "[INFO] singularity_image= $SIF"
echo "[INFO] min_seq_id       = $MIN_SEQ_ID"
echo "[INFO] coverage         = $COVERAGE"
echo "[INFO] cov_mode         = $COV_MODE"
echo "[INFO] threads          = $THREADS"

singularity exec "$SIF" \
  mmseqs easy-cluster \
  "$INPUT_FASTA" \
  "$OUT_PREFIX" \
  "$TMP_DIR" \
  --min-seq-id "$MIN_SEQ_ID" \
  -c "$COVERAGE" \
  --cov-mode "$COV_MODE" \
  --threads "$THREADS"

echo "[OK] wrote mmseqs_results_cluster.tsv"

