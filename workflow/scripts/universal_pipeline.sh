#!/bin/bash
#
# ==============================================================================
# Script Name: universal_pipeline.sh
#
# Description:
#   A fully automated, zero-configuration pipeline for WGS/WGBS data analysis.
#   It runs a fixed workflow and automatically discovers ALL dependencies
#   (tools, scripts, reference files) from a standard directory structure.
#
# Assumed Directory Structure:
#   /project_root/
#   ├── scripts/  (this script and all other tools/scripts)
#   └── reference/ (hg38.fa, hg38.tss.bed, hg38.ctcf.bed, etc.)
#
# Author:
#   Shaogang Li (Merged and Enhanced by AI)
# ==============================================================================

# --- Strict Mode & Error Handling ---
set -e
set -o pipefail

# --- Activate Conda Environment ---
if ! command -v conda &> /dev/null; then echo "Error: 'conda' not found." >&2; exit 1; fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cfDNA

# --- Function Definitions ---
usage() {
    local script_name
    script_name=$(basename "$0")
    echo "Usage: $script_name -y <type> -1 <fastq1> -2 <fastq2> -s <sample_id> -p <output_path> -l <label> [options]"
    echo ""
    echo "This script runs a FIXED pipeline. It automatically finds all required reference files."
    echo ""
    echo "Core Required Arguments:"
    echo "  -y, --type           Analysis type. Must be 'wgs' or 'wgbs'."
    echo "  -1, --fastq1         Path to forward-strand FASTQ file (R1)."
    echo "  -2, --fastq2         Path to reverse-strand FASTQ file (R2)."
    echo "  -s, --sample-id      Unique Sample Identifier."
    echo "  -p, --output-path    Directory where all final result files will be saved."
    echo "  -l, --label          Sample classification label (e.g., Cancer, Normal)."
    echo ""
    echo "General Options:"
    echo "  -t, --threads        Number of threads to use. [Default: 8]"
    echo "  --keep-intermediate  Keep all intermediate files."
    echo "  -h, --help           Display this help message."
    exit 1
}

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] - $1"; }
cleanup() { if [[ "$KEEP_INTERMEDIATE" = false && -d "$TMP_DIR" ]]; then log "Cleaning up temporary directory: $TMP_DIR"; rm -rf "$TMP_DIR"; log "Cleanup finished."; fi; }

# --- Configuration & Default values ---
THREADS=8
KEEP_INTERMEDIATE=false

# --- Automatic Path Discovery ---
# This block now defines ALL paths to tools and reference files automatically.
# No path arguments are needed from the user for standard runs.
log "Discovering dependency paths..."
SCRIPT_DIR=$(cd -P "$(dirname "$0")" && pwd)
PROJECT_ROOT_DIR=$(dirname "$SCRIPT_DIR")

# Tool Paths (assuming they are in the 'scripts' directory)
# <--- 你问的第一个代码块现在只负责定义路径，不再需要用户输入 --->
FASTP_CMD="$SCRIPT_DIR/fastp"
PICARD_JAR="$SCRIPT_DIR/picard.jar"
BWAMETH_CMD="$SCRIPT_DIR/bwameth.py"
BAM_TO_FRAG_CMD="$SCRIPT_DIR/bam_to_frag_bed.py"
PROCESS_FRAG_BED_CMD="$SCRIPT_DIR/process_fragment_bed.py"
TSS_R_SCRIPT="${SCRIPT_DIR}/TSS_normaolize.R"
CTCF_R_SCRIPT="${SCRIPT_DIR}/CTCF_normaolize.R"

# Reference File Paths (assuming they are in the '../reference' directory)
REF_DIR="${PROJECT_ROOT_DIR}/reference"
REF_GENOME="${REF_DIR}/hg38.fa" # <--- NEW: 主参考基因组也自动寻找
TSS_REF_BED="${REF_DIR}/hg38.tss.UD2500bp.txt"
CTCF_REF_BED="${REF_DIR}/CTCF_binding_sites_UD1500_hg38.bed"

# --- Parse Command-line Arguments ---
# <--- 你问的第二个代码块现在变得非常简洁，因为它不再需要处理路径参数 --->
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -y|--type) ANALYSIS_TYPE="$2"; shift ;;
        -1|--fastq1) FASTQ1="$2"; shift ;;
        -2|--fastq2) FASTQ2="$2"; shift ;;
        -s|--sample-id) SAMPLE_ID="$2"; shift ;;
        -p|--output-path) OUTPUT_PATH="$2"; shift ;;
        -l|--label) LABEL="$2"; shift ;;
        -t|--threads) THREADS="$2"; shift ;;
        --keep-intermediate) KEEP_INTERMEDIATE=true ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# --- Validate ALL Paths and Inputs ---
log "Validating all inputs and auto-discovered dependencies..."
if [ -z "$ANALYSIS_TYPE" ] || [ -z "$FASTQ1" ] || [ -z "$FASTQ2" ] || [ -z "$SAMPLE_ID" ] || [ -z "$OUTPUT_PATH" ] || [ -z "$LABEL" ]; then
    log "Error: Missing one or more required arguments." >&2; usage;
fi

# Validate existence of all mandatory files and tools
for f in \
    "$FASTQ1" "$FASTQ2" \
    "$FASTP_CMD" "$PICARD_JAR" "$BWAMETH_CMD" "$BAM_TO_FRAG_CMD" "$PROCESS_FRAG_BED_CMD" \
    "$REF_GENOME" "$TSS_REF_BED" "$CTCF_REF_BED" \
    "$TSS_R_SCRIPT" "$CTCF_R_SCRIPT"; do
    if [ ! -e "$f" ]; then log "Error: Required file or script not found at its expected location: '$f'" >&2; exit 1; fi
done

log "All inputs and dependencies validated successfully."

# --- Setup Environment & Auto-generate Paths ---
log "Pipeline started: $(basename "$0")"
mkdir -p "$OUTPUT_PATH"
mkdir -p "${OUTPUT_PATH}/${LABEL}" # Subdirectory for R script outputs

OUTPUT_BAM="${OUTPUT_PATH}/${SAMPLE_ID}.markdup.bam"
OUTPUT_BED="${OUTPUT_PATH}/${SAMPLE_ID}.fragment.bed"
OUTPUT_PREFIX="${OUTPUT_PATH}/${SAMPLE_ID}"

log "Sample ID: $SAMPLE_ID | Label/Group: $LABEL"

TMP_DIR=$(mktemp -d -p "$OUTPUT_PATH" "tmp_${SAMPLE_ID}_XXXXXXXX")
trap cleanup EXIT

# --- Main Pipeline ---
log "Step 1-6: Running core pipeline from FASTQ to Fragment BED..."
# [ Code for steps 1-6 remains here, unchanged ]
log "Step 1-6: Core pipeline finished. Fragment BED created at ${OUTPUT_BED}"

log "Step 7: Running motif and size analysis..."
PYTHON_CMD_ARGS=("--fragment-bed" "$OUTPUT_BED" "--ref-fasta" "$REF_GENOME" "--output-prefix" "$OUTPUT_PREFIX" "--label" "$LABEL")
if [[ "$KEEP_INTERMEDIATE" = true ]]; then PYTHON_CMD_ARGS+=("--keep-intermediate"); fi
python3 "$PROCESS_FRAG_BED_CMD" "${PYTHON_CMD_ARGS[@]}"
log "Step 7: Motif and size analysis completed."

log "Step 8: Running TSS enrichment analysis..."
log "  - Using TSS reference: $TSS_REF_BED"
log "  - Using TSS R script: $TSS_R_SCRIPT"
TSS_COVERAGE_FILE="${TMP_DIR}/${SAMPLE_ID}_TSS.coverage.txt"
R_OUT_DIR="${OUTPUT_PATH}/${LABEL}"
grep -E '^chr([0-9]{1,2}|X|Y)\b' "$OUTPUT_BED" | bedtools coverage -a "$TSS_REF_BED" -b - -d > "$TSS_COVERAGE_FILE"
Rscript "$TSS_R_SCRIPT" "$SAMPLE_ID" "$LABEL" "$TSS_COVERAGE_FILE" "$R_OUT_DIR"
log "Step 8: TSS analysis finished. Results are in $R_OUT_DIR"

log "Step 9: Running CTCF enrichment analysis..."
log "  - Using CTCF reference: $CTCF_REF_BED"
log "  - Using CTCF R script: $CTCF_R_SCRIPT"
CTCF_COVERAGE_FILE="${TMP_DIR}/${SAMPLE_ID}_CTCF.coverage.txt"
grep -E '^chr([0-9]{1,2}|X|Y)\b' "$OUTPUT_BED" | bedtools coverage -a "$CTCF_REF_BED" -b - -d > "$CTCF_COVERAGE_FILE"
Rscript "$CTCF_R_SCRIPT" "$SAMPLE_ID" "$LABEL" "$CTCF_COVERAGE_FILE" "$R_OUT_DIR"
log "Step 9: CTCF analysis finished. Results are in $R_OUT_DIR"

### Step 10: Window Protection Score (WPS) Analysis ###
log "Step 10: Calculating Window Protection Score (WPS)..."
log "  - Using WPS R script: $WPS_R_SCRIPT"
log "  - Using WPS gene reference: $WPS_GENE_BED"

# The R script creates its own output directories, so we just provide the main path.
Rscript "$WPS_R_SCRIPT" \
    --task calculate \
    --mode Multi_group \
    --Fragment "$OUTPUT_BED" \
    --genebed "$WPS_GENE_BED" \
    --prefix "$SAMPLE_ID" \
    --outdir "$OUTPUT_PATH" \
    --bedtools "bedtools" \
    --cpus "$THREADS" \
    --group "$LABEL" \
    --plot FALSE \
    --tss 230000000 # Using the default from your example command

log "Step 10: WPS analysis finished."

# --- Finalization ---
log "Pipeline successfully completed."
echo "================================================================================"
echo "Final outputs for sample '${SAMPLE_ID}' are located in: ${OUTPUT_PATH}"
echo "  - BAM File: ${OUTPUT_BAM}"
echo "  - Size Ratio File: ${OUTPUT_PREFIX}.size_ratio"
echo "  - MDS File: ${OUTPUT_PREFIX}.4mer.motif.MDS"
echo "  - TSS/CTCF Analysis Results: in directory ${OUTPUT_PATH}/${LABEL}/"
echo "  - WPS Analysis Results: in directory ${OUTPUT_PATH}/Result/${LABEL}/${SAMPLE_ID}/"
if [[ "$KEEP_INTERMEDIATE" = true ]]; then
    echo "Intermediate files were kept."
else
    echo "(Note: Intermediate files were removed.)"
fi
echo "================================================================================"