#!/bin/bash

# ==============================================================================
# download_hg38_build_index.sh
#
# Author: lishaogang
# Version: 1.1 (Revised for Conda dependency and proper bwa-meth indexing)
#
# Description:
#   This script automates the download and indexing of the hg38 reference
#   genome for the FragPipe workflow. It creates all necessary indexes for
#   BWA (WGS), Samtools (required by the pipeline), and bwameth (WGBS).
#
# Usage:
#   1. Activate your project's Conda environment:
#      conda activate FragPipe-env
#
#   2. Place this script in the project's root directory and run it:
#      ./download_hg38_build_index.sh
#
# ==============================================================================

# --- Script Setup ---
# Exit immediately on error, treat unset variables as errors, and handle pipeline failures.
conda activate cfDNA

set -euo pipefail

# --- Configuration ---
# Get the absolute path of the directory where the script is located.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

# Define paths relative to the script's location.
TARGET_DIR="${SCRIPT_DIR}/workflow/reference"
FINAL_FASTA_NAME="hg38.fa"
FINAL_FASTA_PATH="${TARGET_DIR}/${FINAL_FASTA_NAME}"
BWAMETH_PY="${SCRIPT_DIR}/workflow/scripts/bwameth.py"


# Define the download URL.
UCSC_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
DOWNLOADED_FILE="hg38.fa.gz"

# --- Step 1: Dependency Check ---
echo "INFO: Checking for required tools (curl/wget, bwa, samtools, python3)..."
if ! command -v curl &> /dev/null && ! command -v wget &> /dev/null; then
    echo "ERROR: 'curl' or 'wget' is required for downloading. Please install one and try again." >&2
    exit 1
fi
if ! command -v bwa &> /dev/null; then
    echo "ERROR: 'bwa' not found. Please activate the correct Conda environment (e.g., 'conda activate FragPipe-env')." >&2
    exit 1
fi
if ! command -v samtools &> /dev/null; then
    echo "ERROR: 'samtools' not found. Please activate the correct Conda environment." >&2
    exit 1
fi
if ! command -v python3 &> /dev/null; then
    echo "ERROR: 'python3' not found. Please activate the correct Conda environment." >&2
    exit 1
fi
echo "INFO: All required tools are present."

# Check for pigz for optional parallel decompression.
if command -v pigz &> /dev/null; then
    DECOMPRESS_CMD="pigz -d -c"
    echo "INFO: Using 'pigz' for faster, parallel decompression."
else
    DECOMPRESS_CMD="gzip -d -c"
    echo "INFO: Using standard 'gzip'. For faster decompression, consider installing 'pigz' ('conda install pigz')."
fi

# --- Step 2: Download and Extract Reference Genome ---
mkdir -p "${TARGET_DIR}"
if [[ -f "${FINAL_FASTA_PATH}" ]]; then
    echo "INFO: ${FINAL_FASTA_NAME} already exists at ${TARGET_DIR}. Skipping download and extraction."
else
    echo "--- Starting Download and Extraction ---"
    cd "${TARGET_DIR}" # Temporarily change to the target directory
    
    echo "INFO: Downloading hg38 from ${UCSC_URL}..."
    if command -v curl &> /dev/null; then
        curl -L -o "${DOWNLOADED_FILE}" "${UCSC_URL}"
    else
        wget -O "${DOWNLOADED_FILE}" "${UCSC_URL}"
    fi
    
    echo "INFO: Download complete. Decompressing to ${FINAL_FASTA_NAME}..."
    ${DECOMPRESS_CMD} "${DOWNLOADED_FILE}" > "${FINAL_FASTA_NAME}"
    
    echo "INFO: Cleaning up downloaded archive..."
    rm "${DOWNLOADED_FILE}"
    
    echo "INFO: Download and extraction complete."
    cd "${SCRIPT_DIR}" # Return to the original execution directory
fi

# --- Step 3: Indexing ---

# 3a: BWA Indexing (for WGS)
if [[ -f "${FINAL_FASTA_PATH}.bwt" ]]; then
    echo "INFO: BWA index (.bwt) already exists. Skipping BWA indexing."
else
    echo "--- Starting BWA Indexing (for WGS) ---"
    echo "INFO: This step can take a significant amount of time (e.g., 30-60 minutes)..."
    bwa index "${FINAL_FASTA_PATH}"
    echo "INFO: BWA indexing complete."
fi

# 3b: Samtools Faidx Indexing (Required by the pipeline)
if [[ -f "${FINAL_FASTA_PATH}.fai" ]]; then
    echo "INFO: Samtools FASTA index (.fai) already exists. Skipping."
else
    echo "--- Starting Samtools Faidx Indexing ---"
    samtools faidx "${FINAL_FASTA_PATH}"
    echo "INFO: Samtools faidx indexing complete."
fi

# 3c: Bwameth Indexing (for WGBS)
if [[ -f "${FINAL_FASTA_PATH}.bwameth.c2t.bwt" && -f "${FINAL_FASTA_PATH}.bwameth.g2a.bwt" ]]; then
    echo "INFO: bwa-meth index files already seem to exist. Skipping bwa-meth indexing."
else
    echo "--- Starting bwa-meth Indexing (for WGBS) ---"
    echo "WARNING: This is a very long process and can take several hours."
    echo "INFO: It will create C->T and G->A converted genomes and then index them."
    python "${BWAMETH_PY}" index "${FINAL_FASTA_PATH}"
    echo "INFO: bwa-meth indexing complete."
fi

# --- Final Message ---
echo ""
echo "========================================================"
echo "✅ Reference genome setup for hg38 is complete!"
echo "   The FASTA file and all necessary indexes are now"
echo "   located in: ${TARGET_DIR}"
echo "========================================================"