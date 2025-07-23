#!/bin/bash
# =================================================================================
# STEP 1: QUALITY CONTROL
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qc_env

echo "======================================================"
echo "STEP 1: QUALITY CONTROL"
echo "======================================================"

# Create output directories
mkdir -p "$QC_DIR/fastqc"

# Process samples
for sample_id in "${ALL_SAMPLES[@]}"; do
    # Skip if specific sample requested and this isn't it
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$sample_id" ]]; then
        continue
    fi
    
    echo "Processing sample: $sample_id"
    
    R1_PATH="${READ_PATHS[$sample_id]}/${READ1_FILENAMES[$sample_id]}"
    R2_PATH="${READ_PATHS[$sample_id]}/${READ2_FILENAMES[$sample_id]}"
    SAMPLE_QC_DIR="$QC_DIR/fastqc/$sample_id"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$SAMPLE_QC_DIR/.done" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    mkdir -p "$SAMPLE_QC_DIR"
    
    echo "  Running FastQC..."
    fastqc -t "$THREADS" -o "$SAMPLE_QC_DIR" "$R1_PATH" "$R2_PATH"
    
    # Mark as completed
    touch "$SAMPLE_QC_DIR/.done"
    echo "  Completed"
done

# Generate MultiQC report if all samples processed
if [[ -z "$SAMPLE" ]]; then
    echo "Generating MultiQC report..."
    multiqc "$QC_DIR/fastqc" -o "$QC_DIR" -n "multiqc_report" -f --quiet
fi

conda deactivate

echo "======================================================"
echo "Quality control complete"
echo "Reports available in: $QC_DIR/fastqc"
echo "======================================================"