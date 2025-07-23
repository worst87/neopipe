#!/bin/bash
# =================================================================================
# STEP 3: EXPRESSION ANALYSIS (RNA-SPECIFIC)
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"
HELPER_DIR="$(dirname "$0")/helpers"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate analysis_env

echo "======================================================"
echo "STEP 3: EXPRESSION ANALYSIS"
echo "======================================================"

# Build Kallisto index
KALLISTO_INDEX="$RESOURCE_DIR/kallisto_index/transcriptome.idx"
mkdir -p "$(dirname "$KALLISTO_INDEX")"

if [[ ! -f "$KALLISTO_INDEX" ]]; then
    echo "Building Kallisto index..."
    kallisto index -i "$KALLISTO_INDEX" "$REF_TRANSCRIPTOME_FA"
fi

# Create output directories
mkdir -p "$EXPR_DIR/kallisto" "$EXPR_DIR/gene_tpm"

# Process each sample
for sample_id in "${ALL_SAMPLES[@]}"; do
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$sample_id" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing sample: $sample_id"
    
    KALLISTO_OUTPUT="$EXPR_DIR/kallisto/$sample_id"
    GENE_TPM_FILE="$EXPR_DIR/gene_tpm/${sample_id}_gene_tpm.tsv"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$GENE_TPM_FILE" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    # Get read paths
    R1_PATH="${READ_PATHS[$sample_id]}/${READ1_FILENAMES[$sample_id]}"
    R2_PATH="${READ_PATHS[$sample_id]}/${READ2_FILENAMES[$sample_id]}"
    
    # Step 1: Kallisto quantification
    echo "  [1/2] Running Kallisto quantification..."
    mkdir -p "$KALLISTO_OUTPUT"
    
    kallisto quant \
        -i "$KALLISTO_INDEX" \
        -o "$KALLISTO_OUTPUT" \
        --threads="$THREADS" \
        -b "$KALLISTO_BOOTSTRAPS" \
        "$R1_PATH" "$R2_PATH"
    
    # Step 2: Summarize to gene level
    echo "  [2/2] Summarizing to gene-level TPM..."
    Rscript "$HELPER_DIR/summarize_kallisto.R" \
        --gtf "$REF_GENOME_GTF" \
        --kallisto_dir "$KALLISTO_OUTPUT" \
        --sample_name "$sample_id" \
        --output_dir "$EXPR_DIR/gene_tpm"
    
    echo "  Completed expression analysis for $sample_id"
done

conda deactivate

echo ""
echo "======================================================"
echo "Expression analysis complete"
echo "Gene-level TPMs: $EXPR_DIR/gene_tpm"
echo "======================================================"