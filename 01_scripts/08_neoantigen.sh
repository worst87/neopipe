#!/bin/bash
# =================================================================================
# STEP 8: NEOANTIGEN PREDICTION
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pvacseq_env

echo "======================================================"
echo "STEP 8: NEOANTIGEN PREDICTION WITH PVACSEQ"
echo "======================================================"

for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    NORMAL_SAMPLE="${TUMOR_NORMAL_PAIRS[$TUMOR_SAMPLE]}"
    
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "=========================================="
    echo "Processing: $TUMOR_SAMPLE"
    echo "=========================================="
    
    # Define paths
    SAMPLE_VAR_DIR="$VAR_DIR/$TUMOR_SAMPLE"
    INPUT_VCF="$SAMPLE_VAR_DIR/annotated/${TUMOR_SAMPLE}.final.annotated.vcf.gz"
    PHASED_VCF="$SAMPLE_VAR_DIR/phased/${TUMOR_SAMPLE}.phased.annotated.vcf.gz"
    PVAC_OUTPUT="$NEO_DIR/$TUMOR_SAMPLE/pvacseq"
    REPORTS_DIR="$NEO_DIR/$TUMOR_SAMPLE/reports"
    
    # Validate input
    if [[ ! -f "$INPUT_VCF" ]]; then
        echo "ERROR: Annotated VCF not found: $INPUT_VCF"
        continue
    fi
    
    # Check variant counts
    MAIN_COUNT=$(bcftools view -H "$INPUT_VCF" | wc -l)
    echo "  Main VCF variants: $MAIN_COUNT"
    
    USE_PHASED=false
    if [[ -s "$PHASED_VCF" ]]; then
        PHASED_COUNT=$(bcftools view -H "$PHASED_VCF" | wc -l)
        echo "  Phased VCF variants: $PHASED_COUNT"
        if [[ "$PHASED_COUNT" -gt 0 ]]; then
            USE_PHASED=true
            echo "  Will use phased VCF for proximal variant correction"
        fi
    else
        echo "  No phased VCF available - running without proximal variant correction"
    fi
    
    # Clean output directory
    rm -rf "$PVAC_OUTPUT"
    mkdir -p "$PVAC_OUTPUT" "$REPORTS_DIR"
    
    # Build pVACseq command
    echo ""
    echo "  Starting pVACseq run..."
    
    PVACSEQ_CMD="pvacseq run \
        \"$INPUT_VCF\" \
        \"$TUMOR_SAMPLE\" \
        \"$MHC_ALLELES\" \
        $PVAC_ALGORITHMS \
        \"$PVAC_OUTPUT\" \
        --normal-sample-name \"$NORMAL_SAMPLE\" \
        -e1 \"$PEPTIDE_LENGTHS\" \
        -b \"$BINDING_THRESHOLD\" \
        --top-score-metric median \
        --keep-tmp-files \
        -t \"$PVACSEQ_THREADS\" \
        --expn-val \"$EXPRESSION_CUTOFF\" \
        --trna-vaf \"$TUMOR_RNA_VAF\" \
        --trna-cov \"$TUMOR_RNA_COV\" \
        --normal-vaf \"$NORMAL_VAF\" \
        --normal-cov \"$NORMAL_COV\" \
        --maximum-transcript-support-level \"$MAX_TSL\""
    
    if [[ "$USE_PHASED" == true ]]; then
        PVACSEQ_CMD="${PVACSEQ_CMD} --phased-proximal-variants-vcf \"$PHASED_VCF\""
    fi
    
    # Run pVACseq
    eval $PVACSEQ_CMD
    
    # Generate aggregated report
    if [[ -d "$PVAC_OUTPUT/MHC_Class_I" ]]; then
        echo ""
        echo "  Generating aggregated report..."
        RESULTS_TSV="$PVAC_OUTPUT/MHC_Class_I/${TUMOR_SAMPLE}.all_epitopes.tsv"
        
        if [[ -f "$RESULTS_TSV" ]]; then
            AGGREGATED_TSV="$REPORTS_DIR/${TUMOR_SAMPLE}.all_epitopes.aggregated.tsv"
            AGGREGATED_JSON="$REPORTS_DIR/${TUMOR_SAMPLE}.all_epitopes.aggregated.metrics.json"
            
            pvacseq generate_aggregated_report "$RESULTS_TSV" "$AGGREGATED_TSV"
            
            # Copy for convenience
            cp "$PVAC_OUTPUT/MHC_Class_I/${TUMOR_SAMPLE}.all_epitopes.aggregated.metrics.json" "$AGGREGATED_JSON" 2>/dev/null || true
            
            # Count results
            EPITOPE_COUNT=$(tail -n +2 "$RESULTS_TSV" | wc -l)
            FILTERED_COUNT=$(tail -n +2 "$AGGREGATED_TSV" | wc -l)
            
            echo "  Total epitopes predicted: $EPITOPE_COUNT"
            echo "  Filtered epitopes: $FILTERED_COUNT"
        else
            echo "  WARNING: No epitopes predicted"
        fi
    else
        echo "  ERROR: pVACseq output structure not found"
    fi
    
    echo ""
    echo "  pVACseq complete for $TUMOR_SAMPLE"
done

conda deactivate

echo ""
echo "======================================================"
echo "Neoantigen prediction complete"
echo "Results: $NEO_DIR/*/reports"
echo "======================================================"