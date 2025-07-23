#!/bin/bash
# =================================================================================
# STEP 9: INTERACTIVE VISUALIZATION
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pvacdownstream_env

echo "======================================================"
echo "STEP 9: INTERACTIVE ANALYSIS WITH PVACVIEW"
echo "======================================================"

for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "Preparing pVACview for: $TUMOR_SAMPLE"
    echo "====================================="
    
    # Define paths
    PVAC_OUTPUT="$NEO_DIR/$TUMOR_SAMPLE/pvacseq"
    REPORTS_DIR="$NEO_DIR/$TUMOR_SAMPLE/reports"
    
    AGGREGATED_TSV="$REPORTS_DIR/${TUMOR_SAMPLE}.all_epitopes.aggregated.tsv"
    AGGREGATED_JSON="$REPORTS_DIR/${TUMOR_SAMPLE}.all_epitopes.aggregated.metrics.json"
    
    # Check required files
    if [[ ! -f "$AGGREGATED_TSV" ]] || [[ ! -f "$AGGREGATED_JSON" ]]; then
        echo "ERROR: Required pVACview files not found!"
        echo "  Missing TSV: $AGGREGATED_TSV"
        echo "  Missing JSON: $AGGREGATED_JSON"
        echo "  Please ensure neoantigen prediction completed successfully."
        continue
    fi
    
    echo "Found required files:"
    echo "  TSV: $AGGREGATED_TSV"
    echo "  JSON: $AGGREGATED_JSON"
    
    # Count candidates
    CANDIDATE_COUNT=$(tail -n +2 "$AGGREGATED_TSV" | wc -l)
    echo "Total candidates to review: $CANDIDATE_COUNT"
    
    echo ""
    echo "======================================================"
    echo "PVACVIEW INSTRUCTIONS"
    echo "======================================================"
    echo "1. The pVACview server will start in your terminal"
    echo "2. Open your web browser and go to the URL shown"
    echo "3. On the 'PVACtools Output' tab:"
    echo "   - Browse for Aggregate Report: $AGGREGATED_TSV"
    echo "   - Browse for Metrics file: $AGGREGATED_JSON"
    echo "4. Click 'Visualize' to load your data"
    echo "5. Press CTRL+C to stop the server when finished"
    echo "======================================================"
    echo ""
    echo "Starting pVACview... (Press CTRL+C to stop)"
    
    pvacview run "$PVAC_OUTPUT"
done

conda deactivate

echo ""
echo "======================================================"
echo "Visualization complete"
echo "======================================================"