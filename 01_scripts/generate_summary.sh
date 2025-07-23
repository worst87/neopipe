#!/bin/bash
# =================================================================================
# SUMMARY REPORT GENERATOR
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Output file
SUMMARY_REPORT="$OUTPUT_DIR/summary_report.html"

echo "Generating summary report..."

# Start HTML
cat > "$SUMMARY_REPORT" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>RNA-seq Neoantigen Pipeline Summary</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1, h2 { color: #2e4057; }
        table { border-collapse: collapse; width: 100%; margin: 20px 0; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .success { color: green; }
        .warning { color: orange; }
        .error { color: red; }
        .section { margin: 30px 0; }
        .metric { display: inline-block; margin: 10px 20px 10px 0; }
        .metric-value { font-size: 24px; font-weight: bold; color: #2e4057; }
        .metric-label { font-size: 14px; color: #666; }
    </style>
</head>
<body>
EOF

# Add header
echo "<h1>RNA-seq Neoantigen Pipeline Summary</h1>" >> "$SUMMARY_REPORT"
echo "<p>Generated: $(date)</p>" >> "$SUMMARY_REPORT"
echo "<p>Pipeline Version: $PIPELINE_VERSION</p>" >> "$SUMMARY_REPORT"

# Overall statistics
echo "<div class='section'>" >> "$SUMMARY_REPORT"
echo "<h2>Overall Statistics</h2>" >> "$SUMMARY_REPORT"

TOTAL_SAMPLES=${#ALL_SAMPLES[@]}
TUMOR_COUNT=${#TUMOR_SAMPLES[@]}
NORMAL_COUNT=${#NORMAL_SAMPLES[@]}

echo "<div class='metric'><div class='metric-value'>$TOTAL_SAMPLES</div><div class='metric-label'>Total Samples</div></div>" >> "$SUMMARY_REPORT"
echo "<div class='metric'><div class='metric-value'>$TUMOR_COUNT</div><div class='metric-label'>Tumor Samples</div></div>" >> "$SUMMARY_REPORT"
echo "<div class='metric'><div class='metric-value'>$NORMAL_COUNT</div><div class='metric-label'>Normal Samples</div></div>" >> "$SUMMARY_REPORT"
echo "</div>" >> "$SUMMARY_REPORT"

# Sample processing status
echo "<div class='section'>" >> "$SUMMARY_REPORT"
echo "<h2>Sample Processing Status</h2>" >> "$SUMMARY_REPORT"
echo "<table>" >> "$SUMMARY_REPORT"
echo "<tr><th>Sample</th><th>Type</th><th>QC</th><th>Alignment</th><th>Expression</th><th>Status</th></tr>" >> "$SUMMARY_REPORT"

for sample in "${ALL_SAMPLES[@]}"; do
    echo "<tr>" >> "$SUMMARY_REPORT"
    echo "<td>$sample</td>" >> "$SUMMARY_REPORT"
    
    # Sample type
    if [[ " ${TUMOR_SAMPLES[@]} " =~ " $sample " ]]; then
        echo "<td>Tumor</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>Normal</td>" >> "$SUMMARY_REPORT"
    fi
    
    # QC status
    if [[ -f "$QC_DIR/fastqc/$sample/.done" ]]; then
        echo "<td class='success'>✓</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td class='error'>✗</td>" >> "$SUMMARY_REPORT"
    fi
    
    # Alignment status
    if [[ -f "$ALIGN_DIR/processed/$sample/.done" ]]; then
        echo "<td class='success'>✓</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td class='error'>✗</td>" >> "$SUMMARY_REPORT"
    fi
    
    # Expression status
    if [[ -f "$EXPR_DIR/gene_tpm/${sample}_gene_tpm.tsv" ]]; then
        echo "<td class='success'>✓</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td class='error'>✗</td>" >> "$SUMMARY_REPORT"
    fi
    
    # Overall status
    if [[ -f "$ALIGN_DIR/processed/$sample/.done" ]]; then
        echo "<td class='success'>Complete</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td class='warning'>In Progress</td>" >> "$SUMMARY_REPORT"
    fi
    
    echo "</tr>" >> "$SUMMARY_REPORT"
done
echo "</table>" >> "$SUMMARY_REPORT"
echo "</div>" >> "$SUMMARY_REPORT"

# Variant calling results
echo "<div class='section'>" >> "$SUMMARY_REPORT"
echo "<h2>Variant Calling Results</h2>" >> "$SUMMARY_REPORT"
echo "<table>" >> "$SUMMARY_REPORT"
echo "<tr><th>Tumor Sample</th><th>Normal Sample</th><th>Mutect2</th><th>VarScan2</th><th>High-Confidence</th><th>Phased</th></tr>" >> "$SUMMARY_REPORT"

for tumor in "${TUMOR_SAMPLES[@]}"; do
    normal="${TUMOR_NORMAL_PAIRS[$tumor]}"
    echo "<tr>" >> "$SUMMARY_REPORT"
    echo "<td>$tumor</td>" >> "$SUMMARY_REPORT"
    echo "<td>$normal</td>" >> "$SUMMARY_REPORT"
    
    # Mutect2 count
    MUTECT2_VCF="$VAR_DIR/$tumor/raw_calls/mutect2/${tumor}.mutect2.pass.vcf.gz"
    if [[ -f "$MUTECT2_VCF" ]]; then
        MUTECT2_COUNT=$(bcftools view -H "$MUTECT2_VCF" | wc -l)
        echo "<td>$MUTECT2_COUNT</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
    fi
    
    # VarScan2 count
    VARSCAN_VCF="$VAR_DIR/$tumor/raw_calls/varscan/${tumor}.varscan.pass.vcf.gz"
    if [[ -f "$VARSCAN_VCF" ]]; then
        VARSCAN_COUNT=$(bcftools view -H "$VARSCAN_VCF" | wc -l)
        echo "<td>$VARSCAN_COUNT</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
    fi
    
    # High-confidence count
    HC_VCF="$VAR_DIR/$tumor/merged/${tumor}.hc.vcf.gz"
    if [[ -f "$HC_VCF" ]]; then
        HC_COUNT=$(bcftools view -H "$HC_VCF" | wc -l)
        echo "<td>$HC_COUNT</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
    fi
    
    # Phased count
    PHASED_VCF="$VAR_DIR/$tumor/phased/${tumor}.phased.annotated.vcf.gz"
    if [[ -f "$PHASED_VCF" ]]; then
        PHASED_COUNT=$(bcftools view -H -i 'FORMAT/PS!="."' "$PHASED_VCF" 2>/dev/null | wc -l || echo 0)
        echo "<td>$PHASED_COUNT</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
    fi
    
    echo "</tr>" >> "$SUMMARY_REPORT"
done
echo "</table>" >> "$SUMMARY_REPORT"
echo "</div>" >> "$SUMMARY_REPORT"

# Neoantigen results
echo "<div class='section'>" >> "$SUMMARY_REPORT"
echo "<h2>Neoantigen Prediction Results</h2>" >> "$SUMMARY_REPORT"
echo "<table>" >> "$SUMMARY_REPORT"
echo "<tr><th>Tumor Sample</th><th>Total Epitopes</th><th>Filtered Candidates</th><th>Report Location</th></tr>" >> "$SUMMARY_REPORT"

for tumor in "${TUMOR_SAMPLES[@]}"; do
    echo "<tr>" >> "$SUMMARY_REPORT"
    echo "<td>$tumor</td>" >> "$SUMMARY_REPORT"
    
    RESULTS_TSV="$NEO_DIR/$tumor/pvacseq/MHC_Class_I/${tumor}.all_epitopes.tsv"
    AGGREGATED_TSV="$NEO_DIR/$tumor/reports/${tumor}.all_epitopes.aggregated.tsv"
    
    # Total epitopes
    if [[ -f "$RESULTS_TSV" ]]; then
        TOTAL_EPITOPES=$(tail -n +2 "$RESULTS_TSV" | wc -l)
        echo "<td>$TOTAL_EPITOPES</td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
    fi
    
    # Filtered candidates
    if [[ -f "$AGGREGATED_TSV" ]]; then
        FILTERED_COUNT=$(tail -n +2 "$AGGREGATED_TSV" | wc -l)
        echo "<td>$FILTERED_COUNT</td>" >> "$SUMMARY_REPORT"
        echo "<td><a href='file://$AGGREGATED_TSV'>View Report</a></td>" >> "$SUMMARY_REPORT"
    else
        echo "<td>-</td>" >> "$SUMMARY_REPORT"
        echo "<td>Not available</td>" >> "$SUMMARY_REPORT"
    fi
    
    echo "</tr>" >> "$SUMMARY_REPORT"
done
echo "</table>" >> "$SUMMARY_REPORT"
echo "</div>" >> "$SUMMARY_REPORT"

# Alignment statistics
echo "<div class='section'>" >> "$SUMMARY_REPORT"
echo "<h2>Alignment Statistics</h2>" >> "$SUMMARY_REPORT"
echo "<table>" >> "$SUMMARY_REPORT"
echo "<tr><th>Sample</th><th>Total Reads</th><th>Uniquely Mapped</th><th>% Uniquely Mapped</th></tr>" >> "$SUMMARY_REPORT"

for sample in "${ALL_SAMPLES[@]}"; do
    STAR_LOG="$ALIGN_DIR/star/$sample/${sample}.Log.final.out"
    if [[ -f "$STAR_LOG" ]]; then
        TOTAL_READS=$(grep "Number of input reads" "$STAR_LOG" | awk '{print $NF}')
        UNIQUE_MAPPED=$(grep "Uniquely mapped reads number" "$STAR_LOG" | awk '{print $NF}')
        UNIQUE_PERCENT=$(grep "Uniquely mapped reads %" "$STAR_LOG" | awk '{print $NF}')
        
        echo "<tr>" >> "$SUMMARY_REPORT"
        echo "<td>$sample</td>" >> "$SUMMARY_REPORT"
        echo "<td>$TOTAL_READS</td>" >> "$SUMMARY_REPORT"
        echo "<td>$UNIQUE_MAPPED</td>" >> "$SUMMARY_REPORT"
        echo "<td>$UNIQUE_PERCENT</td>" >> "$SUMMARY_REPORT"
        echo "</tr>" >> "$SUMMARY_REPORT"
    fi
done
echo "</table>" >> "$SUMMARY_REPORT"
echo "</div>" >> "$SUMMARY_REPORT"

# Close HTML
cat >> "$SUMMARY_REPORT" << 'EOF'
</body>
</html>
EOF

echo "Summary report generated: $SUMMARY_REPORT"