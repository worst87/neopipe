#!/bin/bash
# =================================================================================
# STEP 5: GERMLINE VARIANT CALLING FOR PHASING
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

echo "======================================================"
echo "STEP 5: GERMLINE CALLING FOR PHASING"
echo "======================================================"

# Process tumor-normal pairs
for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    NORMAL_SAMPLE="${TUMOR_NORMAL_PAIRS[$TUMOR_SAMPLE]}"
    
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing normal sample: $NORMAL_SAMPLE (for tumor: $TUMOR_SAMPLE)"
    
    # Define paths
    NORMAL_BAM="$ALIGN_DIR/processed/$NORMAL_SAMPLE/$NORMAL_SAMPLE.analysis_ready.bam"
    GERMLINE_DIR="$VAR_DIR/$TUMOR_SAMPLE/germline"
    mkdir -p "$GERMLINE_DIR"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$GERMLINE_DIR/.done" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    # Step 1: HaplotypeCaller in GVCF mode
    echo "  [1/3] Running HaplotypeCaller (RNA-optimized)..."
    GVCF_FILE="$GERMLINE_DIR/${NORMAL_SAMPLE}.raw_germline.g.vcf.gz"
    
    gatk HaplotypeCaller \
        -R "$REF_GENOME_FA" \
        -I "$NORMAL_BAM" \
        -O "$GVCF_FILE" \
        --dbsnp "$GERMLINE_COMMON_SNPS_VCF" \
        -ERC GVCF \
        --dont-use-soft-clipped-bases \
        --standard-min-confidence-threshold-for-calling 20.0 \
        --min-base-quality-score 20
    
    # Step 2: Genotype GVCF
    echo "  [2/3] Genotyping GVCF..."
    GENOTYPED_VCF="$GERMLINE_DIR/${NORMAL_SAMPLE}.genotyped_germline.vcf.gz"
    
    gatk GenotypeGVCFs \
        -R "$REF_GENOME_FA" \
        -V "$GVCF_FILE" \
        -O "$GENOTYPED_VCF" \
        --standard-min-confidence-threshold-for-calling 20.0
    
    # Step 3: Filter variants
    echo "  [3/3] Filtering germline variants (RNA-specific)..."
    FILTERED_VCF="$GERMLINE_DIR/${NORMAL_SAMPLE}.filtered_germline.vcf.gz"
    
    gatk VariantFiltration \
        -R "$REF_GENOME_FA" \
        -V "$GENOTYPED_VCF" \
        -O "$FILTERED_VCF" \
        --cluster-window-size 35 --cluster-size 3 \
        --filter-expression "QD < 2.0" --filter-name "QD_Filter" \
        --filter-expression "FS > 30.0" --filter-name "FS_Filter" \
        --filter-expression "SOR > 4.0" --filter-name "SOR_Filter"
    
    # Replace genotyped with filtered version
    mv "$FILTERED_VCF" "$GENOTYPED_VCF"
    mv "$FILTERED_VCF.tbi" "$GENOTYPED_VCF.tbi"
    
    # Report variant count
    GERMLINE_COUNT=$(bcftools view -H "$GENOTYPED_VCF" | wc -l)
    echo "  Germline variants for phasing: $GERMLINE_COUNT"
    
    # Mark as completed
    touch "$GERMLINE_DIR/.done"
done

conda deactivate

echo ""
echo "======================================================"
echo "Germline calling complete"
echo "Germline VCFs: $VAR_DIR/*/germline"
echo "======================================================"