#!/bin/bash
# =================================================================================
# STEP 4: SOMATIC VARIANT CALLING
# Version 1.0
# =================================================================================

set -eo pipefail

# Set locale to prevent VarScan warnings
export LC_ALL=C

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

echo "======================================================"
echo "STEP 4: SOMATIC VARIANT CALLING"
echo "======================================================"

# Process tumor-normal pairs
for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    NORMAL_SAMPLE="${TUMOR_NORMAL_PAIRS[$TUMOR_SAMPLE]}"
    
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing: $TUMOR_SAMPLE (tumor) vs $NORMAL_SAMPLE (normal)"
    echo "================================================"
    
    # Define paths
    TUMOR_BAM="$ALIGN_DIR/processed/$TUMOR_SAMPLE/$TUMOR_SAMPLE.analysis_ready.bam"
    NORMAL_BAM="$ALIGN_DIR/processed/$NORMAL_SAMPLE/$NORMAL_SAMPLE.analysis_ready.bam"
    
    SAMPLE_VAR_DIR="$VAR_DIR/$TUMOR_SAMPLE"
    RAW_CALLS_DIR="$SAMPLE_VAR_DIR/raw_calls"
    mkdir -p "$RAW_CALLS_DIR/mutect2" "$RAW_CALLS_DIR/varscan"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$RAW_CALLS_DIR/.done" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    # Part 1: Mutect2
    echo "  [Part 1] Running Mutect2..."
    MUTECT2_RAW="$RAW_CALLS_DIR/mutect2/${TUMOR_SAMPLE}.raw.vcf.gz"
    MUTECT2_FILTERED="$RAW_CALLS_DIR/mutect2/${TUMOR_SAMPLE}.filtered.vcf.gz"
    MUTECT2_PASS="$RAW_CALLS_DIR/mutect2/${TUMOR_SAMPLE}.mutect2.pass.vcf.gz"
    
    echo "    [1/3] Calling variants..."
    gatk Mutect2 \
        -R "$REF_GENOME_FA" \
        -I "$TUMOR_BAM" \
        -I "$NORMAL_BAM" \
        -tumor "$TUMOR_SAMPLE" \
        -normal "$NORMAL_SAMPLE" \
        --native-pair-hmm-threads "$THREADS" \
        --germline-resource "$GERMLINE_COMMON_SNPS_VCF" \
        --af-of-alleles-not-in-resource "$AF_OF_ALLELES_NOT_IN_RESOURCE" \
        -O "$MUTECT2_RAW"
    
    echo "    [2/3] Filtering calls..."
    gatk FilterMutectCalls \
        -R "$REF_GENOME_FA" \
        -V "$MUTECT2_RAW" \
        -O "$MUTECT2_FILTERED"
    
    echo "    [3/3] Extracting PASS variants..."
    bcftools view -f PASS "$MUTECT2_FILTERED" | \
    bcftools filter -e "FORMAT/AD[0:1] < $MIN_AD" -O z -o "$MUTECT2_PASS"
    bcftools index -f -t "$MUTECT2_PASS"
    
    # Part 2: VarScan2
    echo "  [Part 2] Running VarScan2..."
    VARSCAN_DIR="$RAW_CALLS_DIR/varscan"
    
    echo "    [1/5] Generating mpileup..."
    samtools mpileup -B -f "$REF_GENOME_FA" "$NORMAL_BAM" "$TUMOR_BAM" > "$VARSCAN_DIR/somatic.mpileup"
    
    echo "    [2/5] Calling somatic variants..."
    varscan somatic \
        "$VARSCAN_DIR/somatic.mpileup" \
        "$VARSCAN_DIR/somatic" \
        --mpileup 1 \
        --output-vcf 1 \
        --min-coverage-normal "$VARSCAN_MIN_COV_NORMAL" \
        --min-coverage-tumor "$VARSCAN_MIN_COV_TUMOR" \
        --p-value "$VARSCAN_SOMATIC_PVALUE" \
        --somatic-p-value "$VARSCAN_SOMATIC_SOMATICPVALUE" \
        --strand-filter "$VARSCAN_STRAND_FILTER"
    
    echo "    [3/5] Processing somatic calls..."
    varscan processSomatic "$VARSCAN_DIR/somatic.snp.vcf"
    varscan processSomatic "$VARSCAN_DIR/somatic.indel.vcf"
    
    echo "    [4/5] Formatting output..."
    # Rename samples
    echo -e "TUMOR\t${TUMOR_SAMPLE}\nNORMAL\t${NORMAL_SAMPLE}" > "$VARSCAN_DIR/sample_rename.txt"
    
    # Process SNVs
    bcftools reheader -s "$VARSCAN_DIR/sample_rename.txt" "$VARSCAN_DIR/somatic.snp.Somatic.hc.vcf" | \
    bgzip -c > "$VARSCAN_DIR/somatic.snp.Somatic.hc.vcf.gz"
    bcftools index -f -t "$VARSCAN_DIR/somatic.snp.Somatic.hc.vcf.gz"
    
    # Process Indels
    bcftools reheader -s "$VARSCAN_DIR/sample_rename.txt" "$VARSCAN_DIR/somatic.indel.Somatic.hc.vcf" | \
    bgzip -c > "$VARSCAN_DIR/somatic.indel.Somatic.hc.vcf.gz"
    bcftools index -f -t "$VARSCAN_DIR/somatic.indel.Somatic.hc.vcf.gz"
    
    # Merge SNVs and Indels
    VARSCAN_MERGED="$VARSCAN_DIR/${TUMOR_SAMPLE}.varscan.merged.vcf.gz"
    bcftools concat -a -O z -o "$VARSCAN_MERGED" \
        "$VARSCAN_DIR/somatic.snp.Somatic.hc.vcf.gz" \
        "$VARSCAN_DIR/somatic.indel.Somatic.hc.vcf.gz"
    
    echo "    [5/5] Standardizing header..."
    VARSCAN_PASS="$RAW_CALLS_DIR/varscan/${TUMOR_SAMPLE}.varscan.pass.vcf.gz"
    bcftools reheader -f "${REF_GENOME_FA}.fai" "$VARSCAN_MERGED" -o "$VARSCAN_PASS"
    bcftools index -f -t "$VARSCAN_PASS"
    
    # Cleanup
    if [[ "$KEEP_INTERMEDIATE" != "true" ]]; then
        echo "  Cleaning up intermediate files..."
        rm -f "$VARSCAN_DIR/somatic.mpileup" "$VARSCAN_DIR/sample_rename.txt"
        rm -f "$VARSCAN_DIR"/somatic.snp.* "$VARSCAN_DIR"/somatic.indel.*
        rm -f "$VARSCAN_MERGED"
    fi
    
    # Report variant counts
    MUTECT2_COUNT=$(bcftools view -H "$MUTECT2_PASS" | wc -l)
    VARSCAN_COUNT=$(bcftools view -H "$VARSCAN_PASS" | wc -l)
    echo ""
    echo "  Variant calling complete:"
    echo "    Mutect2: $MUTECT2_COUNT variants"
    echo "    VarScan2: $VARSCAN_COUNT variants"
    
    # Mark as completed
    touch "$RAW_CALLS_DIR/.done"
done

conda deactivate

echo ""
echo "======================================================"
echo "Variant calling complete"
echo "Raw calls: $VAR_DIR/*/raw_calls"
echo "======================================================"