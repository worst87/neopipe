#!/bin/bash
# =================================================================================
# STEP 6: VARIANT PROCESSING AND ANNOTATION
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"
HELPER_DIR="$(dirname "$0")/helpers"

echo "======================================================"
echo "STEP 6: VARIANT PROCESSING AND ANNOTATION"
echo "======================================================"

# Part 1: Merge and normalize variants
echo ""
echo "[Part 1] Creating high-confidence variant callset"
echo "================================================"

# Initialize conda for merging
source ~/miniconda3/etc/profile.d/conda.sh
conda activate pvacseq_env

for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    NORMAL_SAMPLE="${TUMOR_NORMAL_PAIRS[$TUMOR_SAMPLE]}"
    
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing: $TUMOR_SAMPLE"
    
    # Define paths
    SAMPLE_VAR_DIR="$VAR_DIR/$TUMOR_SAMPLE"
    MERGED_DIR="$SAMPLE_VAR_DIR/merged"
    mkdir -p "$MERGED_DIR"
    
    MUTECT2_VCF="$SAMPLE_VAR_DIR/raw_calls/mutect2/${TUMOR_SAMPLE}.mutect2.pass.vcf.gz"
    VARSCAN_VCF="$SAMPLE_VAR_DIR/raw_calls/varscan/${TUMOR_SAMPLE}.varscan.pass.vcf.gz"
    
    # Step 1: Create high-confidence VCF
    echo "  [1/3] Merging variant calls..."
    HC_VCF="$MERGED_DIR/${TUMOR_SAMPLE}.hc.vcf"
    HC_SINGLE_VCF="$MERGED_DIR/${TUMOR_SAMPLE}.hc_single.vcf"
    
    python "$HELPER_DIR/make_hc_vcf.py" \
        --primary "$MUTECT2_VCF" \
        --primary_name "$PRIMARY_CALLER" \
        --confirming "$VARSCAN_VCF" \
        --confirming_names "VS" \
        --out_vcf "$HC_VCF" \
        --out_single_vcf "$HC_SINGLE_VCF" \
        --normal_sample_name "$NORMAL_SAMPLE"
    
    # Step 2: Compress and index
    echo "  [2/3] Compressing VCFs..."
    bgzip -f "$HC_VCF"
    bcftools index -f -t "${HC_VCF}.gz"
    
    if [[ -f "$HC_SINGLE_VCF" ]]; then
        bgzip -f "$HC_SINGLE_VCF"
        bcftools index -f -t "${HC_SINGLE_VCF}.gz"
    fi
    
    # Step 3: Normalize
    echo "  [3/3] Normalizing variants..."
    NORMALIZED_VCF="$MERGED_DIR/${TUMOR_SAMPLE}.hc.normalized.vcf.gz"
    
    bcftools norm -m- -f "$REF_GENOME_FA" -O z -o "$NORMALIZED_VCF" "${HC_VCF}.gz"
    bcftools index -f -t "$NORMALIZED_VCF"
    
    # Report statistics
    TOTAL_HC=$(bcftools view -H "${HC_VCF}.gz" | wc -l)
    if [[ -f "${HC_SINGLE_VCF}.gz" ]]; then
        TOTAL_SINGLE=$(bcftools view -H "${HC_SINGLE_VCF}.gz" | wc -l)
    else
        TOTAL_SINGLE=0
    fi
    
    echo "    High-confidence variants: $TOTAL_HC"
    echo "    Single-caller variants: $TOTAL_SINGLE"
done

# Part 2: VEP annotation
echo ""
echo "[Part 2] Annotating variants with VEP"
echo "====================================="

for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$TUMOR_SAMPLE" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing: $TUMOR_SAMPLE"
    
    # Define paths
    SAMPLE_VAR_DIR="$VAR_DIR/$TUMOR_SAMPLE"
    INPUT_VCF="$SAMPLE_VAR_DIR/merged/${TUMOR_SAMPLE}.hc.normalized.vcf.gz"
    ANNOTATED_DIR="$SAMPLE_VAR_DIR/annotated"
    mkdir -p "$ANNOTATED_DIR"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$ANNOTATED_DIR/.done" ]]; then
        echo "  Skipping annotation - already completed"
        continue
    fi
    
    # Step 1: VEP annotation
    echo "  [1/5] Running VEP..."
    VEP_VCF="$ANNOTATED_DIR/${TUMOR_SAMPLE}.vep.vcf"
    VEP_STATS="$ANNOTATED_DIR/${TUMOR_SAMPLE}.vep_summary.html"
    
    vep \
        --input_file "$INPUT_VCF" \
        --output_file "$VEP_VCF" \
        --stats_file "$VEP_STATS" \
        --force_overwrite \
        --format vcf --vcf --symbol --terms SO --tsl \
        --biotype --hgvs --fasta "$REF_GENOME_FA" \
        --offline --cache --dir_cache "$VEP_CACHE_DIR" \
        --cache_version "$VEP_CACHE_VERSION" \
        --species "$VEP_SPECIES" --assembly "$VEP_ASSEMBLY" \
        --fork "$THREADS" \
        --plugin Frameshift --plugin Wildtype \
        --dir_plugins "$VEP_PLUGINS_DIR" \
        --transcript_version
    
    # Step 2: Add expression data
    echo "  [2/5] Adding gene expression..."
    EXPR_VCF="$ANNOTATED_DIR/${TUMOR_SAMPLE}.vep.gx.vcf"
    TUMOR_TPM="$EXPR_DIR/gene_tpm/${TUMOR_SAMPLE}_gene_tpm.tsv"
    
    vcf-expression-annotator \
        "$VEP_VCF" \
        "$TUMOR_TPM" \
        custom gene \
        --id-column "gene_id" \
        --expression-column "tpm" \
        -s "$TUMOR_SAMPLE" \
        --ignore-ensembl-id-version \
        -o "$EXPR_VCF"
    
    # Step 3: Decompose for readcount annotation
    echo "  [3/5] Decomposing multi-allelic sites..."
    DECOMPOSED_VCF="$ANNOTATED_DIR/${TUMOR_SAMPLE}.decomposed.vcf"
    vt decompose -s -o "$DECOMPOSED_VCF" "$EXPR_VCF"
    
    # Step 4: Add readcount data
    echo "  [4/5] Adding RNA readcounts..."
    TEMP_DIR="$ANNOTATED_DIR/temp"
    mkdir -p "$TEMP_DIR"
    
    # Extract SNV and indel positions
    awk -F'\t' '!/^#/ { if (length($4) == 1 && length($5) == 1) print $1"\t"$2"\t"$2 }' \
        "$DECOMPOSED_VCF" > "$TEMP_DIR/snv.sites"
    awk -F'\t' '!/^#/ { if (length($4) != 1 || length($5) != 1) print $1"\t"$2"\t"$2 }' \
        "$DECOMPOSED_VCF" > "$TEMP_DIR/indel.sites"
    
    # Get BAM paths
    TUMOR_BAM="$ALIGN_DIR/processed/$TUMOR_SAMPLE/$TUMOR_SAMPLE.analysis_ready.bam"
    NORMAL_BAM="$ALIGN_DIR/processed/$NORMAL_SAMPLE/$NORMAL_SAMPLE.analysis_ready.bam"
    
    # Generate readcounts
    RC_VCF="$DECOMPOSED_VCF"
    STEP=1
    
    for SAMPLE_TYPE in tumor normal; do
        if [[ "$SAMPLE_TYPE" == "tumor" ]]; then
            BAM="$TUMOR_BAM"
            SAMPLE_NAME="$TUMOR_SAMPLE"
        else
            BAM="$NORMAL_BAM"
            SAMPLE_NAME="$NORMAL_SAMPLE"
        fi
        
        for VAR_TYPE in snv indel; do
            SITES_FILE="$TEMP_DIR/${VAR_TYPE}.sites"
            RC_FILE="$TEMP_DIR/${SAMPLE_TYPE}_readcount_${VAR_TYPE}.tsv"
            
            if [[ -s "$SITES_FILE" ]]; then
                echo "    Generating $SAMPLE_TYPE $VAR_TYPE readcounts..."
                if [[ "$VAR_TYPE" == "snv" ]]; then
                    bam-readcount -f "$REF_GENOME_FA" -l "$SITES_FILE" -w 1 -b 20 "$BAM" > "$RC_FILE"
                else
                    bam-readcount -f "$REF_GENOME_FA" -l "$SITES_FILE" -w 1 -i -b 20 "$BAM" > "$RC_FILE"
                fi
                
                OUT_VCF="$TEMP_DIR/step${STEP}.vcf"
                vcf-readcount-annotator "$RC_VCF" "$RC_FILE" RNA -s "$SAMPLE_NAME" -t "$VAR_TYPE" -o "$OUT_VCF"
                RC_VCF="$OUT_VCF"
                ((STEP++))
            fi
        done
    done
    
    # Step 5: Create final annotated VCF
    echo "  [5/5] Creating final annotated VCF..."
    FINAL_VCF="$ANNOTATED_DIR/${TUMOR_SAMPLE}.final.annotated.vcf.gz"
    
    bgzip -c "$RC_VCF" > "$FINAL_VCF"
    bcftools index -f -t "$FINAL_VCF"
    
    # Cleanup
    if [[ "$KEEP_INTERMEDIATE" != "true" ]]; then
        echo "  Cleaning up intermediate files..."
        rm -rf "$TEMP_DIR"
        rm -f "$VEP_VCF" "$EXPR_VCF" "$DECOMPOSED_VCF"
    fi
    
    # Mark as completed
    touch "$ANNOTATED_DIR/.done"
    echo "  Annotation complete"
done

conda deactivate

echo ""
echo "======================================================"
echo "Variant processing and annotation complete"
echo "Annotated VCFs: $VAR_DIR/*/annotated"
echo "======================================================"