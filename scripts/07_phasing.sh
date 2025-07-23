#!/bin/bash
# =================================================================================
# STEP 7: VARIANT PHASING
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Define Java 8 path for GATK3
JAVA8_PATH="$HOME/miniconda3/envs/analysis_env/bin/java"

if [[ ! -f "$JAVA8_PATH" ]]; then
    echo "ERROR: Java 8 not found at $JAVA8_PATH"
    exit 1
fi

echo "======================================================"
echo "STEP 7: VARIANT PHASING"
echo "======================================================"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

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
    PHASED_DIR="$SAMPLE_VAR_DIR/phased"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$PHASED_DIR/.done" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    # Clean and create directory
    rm -rf "$PHASED_DIR"
    mkdir -p "$PHASED_DIR"
    
    # Input files
    SOMATIC_VCF="$SAMPLE_VAR_DIR/merged/${TUMOR_SAMPLE}.hc.normalized.vcf.gz"
    GERMLINE_VCF="$SAMPLE_VAR_DIR/germline/${NORMAL_SAMPLE}.genotyped_germline.vcf.gz"
    TUMOR_BAM="$ALIGN_DIR/processed/$TUMOR_SAMPLE/$TUMOR_SAMPLE.analysis_ready.bam"
    
    # Step 1: Standardize headers
    echo "  [1/8] Standardizing VCF headers..."
    SOMATIC_STD="$PHASED_DIR/somatic.std_header.vcf.gz"
    GERMLINE_STD="$PHASED_DIR/germline.std_header.vcf.gz"
    
    bcftools reheader -f "${REF_GENOME_FA}.fai" "$SOMATIC_VCF" -o "$SOMATIC_STD"
    bcftools index -f -t "$SOMATIC_STD"
    
    bcftools reheader -f "${REF_GENOME_FA}.fai" "$GERMLINE_VCF" -o "$GERMLINE_STD"
    bcftools index -f -t "$GERMLINE_STD"
    
    # Step 2: Extract tumor-only somatic
    echo "  [2/8] Extracting tumor-only variants..."
    TUMOR_ONLY_VCF="$PHASED_DIR/tumor_only.vcf.gz"
    gatk SelectVariants \
        -R "$REF_GENOME_FA" \
        -V "$SOMATIC_STD" \
        --sample-name "$TUMOR_SAMPLE" \
        -O "$TUMOR_ONLY_VCF"
    
    # Step 3: Rename germline sample
    echo "  [3/8] Renaming germline sample..."
    GERMLINE_RENAMED="$PHASED_DIR/germline_renamed.vcf.gz"
    gatk RenameSampleInVcf \
        -I "$GERMLINE_STD" \
        -O "$GERMLINE_RENAMED" \
        --NEW_SAMPLE_NAME "$TUMOR_SAMPLE"
    
    # Step 4: Combine variants
    echo "  [4/8] Combining somatic and germline variants..."
    COMBINED_VCF="$PHASED_DIR/combined_for_phasing.vcf.gz"
    gatk MergeVcfs \
        -I "$TUMOR_ONLY_VCF" \
        -I "$GERMLINE_RENAMED" \
        -O "$COMBINED_VCF"
    
    # Step 5: Sort for GATK3
    echo "  [5/8] Sorting combined VCF..."
    SORTED_VCF="$PHASED_DIR/combined_sorted.vcf"
    bcftools sort -O v -o "$SORTED_VCF" "$COMBINED_VCF"
    
    # Step 6: Run ReadBackedPhasing
    echo "  [6/8] Running ReadBackedPhasing..."
    PHASED_RAW="$PHASED_DIR/all_phased.vcf"
    
    "$JAVA8_PATH" -Xmx"$JAVA_MEM" -jar "$GATK3_JAR" \
        -T ReadBackedPhasing \
        -R "$REF_GENOME_FA" \
        -I "$TUMOR_BAM" \
        --variant "$SORTED_VCF" \
        -o "$PHASED_RAW" \
        --phaseQualityThresh 20.0
    
    # Step 7: Re-annotate with VEP
    echo "  [7/8] Re-annotating phased VCF..."
    conda deactivate
    conda activate pvacseq_env
    
    FINAL_PHASED="$PHASED_DIR/${TUMOR_SAMPLE}.phased.annotated.vcf.gz"
    VEP_STATS="$PHASED_DIR/${TUMOR_SAMPLE}.phased.vep_summary.html"
    
    vep \
        --input_file "$PHASED_RAW" \
        --output_file stdout \
        --stats_file "$VEP_STATS" \
        --format vcf --vcf --symbol --terms SO --tsl \
        --biotype --hgvs --fasta "$REF_GENOME_FA" \
        --offline --cache --dir_cache "$VEP_CACHE_DIR" \
        --cache_version "$VEP_CACHE_VERSION" \
        --species "$VEP_SPECIES" --assembly "$VEP_ASSEMBLY" \
        --fork "$THREADS" \
        --force_overwrite \
        --plugin Frameshift --plugin Wildtype \
        --dir_plugins "$VEP_PLUGINS_DIR" \
        --transcript_version | \
    bgzip -c > "$FINAL_PHASED"
    
    bcftools index -f -t "$FINAL_PHASED"
    
    conda deactivate
    conda activate variant_env
    
    # Step 8: Report statistics
    echo "  [8/8] Generating statistics..."
    TOTAL_SOMATIC=$(bcftools view -H "$SOMATIC_STD" | wc -l)
    FINAL_COUNT=$(bcftools view -H "$FINAL_PHASED" 2>/dev/null | wc -l || echo 0)
    PHASED_COUNT=$(bcftools view -H -i 'FORMAT/PS!="."' "$FINAL_PHASED" 2>/dev/null | wc -l || echo 0)
    
    echo "    Original somatic variants: $TOTAL_SOMATIC"
    echo "    Total variants in phased VCF: $FINAL_COUNT"
    echo "    Successfully phased variants: $PHASED_COUNT"
    
    # Cleanup
    if [[ "$KEEP_INTERMEDIATE" != "true" ]]; then
        echo "  Cleaning up intermediate files..."
        rm -f "$SOMATIC_STD" "$GERMLINE_STD" "$TUMOR_ONLY_VCF" "$GERMLINE_RENAMED"
        rm -f "$COMBINED_VCF" "$SORTED_VCF" "$PHASED_RAW"
        rm -f "$SOMATIC_STD.tbi" "$GERMLINE_STD.tbi" "$TUMOR_ONLY_VCF.tbi"
        rm -f "$GERMLINE_RENAMED.tbi" "$COMBINED_VCF.tbi"
    fi
    
    # Mark as completed
    touch "$PHASED_DIR/.done"
    echo "  Phasing complete"
done

conda deactivate

echo ""
echo "======================================================"
echo "Phasing complete"
echo "Phased VCFs: $VAR_DIR/*/phased"
echo "======================================================"