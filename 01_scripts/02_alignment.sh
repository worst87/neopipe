#!/bin/bash
# =================================================================================
# STEP 2: ALIGNMENT AND GATK PREPROCESSING
# Version 1.0
# =================================================================================

set -eo pipefail

# Source configuration
source "$(dirname "$0")/config.sh"

# Initialize conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

echo "======================================================"
echo "STEP 2: ALIGNMENT AND PREPROCESSING"
echo "======================================================"

# Setup reference indices
echo "Checking reference indices..."
REF_PARENT_DIR=$(dirname "$REF_GENOME_FA")
REF_DICT="${REF_GENOME_FA%.fa}.dict"
STAR_INDEX_DIR="${REF_PARENT_DIR}/STAR_genome_index"

# Create reference indices if needed
[[ ! -f "$REF_DICT" ]] && echo "Creating sequence dictionary..." && gatk CreateSequenceDictionary -R "$REF_GENOME_FA" -O "$REF_DICT"
[[ ! -f "${REF_GENOME_FA}.fai" ]] && echo "Creating FASTA index..." && samtools faidx "$REF_GENOME_FA"
[[ ! -f "${GERMLINE_SNPS_VCF}.tbi" ]] && echo "Indexing SNPs VCF..." && tabix -p vcf "$GERMLINE_SNPS_VCF"
[[ ! -f "${GERMLINE_INDELS_VCF}.tbi" ]] && echo "Indexing indels VCF..." && tabix -p vcf "$GERMLINE_INDELS_VCF"

# Build STAR index
if [[ ! -d "$STAR_INDEX_DIR" ]]; then
    echo "Building STAR genome index..."
    mkdir -p "$STAR_INDEX_DIR"
    STAR --runMode genomeGenerate \
         --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$REF_GENOME_FA" \
         --sjdbGTFfile "$REF_GENOME_GTF" \
         --sjdbOverhang "$STAR_OVERHANG"
fi

# Process each sample
for sample_id in "${ALL_SAMPLES[@]}"; do
    # Skip if specific sample requested
    if [[ -n "$SAMPLE" && "$SAMPLE" != "$sample_id" ]]; then
        continue
    fi
    
    echo ""
    echo "Processing sample: $sample_id"
    echo "================================"
    
    # Define paths
    R1_PATH="${READ_PATHS[$sample_id]}/${READ1_FILENAMES[$sample_id]}"
    R2_PATH="${READ_PATHS[$sample_id]}/${READ2_FILENAMES[$sample_id]}"
    STAR_OUTPUT_DIR="$ALIGN_DIR/star/$sample_id"
    PROCESSED_DIR="$ALIGN_DIR/processed/$sample_id"
    SAMPLE_QC_DIR="$QC_DIR/qualimap/$sample_id"
    
    mkdir -p "$STAR_OUTPUT_DIR" "$PROCESSED_DIR" "$SAMPLE_QC_DIR" "$QC_DIR/bqsr_plots"
    
    # Output files
    ALIGNED_BAM="$STAR_OUTPUT_DIR/${sample_id}.Aligned.sortedByCoord.out.bam"
    FINAL_BAM="$PROCESSED_DIR/${sample_id}.analysis_ready.bam"
    
    # Check if already completed
    if [[ "$SKIP_COMPLETED" == "true" && -f "$PROCESSED_DIR/.done" ]]; then
        echo "  Skipping - already completed"
        continue
    fi
    
    # Step 1: STAR Alignment
    echo "  [1/5] Running STAR alignment..."
    STAR --runThreadN "$THREADS" \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$R1_PATH" "$R2_PATH" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$STAR_OUTPUT_DIR/${sample_id}." \
         --outSAMtype BAM SortedByCoordinate \
         --twopassMode Basic \
         --outSAMunmapped Within \
         --outSAMattributes NH HI AS nM NM MD ch \
         --outFilterType BySJout \
         --outFilterMultimapNmax "$STAR_OUTFILTERMULTIMAPNMAX" \
         --outFilterMismatchNmax "$STAR_OUTFILTERMISMATCHNMAX" \
         --outFilterMismatchNoverLmax "$STAR_OUTFILTERMISMATCHNOVERLMAX" \
         --alignIntronMin "$STAR_ALIGNINTRONMIN" \
         --alignIntronMax "$STAR_ALIGNINTRONMAX" \
         --alignMatesGapMax "$STAR_ALIGNMATESGAPMAX" \
         --alignSJoverhangMin "$STAR_ALIGNSJOVERHANGMIN" \
         --alignSJDBoverhangMin "$STAR_ALIGNSJDBOVERHANGMIN"
    
    # Step 2: Add Read Groups
    echo "  [2/5] Adding read groups..."
    RG_BAM="$PROCESSED_DIR/${sample_id}.rg.bam"
    picard AddOrReplaceReadGroups \
        I="$ALIGNED_BAM" \
        O="$RG_BAM" \
        RGID="$sample_id" \
        RGLB="lib1" \
        RGPL="ILLUMINA" \
        RGPU="unit1" \
        RGSM="$sample_id"
    
    # Step 3: Mark Duplicates
    echo "  [3/5] Marking duplicates..."
    DEDUP_BAM="$PROCESSED_DIR/${sample_id}.dedup.bam"
    picard MarkDuplicates \
        I="$RG_BAM" \
        O="$DEDUP_BAM" \
        M="$PROCESSED_DIR/dedup_metrics.txt" \
        CREATE_INDEX=true
    
    # Step 4: GATK Processing
    echo "  [4/5] GATK preprocessing..."
    
    # Set NM, MD, UQ tags
    DEDUP_TAGS_BAM="$PROCESSED_DIR/${sample_id}.dedup_tags.bam"
    gatk SetNmMdAndUqTags \
        -R "$REF_GENOME_FA" \
        -I "$DEDUP_BAM" \
        -O "$DEDUP_TAGS_BAM" \
        --CREATE_INDEX true
    
    # Split N Cigar reads
    SPLIT_BAM="$PROCESSED_DIR/${sample_id}.split.bam"
    gatk SplitNCigarReads \
        -R "$REF_GENOME_FA" \
        -I "$DEDUP_TAGS_BAM" \
        -O "$SPLIT_BAM"
    
    # Base quality score recalibration
    RECAL_TABLE="$PROCESSED_DIR/recal.table"
    gatk BaseRecalibrator \
        -R "$REF_GENOME_FA" \
        -I "$SPLIT_BAM" \
        --known-sites "$GERMLINE_SNPS_VCF" \
        --known-sites "$GERMLINE_INDELS_VCF" \
        -O "$RECAL_TABLE"
    
    gatk ApplyBQSR \
        -R "$REF_GENOME_FA" \
        -I "$SPLIT_BAM" \
        --bqsr-recal-file "$RECAL_TABLE" \
        -O "$FINAL_BAM"
    
    # Generate BQSR plots
    POST_RECAL_TABLE="$PROCESSED_DIR/post_recal.table"
    gatk BaseRecalibrator \
        -R "$REF_GENOME_FA" \
        -I "$FINAL_BAM" \
        --known-sites "$GERMLINE_SNPS_VCF" \
        --known-sites "$GERMLINE_INDELS_VCF" \
        -O "$POST_RECAL_TABLE"
    
    gatk AnalyzeCovariates \
        -before "$RECAL_TABLE" \
        -after "$POST_RECAL_TABLE" \
        -plots "$QC_DIR/bqsr_plots/${sample_id}.bqsr_plots.pdf"
    
    # Step 5: Post-alignment QC
    echo "  [5/5] Running Qualimap..."
    qualimap rnaseq \
        -bam "$FINAL_BAM" \
        -gtf "$REF_GENOME_GTF" \
        -outdir "$SAMPLE_QC_DIR" \
        --java-mem-size="$JAVA_MEM"
    
    # Cleanup intermediate files if requested
    if [[ "$KEEP_INTERMEDIATE" != "true" ]]; then
        echo "  Cleaning up intermediate files..."
        rm -f "$RG_BAM" "$DEDUP_BAM" "$DEDUP_TAGS_BAM" "$SPLIT_BAM"
        rm -f "$RG_BAM.bai" "$DEDUP_BAM.bai" "$DEDUP_TAGS_BAM.bai" "$SPLIT_BAM.bai"
    fi
    
    # Mark as completed
    touch "$PROCESSED_DIR/.done"
    echo "  Completed processing for $sample_id"
done

conda deactivate

echo ""
echo "======================================================"
echo "Alignment and preprocessing complete"
echo "Processed BAMs: $ALIGN_DIR/processed"
echo "QC reports: $QC_DIR"
echo "======================================================"