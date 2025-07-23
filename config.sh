#!/bin/bash
# =================================================================================
# RNA-SEQ VARIANT & NEOANTIGEN ANALYSIS PIPELINE CONFIGURATION
# Version 1.0
# =================================================================================

# =================================================================================
# REFERENCE FILES
# =================================================================================
UNIVERSAL_DATA_DIR="/home/worst/projects/universal_data"
REF_GENOME_FA="${UNIVERSAL_DATA_DIR}/reference/Mus_musculus.GRCm39.105/Mus_musculus.GRCm39.dna.primary_assembly.fa"
REF_GENOME_GTF="${UNIVERSAL_DATA_DIR}/reference/Mus_musculus.GRCm39.105/Mus_musculus.GRCm39.105.gtf"
REF_TRANSCRIPTOME_FA="${UNIVERSAL_DATA_DIR}/reference/Mus_musculus.GRCm39.105/Mus_musculus.GRCm39.cdna.all.fa.gz"
GERMLINE_SNPS_VCF="${UNIVERSAL_DATA_DIR}/known_sites/mousegenomeproject/mgp_REL2021_snps.rsID.vcf.gz"
GERMLINE_INDELS_VCF="${UNIVERSAL_DATA_DIR}/known_sites/mousegenomeproject/mgp_REL2021_indels.rsID.vcf.gz"
GERMLINE_COMMON_SNPS_VCF="${UNIVERSAL_DATA_DIR}/known_sites/mousegenomeproject/mgp_REL2021_common_snps.vcf.gz"

# =================================================================================
# PROJECT CONFIGURATION
# =================================================================================
PROJ_DIR="/home/worst/projects/rnaseq_project"
PIPELINE_VERSION="1.0"
RUN_DATE=$(date +%Y%m%d_%H%M%S)

# =================================================================================
# OUTPUT DIRECTORY STRUCTURE
# =================================================================================
OUTPUT_DIR="$PROJ_DIR/pipeline_output"
QC_DIR="$OUTPUT_DIR/QC"
ALIGN_DIR="$OUTPUT_DIR/alignment"
EXPR_DIR="$OUTPUT_DIR/expression"
VAR_DIR="$OUTPUT_DIR/variants"
NEO_DIR="$OUTPUT_DIR/neoantigens"
RESOURCE_DIR="$OUTPUT_DIR/resources"
LOG_DIR="$OUTPUT_DIR/logs/${RUN_DATE}"

# =================================================================================
# SAMPLE CONFIGURATION
# =================================================================================
ALL_SAMPLES=("ncbi_control" "Hepa1_6")
declare -A TUMOR_NORMAL_PAIRS=(["Hepa1_6"]="ncbi_control")
TUMOR_SAMPLES=("Hepa1_6")
NORMAL_SAMPLES=("ncbi_control")

# =================================================================================
# INPUT DATA PATHS
# =================================================================================
declare -A READ_PATHS=(
    ["ncbi_control"]="/home/worst/projects/rnaseq_project/00_data/raw_reads/ncbi_control"
    ["Hepa1_6"]="/home/worst/projects/rnaseq_project/00_data/raw_reads/Hepa1_6"
)
declare -A READ1_FILENAMES=(
    ["ncbi_control"]="SRR30981057_1.fastq.gz"
    ["Hepa1_6"]="Lib037_Hepa1_6.R1.trimmed.fastq.gz"
)
declare -A READ2_FILENAMES=(
    ["ncbi_control"]="SRR30981057_2.fastq.gz"
    ["Hepa1_6"]="Lib037_Hepa1_6.R2.trimmed.fastq.gz"
)

# =================================================================================
# COMPUTATIONAL PARAMETERS
# =================================================================================
THREADS=16
JAVA_MEM="16g"

# =================================================================================
# PIPELINE BEHAVIOR
# =================================================================================
KEEP_INTERMEDIATE=true    # Keep intermediate files (SAM, unsorted BAMs, etc.)
SKIP_COMPLETED=true       # Skip steps that have already been completed
VERBOSE_LOGGING=true      # Detailed logging

# =================================================================================
# STAR ALIGNMENT PARAMETERS
# =================================================================================
STAR_OVERHANG=149
STAR_OUTFILTERMULTIMAPNMAX=20
STAR_OUTFILTERMISMATCHNMAX=999
STAR_OUTFILTERMISMATCHNOVERLMAX=0.04
STAR_ALIGNINTRONMIN=20
STAR_ALIGNINTRONMAX=1000000
STAR_ALIGNMATESGAPMAX=1000000
STAR_ALIGNSJOVERHANGMIN=8
STAR_ALIGNSJDBOVERHANGMIN=1

# =================================================================================
# EXPRESSION QUANTIFICATION
# =================================================================================
KALLISTO_BOOTSTRAPS=30

# =================================================================================
# VARIANT CALLING PARAMETERS
# =================================================================================
# General thresholds
MIN_COVERAGE_TUMOR=10
MIN_COVERAGE_NORMAL=10
MIN_AD=3
MIN_VAF_TUMOR=0.05
MAX_VAF_NORMAL=0.02
AF_OF_ALLELES_NOT_IN_RESOURCE=0.0000025

# VarScan2 specific
VARSCAN_MIN_COV_NORMAL=8
VARSCAN_MIN_COV_TUMOR=10
VARSCAN_SOMATIC_PVALUE=0.05
VARSCAN_SOMATIC_SOMATICPVALUE=0.05
VARSCAN_STRAND_FILTER=1

# Primary caller for high-confidence calling
PRIMARY_CALLER="M2"

# =================================================================================
# TOOL PATHS
# =================================================================================
GATK3_JAR="/home/worst/projects/universal_data/tools/GenomeAnalysisTK.jar"

# Picard discovery
if [ -n "$CONDA_PREFIX" ]; then
    PICARD_JAR=$(find "$CONDA_PREFIX/share" -name "picard.jar" 2>/dev/null | head -1)
fi
if [ -z "$PICARD_JAR" ] || [ ! -f "$PICARD_JAR" ]; then
    PICARD_JAR=$(find "$HOME/miniconda3/envs/variant_env/share" -name "picard.jar" 2>/dev/null | head -1)
fi

# =================================================================================
# VEP ANNOTATION
# =================================================================================
VEP_CACHE_DIR="${UNIVERSAL_DATA_DIR}/vep/vep_cache"
VEP_PLUGINS_DIR="${UNIVERSAL_DATA_DIR}/vep/VEP_plugins"
VEP_ASSEMBLY="GRCm39"
VEP_SPECIES="mus_musculus"
VEP_CACHE_VERSION="105"

# =================================================================================
# PVACSEQ PARAMETERS
# =================================================================================
MHC_ALLELES="H-2-Db,H-2-Kb"
PVAC_ALGORITHMS="MHCflurry"
PEPTIDE_LENGTHS="8,9,10,11"
BINDING_THRESHOLD=500
PVACSEQ_THREADS=2

# Expression and coverage filters
EXPRESSION_CUTOFF=1.0
TUMOR_RNA_VAF=0.05
TUMOR_RNA_COV=10
NORMAL_COV=10
NORMAL_VAF=0.02
MAX_TSL=1

# =================================================================================
# ENVIRONMENT SETTINGS
# =================================================================================
export TF_CPP_MIN_LOG_LEVEL='3'
export PYTHONWARNINGS="ignore:UserWarning,ignore:FutureWarning"