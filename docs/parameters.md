# Configuration Parameters

This document describes all configurable parameters in `config.sh`.

## Sample Configuration

### ALL_SAMPLES
Array of all sample IDs to process.
```bash
ALL_SAMPLES=("ncbi_control" "Hepa1_6")
```

### TUMOR_NORMAL_PAIRS
Associates tumor samples with their matched normals.
```bash
declare -A TUMOR_NORMAL_PAIRS=(["Hepa1_6"]="ncbi_control")
```

### READ_PATHS
Paths to directories containing FASTQ files.
```bash
declare -A READ_PATHS=(
    ["sample_id"]="/path/to/fastq/directory"
)
```

## Computational Parameters

### THREADS
Number of CPU threads to use for parallel processing.
- Default: 16
- Reduce if memory limited

### JAVA_MEM
Java heap memory allocation.
- Default: "16g"
- Format: NUMBER[g|m] (g=gigabytes, m=megabytes)

## Pipeline Behavior

### KEEP_INTERMEDIATE
Whether to retain intermediate files.
- Default: true
- Set false to save disk space

### SKIP_COMPLETED
Whether to skip already completed steps.
- Default: true
- Set false to force re-run

### VERBOSE_LOGGING
Enable detailed logging output.
- Default: true

## STAR Alignment Parameters

### STAR_OVERHANG
Read length minus 1 for splice junction database.
- Default: 149 (for 150bp reads)
- Formula: read_length - 1

### STAR_OUTFILTERMULTIMAPNMAX
Maximum number of multiple alignments allowed.
- Default: 20
- Standard for ENCODE

### STAR_OUTFILTERMISMATCHNMAX
Maximum number of mismatches per pair.
- Default: 999
- Effectively no limit

### STAR_OUTFILTERMISMATCHNOVERLMAX
Maximum ratio of mismatches to mapped length.
- Default: 0.04
- 4% mismatch rate allowed

## Expression Quantification

### KALLISTO_BOOTSTRAPS
Number of bootstrap samples for abundance estimation.
- Default: 30
- Higher values increase accuracy but take longer

## Variant Calling Parameters

### Coverage Thresholds

#### MIN_COVERAGE_TUMOR
Minimum coverage in tumor sample.
- Default: 10
- Standard for RNA-seq

#### MIN_COVERAGE_NORMAL
Minimum coverage in normal sample.
- Default: 10

#### MIN_AD
Minimum alternate allele depth.
- Default: 3
- Ensures variant support

### VAF Thresholds

#### MIN_VAF_TUMOR
Minimum variant allele frequency in tumor.
- Default: 0.05 (5%)
- Balanced for RNA-seq sensitivity

#### MAX_VAF_NORMAL
Maximum VAF in normal to filter germline.
- Default: 0.02 (2%)
- Filters contamination

### Mutect2 Specific

#### AF_OF_ALLELES_NOT_IN_RESOURCE
Prior probability of novel variants.
- Default: 0.0000025
- Very rare variants

### VarScan2 Specific

#### VARSCAN_MIN_COV_NORMAL
Minimum coverage for normal sample.
- Default: 8

#### VARSCAN_MIN_COV_TUMOR
Minimum coverage for tumor sample.
- Default: 10

#### VARSCAN_SOMATIC_PVALUE
P-value threshold for somatic calls.
- Default: 0.05

#### VARSCAN_STRAND_FILTER
Enable strand bias filter.
- Default: 1 (enabled)

### Variant Merging

#### PRIMARY_CALLER
Which caller to use as primary for merging.
- Default: "M2" (Mutect2)
- Options: "M2", "VS"

## VEP Annotation

### VEP_ASSEMBLY
Reference genome assembly.
- Default: "GRCm39"

### VEP_SPECIES
Species for annotation.
- Default: "mus_musculus"

### VEP_CACHE_VERSION
VEP cache version.
- Default: "105"
- Must match downloaded cache

## pVACseq Parameters

### MHC_ALLELES
MHC alleles for neoantigen prediction.
- Default: "H-2-Db,H-2-Kb" (mouse)
- Format: comma-separated list

### PVAC_ALGORITHMS
Prediction algorithms to use.
- Default: "MHCflurry"
- Options: MHCflurry, MHCflurryEL, NetMHCpan, others

### PEPTIDE_LENGTHS
Epitope lengths to predict.
- Default: "8,9,10,11"
- MHC Class I typical lengths

### BINDING_THRESHOLD
IC50 binding affinity threshold (nM).
- Default: 500
- Standard cutoff for binders

### PVACSEQ_THREADS
Threads for pVACseq (memory intensive).
- Default: 2
- Increase carefully

### Expression Filters

#### EXPRESSION_CUTOFF
Minimum gene expression (TPM).
- Default: 1.0
- Filters non-expressed genes

#### TUMOR_RNA_VAF
Minimum VAF in tumor RNA.
- Default: 0.05 (5%)

#### TUMOR_RNA_COV
Minimum coverage in tumor RNA.
- Default: 10

#### NORMAL_COV
Coverage threshold for normal.
- Default: 10

#### NORMAL_VAF
VAF threshold for normal.
- Default: 0.02 (2%)

#### MAX_TSL
Maximum transcript support level.
- Default: 1
- Ensures well-supported transcripts

## Performance Tuning

### Memory-Intensive Steps

1. **STAR Alignment**: ~30GB for mouse genome
2. **GATK HaplotypeCaller**: Set by JAVA_MEM
3. **pVACseq**: ~4GB per thread

### Disk Space Requirements

With KEEP_INTERMEDIATE=true:
- ~10GB per sample for BAM files
- ~1GB per sample for VCF files
- ~500MB per sample for expression data

### Speed Optimization

To speed up analysis:
1. Increase THREADS (if memory allows)
2. Reduce KALLISTO_BOOTSTRAPS
3. Use fewer PEPTIDE_LENGTHS
4. Limit PVAC_ALGORITHMS to one

## Modifying Parameters

### When to Adjust

Modify parameters when:
- Using different read lengths (adjust STAR_OVERHANG)
- Analyzing different species (update MHC_ALLELES)
- Working with low-coverage data (reduce thresholds)
- Memory constraints (reduce THREADS, JAVA_MEM)

### Best Practices

1. Always backup original config.sh
2. Change one parameter at a time
3. Document your changes
4. Re-run from appropriate step

### Parameter Validation

The pipeline validates critical parameters:
- File paths must exist
- Numeric values must be in valid ranges
- Required tools must be accessible

## Default Configuration Rationale

The default parameters are optimized for:
- Mouse RNA-seq data
- 150bp paired-end reads
- Tumor-normal pairs
- Standard computing resources (16 cores, 32GB RAM)
- Balance between sensitivity and specificity