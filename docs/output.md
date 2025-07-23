# Understanding Pipeline Output

## Output Directory Structure

```
pipeline_output/
├── QC/                          # Quality control reports
│   ├── fastqc/                  # FastQC reports
│   │   └── {sample}/            # Per-sample QC
│   ├── qualimap/                # Post-alignment QC
│   │   └── {sample}/            # RNA-seq specific metrics
│   ├── bqsr_plots/              # Base quality recalibration plots
│   └── multiqc_report.html      # Combined QC summary
│
├── alignment/                   # Alignment files
│   ├── star/                    # STAR aligner output
│   │   └── {sample}/
│   │       ├── Aligned.sortedByCoord.out.bam
│   │       ├── Log.final.out    # Alignment statistics
│   │       └── SJ.out.tab       # Splice junctions
│   └── processed/               # GATK-processed BAMs
│       └── {sample}/
│           ├── {sample}.analysis_ready.bam      # Final BAM
│           ├── {sample}.analysis_ready.bam.bai  # BAM index
│           └── recal.table      # BQSR table
│
├── expression/                  # Expression quantification
│   ├── kallisto/                # Raw Kallisto output
│   │   └── {sample}/
│   │       ├── abundance.tsv    # Transcript-level
│   │       └── run_info.json    # Run metadata
│   └── gene_tpm/                # Gene-level summaries
│       └── {sample}_gene_tpm.tsv
│
├── variants/                    # Variant calling results
│   └── {tumor_sample}/
│       ├── raw_calls/           # Individual caller outputs
│       │   ├── mutect2/         # Mutect2 results
│       │   └── varscan/         # VarScan2 results
│       ├── merged/              # Merged high-confidence calls
│       │   ├── {sample}.hc.vcf.gz              # Confirmed variants
│       │   └── {sample}.hc.normalized.vcf.gz   # Normalized
│       ├── germline/            # Germline variants
│       │   └── {normal}.genotyped_germline.vcf.gz
│       ├── phased/              # Phased variants
│       │   └── {sample}.phased.annotated.vcf.gz
│       └── annotated/           # Fully annotated
│           └── {sample}.final.annotated.vcf.gz
│
├── neoantigens/                 # Neoantigen predictions
│   └── {tumor_sample}/
│       ├── pvacseq/             # Raw pVACseq output
│       │   └── MHC_Class_I/
│       │       └── {sample}.all_epitopes.tsv
│       └── reports/             # Processed reports
│           ├── {sample}.all_epitopes.aggregated.tsv
│           └── {sample}.all_epitopes.aggregated.metrics.json
│
├── resources/                   # Reusable indices
│   ├── star_index/              # STAR genome index
│   └── kallisto_index/          # Kallisto transcriptome index
│
└── logs/                        # Execution logs
    └── {timestamp}/
        ├── master_pipeline.log
        └── {step}.log
```

## Key Output Files

### Quality Control

#### FastQC Reports
- **Location**: `QC/fastqc/{sample}/`
- **Key files**: 
  - `*_fastqc.html` - Interactive QC report
  - `*_fastqc.zip` - Raw QC data
- **Important metrics**:
  - Per base sequence quality
  - Adapter content
  - Overrepresented sequences

#### Qualimap Reports
- **Location**: `QC/qualimap/{sample}/`
- **Key file**: `qualimapReport.html`
- **Important metrics**:
  - Reads alignment statistics
  - Gene coverage profile
  - Junction analysis

### Alignment Files

#### Final BAM
- **File**: `alignment/processed/{sample}/{sample}.analysis_ready.bam`
- **Description**: GATK-processed, analysis-ready alignment
- **Processing**: 
  - Duplicates marked
  - Base qualities recalibrated
  - Split N cigar reads
- **Use for**: Variant calling, visualization in IGV

#### STAR Log
- **File**: `alignment/star/{sample}/Log.final.out`
- **Key metrics**:
  - Total/unique/multi-mapped reads
  - Mismatch and deletion rates
  - Splicing statistics

### Expression Data

#### Gene-level TPM
- **File**: `expression/gene_tpm/{sample}_gene_tpm.tsv`
- **Format**:
  ```
  gene_id         tpm
  ENSMUSG00001    45.231
  ```
- **Use for**: Filtering variants, differential expression

### Variant Files

#### High-confidence Variants
- **File**: `variants/{sample}/merged/{sample}.hc.vcf.gz`
- **Description**: Variants called by multiple callers
- **INFO field**: `VariantCalledBy` shows which callers

#### Final Annotated VCF
- **File**: `variants/{sample}/annotated/{sample}.final.annotated.vcf.gz`
- **Annotations include**:
  - VEP consequences
  - Gene expression (GX tag)
  - RNA read counts
  - Allele frequencies

#### VCF INFO/FORMAT Fields

Key custom fields added by pipeline:
- `GX`: Gene expression in TPM
- `RAF`: RNA allele frequency
- `RAD`: RNA allele depth
- `RDP`: RNA total depth

### Neoantigen Predictions

#### Aggregated Report
- **File**: `neoantigens/{sample}/reports/{sample}.all_epitopes.aggregated.tsv`
- **Key columns**:
  - `Gene Name`: Source gene
  - `Mutation`: Genomic change
  - `HLA Allele`: MHC allele
  - `Peptide Length`: Epitope length
  - `MT Epitope`: Mutant peptide sequence
  - `WT Epitope`: Wild-type peptide
  - `Median MT Score`: Binding affinity (nM)
  - `Median WT Score`: Wild-type affinity
  - `Median Fold Change`: MT/WT ratio
  - `Tumor RNA Depth`: Coverage at variant
  - `Tumor RNA VAF`: Variant frequency
  - `Gene Expression`: TPM value

#### Filtering Recommendations

Good neoantigen candidates typically have:
- Median MT Score < 500 nM (strong binders < 50 nM)
- Median Fold Change < 1 (better binding than WT)
- Gene Expression > 1 TPM
- Tumor RNA VAF > 0.1
- Tumor RNA Depth > 10

## Understanding VCF Annotations

### VEP Consequences

Common RNA-seq relevant consequences:
- `frameshift_variant`: Potential strong neoantigens
- `missense_variant`: Single amino acid changes
- `stop_gained`: Truncated proteins
- `splice_donor/acceptor_variant`: Altered splicing

### Expression Integration

Variants are annotated with expression data:
```
##INFO=<ID=GX,Number=.,Type=Float,Description="Gene expression in TPM">
1  12345  .  A  T  .  PASS  GX=25.5;CSQ=...
```

## Summary Report

The HTML summary report (`summary_report.html`) provides:
- Overall pipeline statistics
- Per-sample processing status
- Variant counts by caller
- Neoantigen prediction summary
- Links to detailed reports

## Tips for Results Interpretation

### Variant Prioritization

1. Check variant caller agreement:
   ```bash
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VariantCalledBy\n' \
     variants/tumor/merged/tumor.hc.vcf.gz
   ```

2. Filter for expressed variants:
   ```bash
   bcftools view -i 'GX>1' variants/tumor/annotated/tumor.final.annotated.vcf.gz
   ```

### Neoantigen Selection

1. Sort by binding affinity:
   ```bash
   sort -t$'\t' -k15,15n neoantigens/tumor/reports/tumor.all_epitopes.aggregated.tsv
   ```

2. Filter for strong binders with expression:
   ```bash
   awk -F'\t' '$15<50 && $21>10' neoantigens/tumor/reports/tumor.all_epitopes.aggregated.tsv
   ```

### Quality Checks

Always verify:
1. Alignment rate > 80%
2. Sufficient coverage at variant sites
3. Expression of neoantigen source genes
4. Proper phasing of proximal variants

## Downstream Analysis

Results can be used for:
- Neoantigen vaccine design
- T-cell epitope validation
- Tumor mutation burden calculation
- Expression-weighted neoantigen load
- Clonality analysis (with additional tools)