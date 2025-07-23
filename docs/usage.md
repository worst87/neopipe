# Usage Guide

## Running the Complete Pipeline

### Basic Usage

The simplest way to run the pipeline is to execute the master script:

```bash
./master_pipeline.sh
```

This will process all samples defined in `config.sh` through all pipeline steps.

### Running Individual Steps

You can also run individual steps:

```bash
# Step 1: Quality control only
./01_quality_control.sh

# Step 2: Alignment and preprocessing
./02_alignment_preprocessing.sh

# Continue with other steps...
```

## Advanced Options

### Process Specific Sample

To process only one sample:

```bash
./master_pipeline.sh --sample Hepa1_6
```

### Start from Specific Step

If you need to restart from a particular step:

```bash
./master_pipeline.sh --start-from variant_calling
```

Available step names:
- `qc`
- `alignment`
- `expression`
- `variant_calling`
- `germline_calling`
- `variant_processing`
- `phasing`
- `neoantigen_prediction`
- `visualization`

### Run Specific Steps Only

To run only certain steps:

```bash
# Run only expression analysis and neoantigen prediction
./master_pipeline.sh --only expression,neoantigen_prediction
```

### Run Range of Steps

To run from one step to another:

```bash
# Run from variant calling through phasing
./master_pipeline.sh --start-from variant_calling --end-at phasing
```

## Monitoring Progress

### Log Files

All pipeline steps create detailed logs in:
```
pipeline_output/logs/YYYYMMDD_HHMMSS/
```

Monitor the master log:
```bash
tail -f pipeline_output/logs/*/master_pipeline.log
```

### Check Step Completion

The pipeline creates `.done` files to track completed steps:

```bash
# Check which samples have completed alignment
ls pipeline_output/alignment/processed/*/.done

# Check variant calling status
ls pipeline_output/variants/*/raw_calls/.done
```

### Generate Summary Report

At any time, generate an HTML summary:

```bash
./generate_summary.sh
```

View the report:
```bash
firefox pipeline_output/summary_report.html
```

## Interactive Analysis

### pVACview Visualization

After neoantigen prediction, launch interactive visualization:

```bash
./09_visualization.sh
```

This starts a web server. Follow the on-screen instructions to:
1. Open the provided URL in your browser
2. Load the aggregated TSV and JSON files
3. Review and filter neoantigen candidates
4. Export selected candidates

## Common Workflows

### Re-run After Parameter Changes

If you modify filtering parameters in `config.sh`:

```bash
# Re-run from variant processing onward
./master_pipeline.sh --start-from variant_processing
```

### Add New Sample

1. Add sample information to `config.sh`:
   - Add to `ALL_SAMPLES` array
   - Add to `READ_PATHS` and filename arrays
   - Add to `TUMOR_NORMAL_PAIRS` if applicable

2. Run pipeline for new sample only:
   ```bash
   ./master_pipeline.sh --sample new_sample_name
   ```

### Skip Expression Analysis

For non-RNA data or if expression data exists:

```bash
# Run all steps except expression
./master_pipeline.sh --only qc,alignment,variant_calling,germline_calling,variant_processing,phasing,neoantigen_prediction
```

## Output Files

### Key Output Locations

```bash
# Final annotated variants
pipeline_output/variants/{sample}/annotated/{sample}.final.annotated.vcf.gz

# Phased variants for pVACseq
pipeline_output/variants/{sample}/phased/{sample}.phased.annotated.vcf.gz

# Neoantigen predictions
pipeline_output/neoantigens/{sample}/reports/{sample}.all_epitopes.aggregated.tsv

# QC reports
pipeline_output/QC/multiqc_report.html
pipeline_output/QC/qualimap/{sample}/qualimapReport.html
```

### Understanding File Names

- `.hc.vcf.gz` - High-confidence variants (confirmed by multiple callers)
- `.normalized.vcf.gz` - Normalized variants (left-aligned, trimmed)
- `.vep.vcf` - VEP-annotated variants
- `.final.annotated.vcf.gz` - Fully annotated with expression and readcounts
- `.phased.annotated.vcf.gz` - Phased variants for proximal variant correction

## Performance Tips

### Memory Management

If running out of memory:
1. Reduce `THREADS` in `config.sh`
2. Process samples individually
3. Increase `JAVA_MEM` if system allows

### Disk Space

The pipeline generates large intermediate files. To save space:
1. Set `KEEP_INTERMEDIATE=false` in `config.sh`
2. Manually clean up after confirming results:
   ```bash
   # Remove STAR output after GATK processing
   rm -rf pipeline_output/alignment/star/*/
   ```

### Resume Failed Runs

The pipeline automatically skips completed steps. To force re-run:
```bash
# Remove .done file for specific step
rm pipeline_output/variants/Hepa1_6/raw_calls/.done

# Re-run
./master_pipeline.sh --sample Hepa1_6 --start-from variant_calling
```