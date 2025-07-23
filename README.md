# RNA-seq Neoantigen Pipeline v1.0 (incomplete)

## Overview

This pipeline processes RNA-seq data to identify somatic variants and predict neoantigens for cancer immunotherapy research. It integrates multiple variant callers, performs read-backed phasing, and uses pVACseq for neoantigen prediction. This is not yet complete, the configs and file paths are local to my desktop. The uploaded files serve as a backup copy. Do not use.

## Key Features

- **Dual variant calling** with Mutect2 and VarScan2 for high confidence
- **RNA-specific optimizations** throughout the workflow
- **Read-backed phasing** for accurate neoantigen sequences
- **Expression-based filtering** using Kallisto quantification
- **Interactive visualization** with pVACview

## Pipeline Steps

1. **Quality Control** - FastQC analysis of raw reads
2. **Alignment** - STAR 2-pass alignment with GATK preprocessing
3. **Expression Analysis** - Kallisto quantification to gene-level TPM
4. **Variant Calling** - Somatic variant detection with multiple callers
5. **Germline Calling** - HaplotypeCaller for phasing information
6. **Variant Processing** - Merging, normalization, and VEP annotation
7. **Phasing** - ReadBackedPhasing for proximal variant correction
8. **Neoantigen Prediction** - pVACseq analysis with MHCflurry
9. **Visualization** - Interactive review with pVACview

## Quick Start

```bash
# Run full pipeline
./master_pipeline.sh

# Run specific sample
./master_pipeline.sh --sample Hepa1_6

# Run from specific step
./master_pipeline.sh --start-from variant_calling

# Run only certain steps
./master_pipeline.sh --only expression,neoantigen_prediction
```

## Requirements

- Linux system (tested on WSL Ubuntu)
- Conda package manager
- At least 32GB RAM
- GATK3 jar file (manual download required)

## Output Structure

```
pipeline_output/
├── QC/                    # Quality control reports
├── alignment/             # STAR alignments and GATK-processed BAMs
├── expression/            # Kallisto quantification results
├── variants/              # Variant calling outputs
│   └── {sample}/
│       ├── raw_calls/     # Individual caller outputs
│       ├── merged/        # High-confidence variants
│       ├── germline/      # Germline variants for phasing
│       ├── phased/        # Phased variants
│       └── annotated/     # VEP-annotated VCFs
├── neoantigens/          # pVACseq predictions
└── logs/                 # Pipeline execution logs
```

## Configuration

All pipeline parameters are defined in `config.sh`:
- Sample definitions and pairing
- Reference file paths
- Tool parameters
- Filtering thresholds
- MHC alleles for prediction

## Documentation

- [SETUP.md](docs/setup.md) - Detailed installation instructions
- [USAGE.md](docs/usage.md) - Running the pipeline
- [PARAMETERS.md](docs/parameters.md) - Configuration options
- [OUTPUT.md](docs/output.md) - Understanding the results

## Citation

This pipeline integrates multiple published tools. Please cite:
- STAR: Dobin et al. (2013) Bioinformatics
- GATK: McKenna et al. (2010) Genome Research
- Mutect2: Benjamin et al. (2019) bioRxiv
- VarScan2: Koboldt et al. (2012) Genome Research
- VEP: McLaren et al. (2016) Genome Biology
- pVACseq: Hundal et al. (2020) Cancer Immunology Research
- MHCflurry: O'Donnell et al. (2018) Cell Systems
- Reference nextNEOpi's code

## Author

Prototype created for RNA-seq variant analysis and neoantigen discovery research.

## Version

Pipeline Version: 1.0
Release Date: 2025
