# Installation and Setup Guide

## Prerequisites

### System Requirements
- Linux operating system (tested on WSL Ubuntu 20.04)
- Minimum 32GB RAM (64GB recommended)
- Minimum 500GB free disk space
- Internet connection for downloading tools and references

### Required Software
- Conda or Mamba package manager
- Git (for cloning the repository)

## Step 1: Clone Repository

```bash
git clone https://github.com/yourusername/rnaseq-neoantigen-pipeline.git
cd rnaseq-neoantigen-pipeline
```

## Step 2: Create Conda Environments

The pipeline uses multiple conda environments to manage dependencies:

```bash
# Create all environments
conda env create -f envs/qc_env.yml
conda env create -f envs/variant_env.yml
conda env create -f envs/analysis_env.yml
conda env create -f envs/pvacseq_env.yml
conda env create -f envs/pvacdownstream_env.yml
```

This will take 30-60 minutes depending on your internet connection.

## Step 3: Download GATK3

GATK3 must be downloaded manually due to licensing:

1. Go to https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk
2. Download `GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2`
3. Extract and place the jar file:

```bash
tar -xjf GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2
mkdir -p /home/worst/projects/universal_data/tools
cp GenomeAnalysisTK.jar /home/worst/projects/universal_data/tools/
```

## Step 4: Download Reference Files

### Mouse Reference Genome (GRCm39)

```bash
# Create directories
mkdir -p /home/worst/projects/universal_data/reference/Mus_musculus.GRCm39.105

# Download reference files
cd /home/worst/projects/universal_data/reference/Mus_musculus.GRCm39.105
wget https://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
wget https://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

# Decompress genome (keep others compressed)
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

### Mouse Genome Project Variants

```bash
mkdir -p /home/worst/projects/universal_data/known_sites/mousegenomeproject
cd /home/worst/projects/universal_data/known_sites/mousegenomeproject

# Download from Mouse Genome Project
# These files need to be obtained from MGP website
# Place the following files here:
# - mgp_REL2021_snps.rsID.vcf.gz
# - mgp_REL2021_indels.rsID.vcf.gz
# - mgp_REL2021_common_snps.vcf.gz
```

## Step 5: Install VEP Cache and Plugins

```bash
# Activate pvacseq environment
conda activate pvacseq_env

# Create VEP directories
mkdir -p /home/worst/projects/universal_data/vep/vep_cache
mkdir -p /home/worst/projects/universal_data/vep/VEP_plugins

# Download VEP cache for mouse
cd /home/worst/projects/universal_data/vep
curl -O https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/mus_musculus_vep_105_GRCm39.tar.gz
tar -xzf mus_musculus_vep_105_GRCm39.tar.gz -C vep_cache/

# Download VEP plugins
cd VEP_plugins
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/Frameshift.pm
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/Wildtype.pm

conda deactivate
```

## Step 6: Download MHCflurry Models

```bash
conda activate pvacseq_env
mhcflurry-downloads fetch models_class1_presentation
conda deactivate
```

## Step 7: Set Up Directory Structure

```bash
# Create project directories
mkdir -p /home/worst/projects/rnaseq_project/pipeline_output
mkdir -p /home/worst/projects/rnaseq_project/00_data/raw_reads/ncbi_control
mkdir -p /home/worst/projects/rnaseq_project/00_data/raw_reads/Hepa1_6
```

## Step 8: Place Your Data

Copy your FASTQ files to the appropriate directories:
- Control sample: `/home/worst/projects/rnaseq_project/00_data/raw_reads/ncbi_control/`
- Tumor sample: `/home/worst/projects/rnaseq_project/00_data/raw_reads/Hepa1_6/`

## Step 9: Make Scripts Executable

```bash
cd /path/to/rnaseq-neoantigen-pipeline
chmod +x *.sh
chmod +x helpers/*.py
chmod +x helpers/*.R
```

## Step 10: Test Installation

Run a quick test to ensure everything is set up correctly:

```bash
# Test conda environments
conda activate qc_env && fastqc --version && conda deactivate
conda activate variant_env && gatk --version && conda deactivate
conda activate pvacseq_env && pvacseq --version && conda deactivate

# Test reference files
ls -la /home/worst/projects/universal_data/reference/Mus_musculus.GRCm39.105/
ls -la /home/worst/projects/universal_data/tools/GenomeAnalysisTK.jar
```

## Troubleshooting

### Conda environment creation fails
- Try using mamba instead: `mamba env create -f envs/xxx.yml`
- Update conda: `conda update -n base conda`

### VEP cache download fails
- Check disk space
- Try downloading via browser and transferring manually

### Out of memory errors
- Reduce thread count in config.sh
- Increase Java heap size if system allows

## Next Steps

After successful setup:
1. Review and adjust parameters in `config.sh`
2. Run initial QC: `./01_quality_control.sh`
3. If QC looks good, run full pipeline: `./master_pipeline.sh`