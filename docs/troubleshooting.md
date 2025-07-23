# Troubleshooting Guide

## Common Issues and Solutions

### Installation Issues

#### Conda environment creation fails

**Error**: `ResolvePackageNotFound` or dependency conflicts

**Solutions**:
1. Update conda:
   ```bash
   conda update -n base conda
   ```

2. Use mamba (faster resolver):
   ```bash
   conda install -n base mamba
   mamba env create -f envs/variant_env.yml
   ```

3. Create environment step by step:
   ```bash
   conda create -n variant_env python=3.9
   conda activate variant_env
   conda install -c bioconda -c conda-forge gatk4 star samtools
   ```

#### VEP plugins not found

**Error**: `Can't locate Frameshift.pm in @INC`

**Solution**:
```bash
# Verify plugin location
ls /home/worst/projects/universal_data/vep/VEP_plugins/

# Re-download if missing
cd /home/worst/projects/universal_data/vep/VEP_plugins
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/Frameshift.pm
wget https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/105/Wildtype.pm
```

### Runtime Errors

#### Out of memory during STAR alignment

**Error**: `EXITING because of FATAL ERROR: not enough memory`

**Solutions**:
1. Reduce thread count in config.sh:
   ```bash
   THREADS=8  # Instead of 16
   ```

2. Clear system memory:
   ```bash
   sync && echo 3 > /proc/sys/vm/drop_caches
   ```

3. Use --genomeLoad NoSharedMemory option (add to script)

#### Java heap space error

**Error**: `java.lang.OutOfMemoryError: Java heap space`

**Solution**:
Increase Java memory in config.sh:
```bash
JAVA_MEM="24g"  # Increase from 16g
```

#### GATK3 not found

**Error**: `ERROR: Java 8 executable not found`

**Solution**:
```bash
# Verify Java 8 in analysis_env
conda activate analysis_env
which java
java -version  # Should show 1.8.x

# Update path in script if needed
JAVA8_PATH="$CONDA_PREFIX/bin/java"
```

### Variant Calling Issues

#### No variants called

**Possible causes**:
1. Low coverage
2. Strict filtering
3. Sample swap

**Debugging**:
```bash
# Check BAM coverage
samtools depth alignment/processed/tumor/tumor.analysis_ready.bam | \
  awk '{sum+=$3} END {print "Average coverage:", sum/NR}'

# Check raw variant calls before filtering
bcftools view -H variants/tumor/raw_calls/mutect2/tumor.raw.vcf.gz | wc -l

# Verify sample pairing
samtools view -H alignment/processed/tumor/tumor.analysis_ready.bam | grep @RG
```

#### VarScan mpileup empty

**Error**: `No output from VarScan`

**Solution**:
```bash
# Test mpileup generation manually
samtools mpileup -B -f $REF_GENOME_FA $NORMAL_BAM $TUMOR_BAM | head

# Check BAM file integrity
samtools quickcheck $TUMOR_BAM
```

### Expression Analysis Issues

#### Kallisto index build fails

**Error**: `Error: could not open file`

**Solution**:
```bash
# Verify transcriptome file
zcat $REF_TRANSCRIPTOME_FA | head

# Rebuild index with full path
kallisto index -i /full/path/to/index /full/path/to/transcriptome.fa.gz
```

#### Zero gene expression

**Possible causes**:
1. GTF/transcriptome mismatch
2. Wrong species
3. Failed Kallisto run

**Debugging**:
```bash
# Check Kallisto output
cat expression/kallisto/sample/run_info.json

# Verify GTF format
head -20 $REF_GENOME_GTF

# Test R script manually
conda activate analysis_env
Rscript helpers/summarize_kallisto.R --help
```

### pVACseq Issues

#### No epitopes predicted

**Common causes**:
1. No expressed variants
2. Wrong MHC alleles
3. Too strict filtering

**Solutions**:
```bash
# Check input VCF has variants
bcftools view -H variants/tumor/annotated/tumor.final.annotated.vcf.gz | wc -l

# Verify expression annotation
bcftools query -f '%INFO/GX\n' variants/tumor/annotated/tumor.final.annotated.vcf.gz | \
  grep -v "^\\.$" | head

# Run with relaxed filters
# In config.sh:
BINDING_THRESHOLD=1000  # Instead of 500
EXPRESSION_CUTOFF=0.1   # Instead of 1.0
```

#### MHCflurry models not found

**Error**: `FileNotFoundError: MHCflurry models not found`

**Solution**:
```bash
conda activate pvacseq_env
mhcflurry-downloads fetch models_class1_presentation
mhcflurry-downloads info
```

### File Format Issues

#### VCF header errors

**Error**: `Contig 'chr1' not found in reference`

**Solution**:
```bash
# Standardize VCF headers
bcftools reheader -f ${REF_GENOME_FA}.fai input.vcf.gz -o output.vcf.gz
```

#### Compressed file errors

**Error**: `[E::hts_open_format] Failed to open file`

**Solutions**:
```bash
# Recompress file
gunzip -c file.vcf.gz | bgzip > file.new.vcf.gz
mv file.new.vcf.gz file.vcf.gz

# Reindex
bcftools index -f file.vcf.gz
```

### Performance Issues

#### Pipeline running slowly

**Optimization strategies**:

1. Check disk I/O:
   ```bash
   iostat -x 1
   ```

2. Use local disk instead of network storage

3. Reduce parallel jobs if memory constrained

4. Skip unnecessary steps:
   ```bash
   ./master_pipeline.sh --start-from variant_calling
   ```

#### Hanging processes

**Diagnosis**:
```bash
# Check CPU usage
top -u $USER

# Check specific process
ps aux | grep STAR

# Monitor memory
free -h
```

**Solution**:
Kill stuck process and restart:
```bash
kill -9 <PID>
rm pipeline_output/alignment/processed/sample/.done
./master_pipeline.sh --sample sample --start-from alignment
```

## Error Diagnosis Checklist

When encountering errors:

1. **Check the log files**:
   ```bash
   tail -50 pipeline_output/logs/*/04_variant_calling.log
   ```

2. **Verify input files exist**:
   ```bash
   ls -la alignment/processed/*/
   ```

3. **Check available disk space**:
   ```bash
   df -h
   ```

4. **Verify environment activation**:
   ```bash
   conda info --envs
   which gatk
   ```

5. **Test individual commands**:
   Copy the failing command from the log and run manually

## Getting Help

If issues persist:

1. Check tool-specific documentation:
   - GATK forums: https://gatk.broadinstitute.org/
   - pVACtools: https://pvactools.readthedocs.io/
   - Biostars: https://www.biostars.org/

2. Verify versions match documentation:
   ```bash
   gatk --version
   pvacseq --version
   vep --version
   ```

3. Create minimal test case to isolate issue

4. Check GitHub issues for similar problems