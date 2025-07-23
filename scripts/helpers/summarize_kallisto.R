#!/usr/bin/env Rscript
# =================================================================================
# Helper Script: Summarize Kallisto to Gene-Level TPM
# Version 3.0 - Simplified for expression annotation
# =================================================================================

# --- Load Libraries ---
suppressPackageStartupMessages({
    library(argparse)
    library(tximport)
    library(GenomicFeatures)
    library(txdbmaker)
    library(readr)
})

# --- Argument Parsing ---
parser <- ArgumentParser(description="Convert Kallisto transcript TPMs to gene-level TPMs")
parser$add_argument("--gtf", type="character", required=TRUE, 
                    help="Path to the reference GTF file")
parser$add_argument("--kallisto_dir", type="character", required=TRUE,
                    help="Directory containing Kallisto abundance.tsv file")
parser$add_argument("--sample_name", type="character", required=TRUE,
                    help="Sample name for output")
parser$add_argument("--output_dir", type="character", required=TRUE,
                    help="Directory to save output files")

args <- parser$parse_args()

# Create output directory if needed
dir.create(args$output_dir, showWarnings=FALSE, recursive=TRUE)

# =================================================================================
# Create transcript-to-gene mapping
# =================================================================================
cat("Creating transcript-to-gene mapping from GTF...\n")
tryCatch({
    # Suppress warnings about unsupported metadata columns
    suppressWarnings({
        txdb <- makeTxDbFromGFF(args$gtf, format="gtf")
    })
    
    k <- keys(txdb, keytype="TXNAME")
    tx2gene <- select(txdb, keys=k, columns="GENEID", keytype="TXNAME")
    colnames(tx2gene) <- c("TXNAME", "GENEID")
    
    cat(sprintf("✓ Created mapping for %d transcripts.\n", nrow(tx2gene)))
}, error = function(e) {
    cat(sprintf("ERROR: Failed to create TxDb from GTF: %s\n", e$message))
    quit(status=1)
})

# =================================================================================
# Import Kallisto results
# =================================================================================
cat("Importing Kallisto results...\n")
kallisto_file <- file.path(args$kallisto_dir, "abundance.tsv")

if (!file.exists(kallisto_file)) {
    cat(sprintf("ERROR: Kallisto file not found: %s\n", kallisto_file))
    quit(status=1)
}

# Import and aggregate to gene level
files <- kallisto_file
names(files) <- args$sample_name

txi <- tximport(files, type="kallisto", tx2gene=tx2gene, 
                ignoreTxVersion=TRUE, txOut=FALSE)

# =================================================================================
# Write gene-level TPM
# =================================================================================
cat("Writing gene-level TPM file...\n")
gene_tpm <- data.frame(
    gene_id = rownames(txi$abundance),
    tpm = round(txi$abundance[, args$sample_name], 6)
)

output_file <- file.path(args$output_dir, 
                        sprintf("%s_gene_tpm.tsv", args$sample_name))
write_tsv(gene_tpm, output_file)
cat(sprintf("✓ Gene-level TPM written to: %s\n", output_file))

# Summary statistics
cat(sprintf("\nSummary:\n"))
cat(sprintf("  Total genes: %d\n", nrow(gene_tpm)))
cat(sprintf("  Expressed genes (TPM > 1): %d\n", sum(gene_tpm$tpm > 1)))
cat(sprintf("  Median TPM (non-zero): %.2f\n", 
            median(gene_tpm$tpm[gene_tpm$tpm > 0])))

cat("\n✓ Complete\n")