#!/bin/bash
# =================================================================================
# RNA-SEQ NEOANTIGEN PIPELINE - MASTER SCRIPT
# Version 1.0
# =================================================================================

set -eo pipefail

# Parse command line arguments
SAMPLE=""
START_FROM=""
END_AT=""
STEPS_TO_RUN=""

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  --sample SAMPLE        Run pipeline for specific sample only"
    echo "  --start-from STEP      Start from specific step (e.g., variant_calling)"
    echo "  --end-at STEP          End at specific step"
    echo "  --only STEPS           Run only specific steps (comma-separated)"
    echo "  --help                 Show this help message"
    echo ""
    echo "Available steps:"
    echo "  qc, alignment, expression, variant_calling, germline_calling,"
    echo "  variant_processing, phasing, neoantigen_prediction, visualization"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE="$2"
            shift 2
            ;;
        --start-from)
            START_FROM="$2"
            shift 2
            ;;
        --end-at)
            END_AT="$2"
            shift 2
            ;;
        --only)
            STEPS_TO_RUN="$2"
            shift 2
            ;;
        --help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Source configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

# Create log directory
mkdir -p "$LOG_DIR"
MASTER_LOG="$LOG_DIR/master_pipeline.log"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MASTER_LOG"
}

# Error handling
error_exit() {
    log "ERROR: $1"
    exit 1
}

# Check if step should be run
should_run_step() {
    local step=$1
    
    # If specific steps requested, check if this step is included
    if [[ -n "$STEPS_TO_RUN" ]]; then
        if [[ ",$STEPS_TO_RUN," == *",$step,"* ]]; then
            return 0
        else
            return 1
        fi
    fi
    
    # Check start/end boundaries
    if [[ -n "$START_FROM" || -n "$END_AT" ]]; then
        local steps_order=(qc alignment expression variant_calling germline_calling variant_processing phasing neoantigen_prediction visualization)
        local start_idx=0
        local end_idx=${#steps_order[@]}
        local current_idx=-1
        
        for i in "${!steps_order[@]}"; do
            if [[ "${steps_order[$i]}" == "$step" ]]; then
                current_idx=$i
            fi
            if [[ "${steps_order[$i]}" == "$START_FROM" ]]; then
                start_idx=$i
            fi
            if [[ "${steps_order[$i]}" == "$END_AT" ]]; then
                end_idx=$i
            fi
        done
        
        if [[ $current_idx -ge $start_idx && $current_idx -le $end_idx ]]; then
            return 0
        else
            return 1
        fi
    fi
    
    return 0
}

# Check dependencies
check_dependencies() {
    log "Checking dependencies..."
    
    # Check if conda environments exist
    local envs=("qc_env" "variant_env" "analysis_env" "pvacseq_env" "pvacdownstream_env")
    for env in "${envs[@]}"; do
        if ! conda env list | grep -q "^$env "; then
            error_exit "Conda environment '$env' not found. Please run setup script."
        fi
    done
    
    # Check reference files
    [[ ! -f "$REF_GENOME_FA" ]] && error_exit "Reference genome not found: $REF_GENOME_FA"
    [[ ! -f "$REF_GENOME_GTF" ]] && error_exit "GTF file not found: $REF_GENOME_GTF"
    [[ ! -f "$GATK3_JAR" ]] && error_exit "GATK3 JAR not found: $GATK3_JAR"
    
    log "All dependencies satisfied."
}

# Main pipeline execution
main() {
    log "=========================================="
    log "RNA-SEQ NEOANTIGEN PIPELINE v${PIPELINE_VERSION}"
    log "=========================================="
    log "Run date: ${RUN_DATE}"
    log "Output directory: $OUTPUT_DIR"
    
    if [[ -n "$SAMPLE" ]]; then
        log "Running for sample: $SAMPLE"
    else
        log "Running for all samples"
    fi
    
    check_dependencies
    
    # Step 1: Quality Control
    if should_run_step "qc"; then
        log "Starting Step 1: Quality Control"
        "$SCRIPT_DIR/01_quality_control.sh" 2>&1 | tee -a "$LOG_DIR/01_qc.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Quality control failed"
    fi
    
    # Step 2: Alignment and Preprocessing
    if should_run_step "alignment"; then
        log "Starting Step 2: Alignment and Preprocessing"
        "$SCRIPT_DIR/02_alignment_preprocessing.sh" 2>&1 | tee -a "$LOG_DIR/02_alignment.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Alignment failed"
    fi
    
    # Step 3: Expression Analysis (RNA-specific)
    if should_run_step "expression"; then
        log "Starting Step 3: Expression Analysis"
        "$SCRIPT_DIR/03_expression_analysis.sh" 2>&1 | tee -a "$LOG_DIR/03_expression.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Expression analysis failed"
    fi
    
    # Step 4: Variant Calling
    if should_run_step "variant_calling"; then
        log "Starting Step 4: Variant Calling"
        "$SCRIPT_DIR/04_variant_calling.sh" 2>&1 | tee -a "$LOG_DIR/04_variant_calling.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Variant calling failed"
    fi
    
    # Step 5: Germline Calling
    if should_run_step "germline_calling"; then
        log "Starting Step 5: Germline Calling for Phasing"
        "$SCRIPT_DIR/05_germline_calling.sh" 2>&1 | tee -a "$LOG_DIR/05_germline.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Germline calling failed"
    fi
    
    # Step 6: Variant Processing
    if should_run_step "variant_processing"; then
        log "Starting Step 6: Variant Processing and Annotation"
        "$SCRIPT_DIR/06_variant_processing.sh" 2>&1 | tee -a "$LOG_DIR/06_processing.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Variant processing failed"
    fi
    
    # Step 7: Phasing
    if should_run_step "phasing"; then
        log "Starting Step 7: Variant Phasing"
        "$SCRIPT_DIR/07_phasing.sh" 2>&1 | tee -a "$LOG_DIR/07_phasing.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Phasing failed"
    fi
    
    # Step 8: Neoantigen Prediction
    if should_run_step "neoantigen_prediction"; then
        log "Starting Step 8: Neoantigen Prediction"
        "$SCRIPT_DIR/08_neoantigen_prediction.sh" 2>&1 | tee -a "$LOG_DIR/08_neoantigens.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Neoantigen prediction failed"
    fi
    
    # Step 9: Visualization
    if should_run_step "visualization"; then
        log "Starting Step 9: Interactive Visualization"
        "$SCRIPT_DIR/09_visualization.sh" 2>&1 | tee -a "$LOG_DIR/09_visualization.log"
        [[ ${PIPESTATUS[0]} -eq 0 ]] || error_exit "Visualization setup failed"
    fi
    
    # Generate summary report
    log "Generating summary report..."
    "$SCRIPT_DIR/generate_summary.sh" 2>&1 | tee -a "$LOG_DIR/summary.log"
    
    log "=========================================="
    log "PIPELINE COMPLETED SUCCESSFULLY"
    log "Summary report: $OUTPUT_DIR/summary_report.html"
    log "=========================================="
}

# Run main function
main