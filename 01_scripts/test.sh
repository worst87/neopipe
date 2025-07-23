#!/bin/bash
# =================================================================================
# DIAGNOSTIC SCRIPT: Find Complex Somatic + Proximal Frameshift Events
#
# This script identifies somatic variants that have a phased frameshift variant
# nearby, which are the cases pVACseq cannot handle correctly.
# =================================================================================

set -e

# --- Initialize Conda ---
echo "--- Activating Conda environment: variant_env ---"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate variant_env

# --- Source Config ---
source "$(dirname "$0")/config.sh"

echo "======================================================"
echo "IDENTIFYING SOMATIC VARIANTS WITH PROXIMAL FRAMESHIFTS"
echo "======================================================"

for TUMOR_SAMPLE in "${!TUMOR_NORMAL_PAIRS[@]}"; do
    echo -e "\n--- Analyzing VCF for: $TUMOR_SAMPLE ---"

    # --- Define Paths ---
    VARIANT_DIR="$PROJ_DIR/03_pipeline_output/04_variants_vcf/$TUMOR_SAMPLE"
    PHASED_VCF="$VARIANT_DIR/phased/${TUMOR_SAMPLE}.phased.annotated.vcf.gz"
    DIAGNOSTIC_DIR="$VARIANT_DIR/phased/diagnostics"
    mkdir -p "$DIAGNOSTIC_DIR"
    
    OUTPUT_REPORT="$DIAGNOSTIC_DIR/complex_frameshift_events.tsv"

    if [[ ! -f "$PHASED_VCF" ]]; then
        echo "ERROR: Phased VCF not found: $PHASED_VCF"
        continue
    fi

    echo "Parsing $PHASED_VCF to find complex events. This may take a minute..."

    # Use bcftools and awk to find these events.
    # 1. Query for relevant fields: Chrom, Pos, Ref, Alt, VCB, CSQ, and the PS (Phase Set) from the FORMAT field.
    # 2. Filter for lines that are actually phased (PS is not '.').
    # 3. Use awk to process the phased variants:
    #    - Store variants in an array keyed by their Phase Set ID.
    #    - When a new variant with the same Phase Set ID is found, iterate through the stored variants.
    #    - If one variant is somatic (from VCB) and the other is a frameshift (from CSQ), print a report.
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/VariantCalledBy\t%INFO/CSQ\t[%PS]\n' "$PHASED_VCF" | \
    grep -v "\[\.$" | \
    awk -F'\t' '
    BEGIN { OFS="\t"; print "PhaseSetID", "SomaticVar_Chrom", "SomaticVar_Pos", "SomaticVar_Ref", "SomaticVar_Alt", "Somatic_Consequence", "ProximalVar_Chrom", "ProximalVar_Pos", "ProximalVar_Ref", "ProximalVar_Alt", "Proximal_Consequence" }
    {
        ps_id = $1 ":" $7; # Create a unique phase set ID (chr:ps_value)
        
        # Extract the primary consequence
        split($6, csq_parts, "|");
        consequence = csq_parts[2];
        
        # Store current variant info
        line_info = $1 FS $2 FS $3 FS $4 FS consequence;
        is_somatic = ($5 != ".");
        is_frameshift = (consequence ~ /frameshift_variant/);

        # Check against already stored variants in the same phase set
        if (ps_id in phase_sets) {
            split(phase_sets[ps_id], stored_variants, ";");
            for (i in stored_variants) {
                split(stored_variants[i], fields, FS);
                stored_consequence = fields[5];
                
                stored_is_somatic = (fields[6] == "1");
                stored_is_frameshift = (stored_consequence ~ /frameshift_variant/);

                # Case 1: Current is somatic, stored is frameshift
                if (is_somatic && !is_frameshift && stored_is_frameshift) {
                    print ps_id, $1, $2, $3, $4, consequence, fields[1], fields[2], fields[3], fields[4], stored_consequence;
                }
                # Case 2: Current is frameshift, stored is somatic
                if (is_frameshift && stored_is_somatic) {
                    print ps_id, fields[1], fields[2], fields[3], fields[4], stored_consequence, $1, $2, $3, $4, consequence;
                }
            }
        }
        
        # Add current variant to the phase set list for future comparisons
        somatic_flag = (is_somatic ? "1" : "0");
        phase_sets[ps_id] = phase_sets[ps_id] (phase_sets[ps_id] ? ";" : "") line_info FS somatic_flag;
    }' > "$OUTPUT_REPORT"

    # --- Report Findings ---
    EVENT_COUNT=$( (tail -n +2 "$OUTPUT_REPORT" | wc -l) || echo 0)
    
    echo -e "\n--- COMPLEX EVENT SUMMARY for $TUMOR_SAMPLE ---"
    if [[ "$EVENT_COUNT" -gt 0 ]]; then
        echo "Found $EVENT_COUNT complex events (somatic variant + phased proximal frameshift)."
        echo "These events require special handling as pVACseq will misinterpret them."
        echo "A detailed report has been saved to: $OUTPUT_REPORT"
        echo -e "\nFirst 5 events found:"
        head -n 6 "$OUTPUT_REPORT" | column -t -s $'\t'
    else
        echo "No complex somatic + proximal frameshift events were found. The warnings you see are likely for other non-missense types (intronic, UTR, etc.) which are less critical to model."
    fi
    echo "--------------------------------------------------"

done

conda deactivate
echo -e "\n======================================================"
echo "DIAGNOSTICS COMPLETE"
echo "======================================================"
