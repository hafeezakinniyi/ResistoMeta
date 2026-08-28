#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 07 - Contig Annotation
#
# Purpose:
#   Identifies antimicrobial resistance genes and plasmid replicons from
#   assembled metagenomic contigs using:
#
#       • ABRicate (CARD)
#       • ABRicate (PlasmidFinder)
#
#   ABRicate-CARD provides contig-level antimicrobial resistance gene
#   annotations, while ABRicate-PlasmidFinder identifies plasmid replicons.
#
#   These outputs are used by the downstream R analysis for:
#       • contig-associated ARG profiling
#       • plasmid replicon profiling
#       • ARG–plasmid co-occurrence
#       • contig-level ARG/plasmid co-localization
#
# Usage:
#     bash 07_contig_annotation.sh SAMPLE_ID
#
###############################################################################

set -euo pipefail


###############################################################################
# DETERMINE PROJECT ROOT
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(dirname "$(dirname "$(dirname "${SCRIPT_DIR}")")")"


###############################################################################
# LOAD CONFIGURATION
###############################################################################

source "${PROJECT_ROOT}/config/config.sh"

source "${SCRIPT_DIR_BASH}/utils.sh"


###############################################################################
# MODULE DIRECTORIES
###############################################################################

MEGAHIT_DIR="${ASSEMBLY_DIR}/megahit"

ABRICATE_CARD_DIR="${CONTIG_ANNOTATION_DIR}/abricate_card"

PLASMIDFINDER_DIR="${CONTIG_ANNOTATION_DIR}/plasmidfinder"


###############################################################################
# RUN CONTIG ANNOTATION
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 07 - Contig Annotation (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILE
    ###########################################################################

    local CONTIGS="${MEGAHIT_DIR}/${SAMPLE_ID}/final.contigs.fa"

    check_file_exists "${CONTIGS}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${ABRICATE_CARD_DIR}"

    create_directory "${PLASMIDFINDER_DIR}"


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local CARD_OUTPUT="${ABRICATE_CARD_DIR}/${SAMPLE_ID}_card.tsv"

    local PLASMID_OUTPUT="${PLASMIDFINDER_DIR}/${SAMPLE_ID}_plasmidfinder.tsv"


    ###########################################################################
    # SKIP COMPLETED SAMPLE
    #
    # Both annotation outputs must exist before the module is considered
    # complete.
    ###########################################################################

    if skip_step \
        "Contig annotation (${SAMPLE_ID})" \
        "${CARD_OUTPUT}" \
        "${PLASMID_OUTPUT}"
    then

        return 0

    fi


    ###########################################################################
    # ACTIVATE MAIN ENVIRONMENT
    ###########################################################################

    activate_env "${MAIN_ENV}"

    check_tool abricate


    ###########################################################################
    # ABRICATE - CARD
    #
    # Identifies antimicrobial resistance genes on assembled contigs.
    ###########################################################################

    log_message \
        "INFO" \
        "Running ABRicate (CARD)..."


        abricate \
        --db "${ABRICATE_CARD_DB}" \
        --minid "${ABRICATE_MIN_IDENTITY}" \
        --mincov "${ABRICATE_MIN_COVERAGE}" \
        "${CONTIGS}" \
        > "${CARD_OUTPUT}"


    ###########################################################################
    # VALIDATE CARD OUTPUT
    ###########################################################################

    check_file_exists "${CARD_OUTPUT}"


    ###########################################################################
    # ABRICATE - PLASMIDFINDER
    #
    # Identifies plasmid replicons on assembled contigs.
    ###########################################################################

    log_message \
        "INFO" \
        "Running ABRicate (PlasmidFinder)..."


        abricate \
        --db "${ABRICATE_PLASMID_DB}" \
        --minid "${ABRICATE_MIN_IDENTITY}" \
        --mincov "${ABRICATE_MIN_COVERAGE}" \
        "${CONTIGS}" \
        > "${PLASMID_OUTPUT}"


    ###########################################################################
    # VALIDATE PLASMIDFINDER OUTPUT
    ###########################################################################

    check_file_exists "${PLASMID_OUTPUT}"


    ###########################################################################
    # COMPLETION
    ###########################################################################

    log_message \
        "SUCCESS" \
        "Module 07 completed for sample ${SAMPLE_ID}."

}


###############################################################################
# MAIN
###############################################################################

if [[ $# -eq 0 ]]; then

    log_message \
        "ERROR" \
        "Usage: $0 SAMPLE_ID"

    exit 1

fi


run_sample "$1"