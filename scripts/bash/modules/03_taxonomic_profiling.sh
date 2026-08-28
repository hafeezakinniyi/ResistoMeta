#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 03 - Taxonomic Profiling
#
# Purpose:
#   Performs taxonomic classification of host-depleted reads using Kraken2,
#   estimates taxonomic abundance using Bracken, and produces sample-level
#   taxonomic profiles for downstream analyses.
#
# Usage:
#     bash 03_taxonomic_profiling.sh SAMPLE_ID
###############################################################################

set -euo pipefail


###############################################################################
# DETERMINE PROJECT ROOT
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(dirname "$(dirname "$(dirname "${SCRIPT_DIR}")")")"


###############################################################################
# LOAD CONFIGURATION AND UTILITIES
###############################################################################

source "${PROJECT_ROOT}/config/config.sh"

source "${SCRIPT_DIR_BASH}/utils.sh"


###############################################################################
# MODULE DIRECTORIES
###############################################################################

HOST_DEPLETED_READS_DIR="${HOST_DEPLETION_DIR}/depleted_reads"

KRAKEN_DIR="${TAXONOMY_DIR}/kraken2"

BRACKEN_DIR="${TAXONOMY_DIR}/bracken"


###############################################################################
# RUN TAXONOMIC PROFILING
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 03 - Taxonomic Profiling (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILES
    ###########################################################################

    local R1="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R1_depleted.fastq.gz"

    local R2="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R2_depleted.fastq.gz"

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${KRAKEN_DIR}"

    create_directory "${BRACKEN_DIR}"


    ###########################################################################
    # SOFTWARE
    ###########################################################################

    activate_env "${MAIN_ENV}"

    check_tool kraken2

    check_tool bracken


    ###########################################################################
    # DATABASE
    ###########################################################################

    check_directory_exists "${KRAKEN_DATABASE}"


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local KRAKEN_REPORT="${KRAKEN_DIR}/${SAMPLE_ID}.report"

    local KRAKEN_OUTPUT="${KRAKEN_DIR}/${SAMPLE_ID}.kraken"

    local BRACKEN_SPECIES="${BRACKEN_DIR}/${SAMPLE_ID}_species.bracken"

    local BRACKEN_GENUS="${BRACKEN_DIR}/${SAMPLE_ID}_genus.bracken"

    local BRACKEN_PHYLUM="${BRACKEN_DIR}/${SAMPLE_ID}_phylum.bracken"


    ###########################################################################
    # SKIP COMPLETED SAMPLE
    ###########################################################################

    if skip_step \
        "Taxonomic profiling (${SAMPLE_ID})" \
        "${KRAKEN_REPORT}" \
        "${KRAKEN_OUTPUT}" \
        "${BRACKEN_SPECIES}" \
        "${BRACKEN_GENUS}" \
        "${BRACKEN_PHYLUM}"
    then

        return 0

    fi


    ###########################################################################
    # KRAKEN2 CLASSIFICATION
    ###########################################################################

    log_message "INFO" \
        "Running Kraken2..."

    run_command \
        kraken2 \
        --db "${KRAKEN_DATABASE}" \
        --paired \
        "${R1}" \
        "${R2}" \
        --threads "${THREADS}" \
        --report "${KRAKEN_REPORT}" \
        --output "${KRAKEN_OUTPUT}"


    ###########################################################################
    # VALIDATE KRAKEN2 OUTPUTS
    ###########################################################################

    check_output "${KRAKEN_REPORT}"

    check_output "${KRAKEN_OUTPUT}"


    ###########################################################################
    # BRACKEN SPECIES ABUNDANCE
    ###########################################################################

    log_message "INFO" \
        "Estimating species abundance using Bracken..."

    run_command \
        bracken \
        -d "${KRAKEN_DATABASE}" \
        -i "${KRAKEN_REPORT}" \
        -o "${BRACKEN_SPECIES}" \
        -r "${BRACKEN_READ_LENGTH}" \
        -l S


    ###########################################################################
    # BRACKEN GENUS ABUNDANCE
    ###########################################################################

    log_message "INFO" \
        "Estimating genus abundance using Bracken..."

    run_command \
        bracken \
        -d "${KRAKEN_DATABASE}" \
        -i "${KRAKEN_REPORT}" \
        -o "${BRACKEN_GENUS}" \
        -r "${BRACKEN_READ_LENGTH}" \
        -l G


    ###########################################################################
    # BRACKEN PHYLUM ABUNDANCE
    ###########################################################################

    log_message "INFO" \
        "Estimating phylum abundance using Bracken..."

    run_command \
        bracken \
        -d "${KRAKEN_DATABASE}" \
        -i "${KRAKEN_REPORT}" \
        -o "${BRACKEN_PHYLUM}" \
        -r "${BRACKEN_READ_LENGTH}" \
        -l P


    ###########################################################################
    # VALIDATE BRACKEN OUTPUTS
    ###########################################################################

    check_output "${BRACKEN_SPECIES}"

    check_output "${BRACKEN_GENUS}"

    check_output "${BRACKEN_PHYLUM}"


    log_message \
        "SUCCESS" \
        "Module 03 completed for sample ${SAMPLE_ID}."

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