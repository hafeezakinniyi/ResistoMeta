#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 05 - Metagenome Assembly
#
# Purpose:
#   Performs de novo metagenome assembly of host-depleted paired-end reads
#   using MEGAHIT.
#
# Usage:
#     bash 05_metagenome_assembly.sh SAMPLE_ID
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

HOST_DEPLETED_READS_DIR="${HOST_DEPLETION_DIR}/depleted_reads"

MEGAHIT_DIR="${ASSEMBLY_DIR}/megahit"


###############################################################################
# RUN METAGENOME ASSEMBLY
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 05 - Metagenome Assembly (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILES
    ###########################################################################

    local R1="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R1_depleted.fastq.gz"

    local R2="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R2_depleted.fastq.gz"

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${MEGAHIT_DIR}"


    ###########################################################################
    # SOFTWARE
    ###########################################################################

    activate_env "${MAIN_ENV}"

    check_tool megahit


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local SAMPLE_OUTPUT_DIR="${MEGAHIT_DIR}/${SAMPLE_ID}"

    local CONTIGS="${SAMPLE_OUTPUT_DIR}/final.contigs.fa"


    ###########################################################################
    # SKIP COMPLETED SAMPLE
    ###########################################################################

    if skip_step \
        "Metagenomic assembly (${SAMPLE_ID})" \
        "${CONTIGS}"
    then

        return 0

    fi


    ###########################################################################
    # MEGAHIT
    ###########################################################################

    log_message "INFO" \
        "Running MEGAHIT..."

    run_command \
        megahit \
        -1 "${R1}" \
        -2 "${R2}" \
        -o "${SAMPLE_OUTPUT_DIR}" \
        -t "${THREADS}"


    ###########################################################################
    # VALIDATE OUTPUT
    ###########################################################################

    check_output "${CONTIGS}"


    ###########################################################################
    # COMPLETION
    ###########################################################################

    log_message \
        "SUCCESS" \
        "Module 05 completed for sample ${SAMPLE_ID}."

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