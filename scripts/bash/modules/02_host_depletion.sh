#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 02 - Host Depletion
#
# Purpose:
#   Removes host-derived reads from quality-controlled paired-end reads using
#   Bowtie2 and evaluates the resulting host-depleted reads using FastQC.
#
# Usage:
#     bash 02_host_depletion.sh SAMPLE_ID
#     bash 02_host_depletion.sh summary
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

TRIMMED_READS_DIR="${QC_DIR}/trimmed_reads"

HOST_DEPLETED_READS_DIR="${HOST_DEPLETION_DIR}/depleted_reads"

HOST_DEPLETED_FASTQC_DIR="${HOST_DEPLETION_DIR}/fastqc"

MULTIQC_HOST_DEPLETED_DIR="${HOST_DEPLETION_DIR}/multiqc"


###############################################################################
# RUN HOST DEPLETION FOR A SINGLE SAMPLE
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 02 - Host Depletion (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILES
    ###########################################################################

    local R1="${TRIMMED_READS_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"

    local R2="${TRIMMED_READS_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz"

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${HOST_DEPLETED_READS_DIR}"

    create_directory "${HOST_DEPLETED_FASTQC_DIR}"


    ###########################################################################
    # SOFTWARE
    ###########################################################################

    activate_env "${MAIN_ENV}"

    check_tool bowtie2

    check_tool fastqc

    check_tool multiqc


    ###########################################################################
    # HOST REFERENCE
    ###########################################################################

    check_bowtie2_index "${HOST_REFERENCE_INDEX}"


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local HOST_R1="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R1_depleted.fastq.gz"

    local HOST_R2="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R2_depleted.fastq.gz"


    local HOST_FASTQC_R1="${HOST_DEPLETED_FASTQC_DIR}/${SAMPLE_ID}_R1_depleted_fastqc.html"

    local HOST_FASTQC_R2="${HOST_DEPLETED_FASTQC_DIR}/${SAMPLE_ID}_R2_depleted_fastqc.html"


    ###########################################################################
    # SKIP COMPLETED SAMPLE
    ###########################################################################

    if skip_step \
        "Host depletion (${SAMPLE_ID})" \
        "${HOST_R1}" \
        "${HOST_R2}" \
        "${HOST_FASTQC_R1}" \
        "${HOST_FASTQC_R2}"
    then

        return 0

    fi


    ###########################################################################
    # REMOVE HOST-DERIVED READS
    ###########################################################################

    log_message "INFO" \
        "Removing host-derived reads using Bowtie2..."

    run_command \
        bowtie2 \
        --threads "${THREADS}" \
        -x "${HOST_REFERENCE_INDEX}" \
        -1 "${R1}" \
        -2 "${R2}" \
        --un-conc-gz "${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R%_depleted.fastq.gz" \
        -S /dev/null


    ###########################################################################
    # VALIDATE HOST-DEPLETED READS
    ###########################################################################

    check_output "${HOST_R1}"

    check_output "${HOST_R2}"


    ###########################################################################
    # FASTQC
    ###########################################################################

    log_message "INFO" \
        "Running FastQC on host-depleted reads..."

    run_command \
        fastqc \
        --threads "${THREADS}" \
        --outdir "${HOST_DEPLETED_FASTQC_DIR}" \
        "${HOST_R1}" \
        "${HOST_R2}"


    ###########################################################################
    # VALIDATE FASTQC OUTPUTS
    ###########################################################################

    check_output "${HOST_FASTQC_R1}"

    check_output "${HOST_FASTQC_R2}"


    log_message \
        "SUCCESS" \
        "Module 02 completed for sample ${SAMPLE_ID}."

}


###############################################################################
# GENERATE PROJECT-LEVEL MULTIQC REPORT
###############################################################################

run_summary() {

    print_step_header \
        "Module 02 - MultiQC Summary"


    create_directory "${MULTIQC_HOST_DEPLETED_DIR}"


    activate_env "${MAIN_ENV}"

    check_tool multiqc


    ###########################################################################
    # MULTIQC
    ###########################################################################

    log_message "INFO" \
        "Generating MultiQC report for host-depleted reads..."

    run_command \
        multiqc \
        "${HOST_DEPLETED_FASTQC_DIR}" \
        --outdir "${MULTIQC_HOST_DEPLETED_DIR}" \
        --force


    ###########################################################################
    # VALIDATE MULTIQC OUTPUT
    ###########################################################################

    check_output \
        "${MULTIQC_HOST_DEPLETED_DIR}/multiqc_report.html"


    log_message \
        "SUCCESS" \
        "Module 02 MultiQC summary generated."

}


###############################################################################
# MAIN
###############################################################################

if [[ $# -eq 0 ]]; then

    log_message \
        "ERROR" \
        "Usage: $0 SAMPLE_ID | summary"

    exit 1

fi


case "$1" in

    summary)

        run_summary

        ;;

    *)

        run_sample "$1"

        ;;

esac