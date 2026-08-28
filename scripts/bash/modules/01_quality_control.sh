#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 01 - Quality Control
#
# Purpose:
#   Performs quality assessment of paired-end sequencing reads by:
#     1. Running FastQC on raw reads
#     2. Trimming adapters and low-quality bases using fastp
#     3. Running FastQC on trimmed reads
#     4. Generating project-level MultiQC reports
#
# Usage:
#     bash 01_quality_control.sh SAMPLE_ID
#     bash 01_quality_control.sh summary
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
# MODULE OUTPUT DIRECTORIES
###############################################################################

RAW_FASTQC_DIR="${QC_DIR}/raw_fastqc"

TRIMMED_FASTQC_DIR="${QC_DIR}/trimmed_fastqc"

FASTP_DIR="${QC_DIR}/fastp"

TRIMMED_READS_DIR="${QC_DIR}/trimmed_reads"

MULTIQC_DIR="${QC_DIR}/multiqc"


###############################################################################
# RUN QUALITY CONTROL FOR A SINGLE SAMPLE
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 01 - Quality Control (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILES
    ###########################################################################

    local R1="${RAW_READS_DIR}/${SAMPLE_ID}_R1_001.fastq.gz"

    local R2="${RAW_READS_DIR}/${SAMPLE_ID}_R2_001.fastq.gz"

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${RAW_FASTQC_DIR}"

    create_directory "${TRIMMED_FASTQC_DIR}"

    create_directory "${FASTP_DIR}"

    create_directory "${TRIMMED_READS_DIR}"


    ###########################################################################
    # SOFTWARE
    ###########################################################################

    activate_env "${MAIN_ENV}"

    check_tool fastqc

    check_tool fastp

    check_tool multiqc


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local TRIMMED_R1="${TRIMMED_READS_DIR}/${SAMPLE_ID}_R1_trimmed.fastq.gz"

    local TRIMMED_R2="${TRIMMED_READS_DIR}/${SAMPLE_ID}_R2_trimmed.fastq.gz"

    local FASTP_HTML="${FASTP_DIR}/${SAMPLE_ID}.html"

    local FASTP_JSON="${FASTP_DIR}/${SAMPLE_ID}.json"


    local RAW_FASTQC_R1="${RAW_FASTQC_DIR}/${SAMPLE_ID}_R1_001_fastqc.html"

    local RAW_FASTQC_R2="${RAW_FASTQC_DIR}/${SAMPLE_ID}_R2_001_fastqc.html"


    local TRIMMED_FASTQC_R1="${TRIMMED_FASTQC_DIR}/${SAMPLE_ID}_R1_trimmed_fastqc.html"

    local TRIMMED_FASTQC_R2="${TRIMMED_FASTQC_DIR}/${SAMPLE_ID}_R2_trimmed_fastqc.html"


    ###########################################################################
    # CHECK WHETHER MODULE IS ALREADY COMPLETE
    ###########################################################################

    if skip_step \
        "Quality Control (${SAMPLE_ID})" \
        "${TRIMMED_R1}" \
        "${TRIMMED_R2}" \
        "${FASTP_HTML}" \
        "${FASTP_JSON}" \
        "${RAW_FASTQC_R1}" \
        "${RAW_FASTQC_R2}" \
        "${TRIMMED_FASTQC_R1}" \
        "${TRIMMED_FASTQC_R2}"
    then

        return 0

    fi


    ###########################################################################
    # FASTQC - RAW READS
    ###########################################################################

    log_message "INFO" \
        "Running FastQC on raw reads..."

    run_command \
        fastqc \
        "${R1}" \
        "${R2}" \
        --threads "${THREADS}" \
        --outdir "${RAW_FASTQC_DIR}"


    ###########################################################################
    # FASTP
    ###########################################################################

    log_message "INFO" \
        "Running fastp..."

    run_command \
        fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${TRIMMED_R1}" \
        --out2 "${TRIMMED_R2}" \
        --detect_adapter_for_pe \
        --cut_tail \
        --cut_tail_mean_quality "${MIN_QUALITY}" \
        --length_required "${MIN_READ_LENGTH}" \
        --thread "${THREADS}" \
        --html "${FASTP_HTML}" \
        --json "${FASTP_JSON}"


    ###########################################################################
    # VALIDATE FASTP OUTPUTS
    ###########################################################################

    check_output "${TRIMMED_R1}"

    check_output "${TRIMMED_R2}"

    check_output "${FASTP_HTML}"

    check_output "${FASTP_JSON}"


    ###########################################################################
    # FASTQC - TRIMMED READS
    ###########################################################################

    log_message "INFO" \
        "Running FastQC on trimmed reads..."

    run_command \
        fastqc \
        "${TRIMMED_R1}" \
        "${TRIMMED_R2}" \
        --threads "${THREADS}" \
        --outdir "${TRIMMED_FASTQC_DIR}"


    ###########################################################################
    # VALIDATE FASTQC OUTPUTS
    ###########################################################################

    check_output "${RAW_FASTQC_R1}"

    check_output "${RAW_FASTQC_R2}"

    check_output "${TRIMMED_FASTQC_R1}"

    check_output "${TRIMMED_FASTQC_R2}"


    log_message \
        "SUCCESS" \
        "Module 01 completed for sample ${SAMPLE_ID}."

}


###############################################################################
# GENERATE PROJECT-LEVEL MULTIQC REPORTS
###############################################################################

run_summary() {

    print_step_header \
        "Module 01 - MultiQC Summary"


    local RAW_MULTIQC_DIR="${MULTIQC_DIR}/raw"

    local TRIMMED_MULTIQC_DIR="${MULTIQC_DIR}/trimmed"


    create_directory "${RAW_MULTIQC_DIR}"

    create_directory "${TRIMMED_MULTIQC_DIR}"


    activate_env "${MAIN_ENV}"

    check_tool multiqc


    ###########################################################################
    # RAW FASTQC SUMMARY
    ###########################################################################

    log_message "INFO" \
        "Generating MultiQC report for raw reads..."

    run_command \
        multiqc \
        "${RAW_FASTQC_DIR}" \
        --outdir "${RAW_MULTIQC_DIR}" \
        --force


    check_output \
        "${RAW_MULTIQC_DIR}/multiqc_report.html"


    ###########################################################################
    # TRIMMED FASTQC + FASTP SUMMARY
    ###########################################################################

    log_message "INFO" \
        "Generating MultiQC report for trimmed reads..."

    run_command \
        multiqc \
        "${TRIMMED_FASTQC_DIR}" \
        "${FASTP_DIR}" \
        --outdir "${TRIMMED_MULTIQC_DIR}" \
        --force


    check_output \
        "${TRIMMED_MULTIQC_DIR}/multiqc_report.html"


    log_message \
        "SUCCESS" \
        "MultiQC summaries generated."

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