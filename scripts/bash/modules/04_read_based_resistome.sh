#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 04 - Read-based Resistome
#
# Purpose:
#   Detects antimicrobial resistance genes directly from host-depleted reads
#   using the CARD Resistance Gene Identifier (RGI) BWT workflow.
#
#   This module performs:
#       • read-level ARG detection using RGI BWT
#       • sequencing-depth calculation from the actual RGI input reads
#       • R1/R2 paired-read equality validation
#       • project-level sequencing-depth summary generation
#
#   Sequencing depth definition:
#
#       For paired-end sequencing:
#
#           R1 reads = R2 reads = total_sequences
#
#       total_sequences therefore represents the number of paired-end
#       sequences (read pairs), rather than the sum of R1 and R2 reads.
#
#   Normalization, richness, prevalence, and drug-class analyses are performed
#   downstream by resistome_analysis.R.
#
# Usage:
#     bash 04_read_based_resistome.sh SAMPLE_ID
#     bash 04_read_based_resistome.sh summary
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

RGI_READS_DIR="${READ_RESISTOME_DIR}/rgi"


###############################################################################
# COUNT READS IN A GZIPPED FASTQ FILE
#
# FASTQ format contains four lines per sequencing read:
#
#     @header
#     sequence
#     +
#     quality
#
# Therefore:
#
#     number of reads = number of lines / 4
#
# The function also validates that the FASTQ contains a complete number of
# four-line records.
###############################################################################

count_fastq_reads() {

    local FASTQ="$1"

    local LINE_COUNT
    local READ_COUNT


    ###########################################################################
    # VALIDATE INPUT
    ###########################################################################

    check_file_exists "${FASTQ}"


    ###########################################################################
    # COUNT FASTQ LINES
    ###########################################################################

    LINE_COUNT=$(gzip -cd "${FASTQ}" | wc -l)


    ###########################################################################
    # VALIDATE FASTQ STRUCTURE
    ###########################################################################

    if (( LINE_COUNT % 4 != 0 )); then

        log_message \
            "ERROR" \
            "Invalid FASTQ structure detected in ${FASTQ}. " \
            "Line count (${LINE_COUNT}) is not divisible by 4."

        return 1

    fi


    ###########################################################################
    # CALCULATE NUMBER OF READS
    ###########################################################################

    READ_COUNT=$((LINE_COUNT / 4))


    if (( READ_COUNT <= 0 )); then

        log_message \
            "ERROR" \
            "No sequencing reads detected in ${FASTQ}."

        return 1

    fi


    echo "${READ_COUNT}"

}


###############################################################################
# CALCULATE SAMPLE SEQUENCING DEPTH
#
# Sequencing depth is defined here as the number of paired-end sequences
# present in the host-depleted FASTQ files supplied directly to RGI BWT.
#
# For valid paired-end data:
#
#     R1 reads = R2 reads
#
# and:
#
#     total_sequences = R1 reads = R2 reads
#
# IMPORTANT:
#
#     total_sequences is NOT R1 + R2.
#
# For example:
#
#     R1 = 25,000,000 reads
#     R2 = 25,000,000 reads
#
#     total_sequences = 25,000,000 read pairs
#
###############################################################################

calculate_sample_depth() {

    local SAMPLE_ID="$1"

    local R1="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R1_depleted.fastq.gz"

    local R2="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R2_depleted.fastq.gz"

    local R1_READS
    local R2_READS
    local TOTAL_SEQUENCES


    ###########################################################################
    # VALIDATE INPUT READ FILES
    ###########################################################################

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # COUNT R1 AND R2
    ###########################################################################

    log_message \
        "INFO" \
        "Counting sequencing reads for ${SAMPLE_ID}..."


    R1_READS=$(count_fastq_reads "${R1}")

    R2_READS=$(count_fastq_reads "${R2}")


    ###########################################################################
    # VALIDATE R1/R2 READ COUNTS
    ###########################################################################

    if [[ "${R1_READS}" -ne "${R2_READS}" ]]; then

        log_message \
            "ERROR" \
            "R1/R2 read-count mismatch for ${SAMPLE_ID}: " \
            "R1=${R1_READS}, R2=${R2_READS}."

        return 1

    fi


    ###########################################################################
    # DEFINE TOTAL SEQUENCING DEPTH
    #
    # For paired-end sequencing, one R1 read and one R2 read constitute one
    # paired-end sequence.
    #
    # Therefore total_sequences equals the number of read pairs.
    ###########################################################################

    TOTAL_SEQUENCES="${R1_READS}"


    ###########################################################################
    # VALIDATE TOTAL
    ###########################################################################

    if [[ "${TOTAL_SEQUENCES}" -le 0 ]]; then

        log_message \
            "ERROR" \
            "No sequencing reads detected for ${SAMPLE_ID}."

        return 1

    fi


    ###########################################################################
    # REPORT
    ###########################################################################

    log_message \
        "INFO" \
        "${SAMPLE_ID}: R1_reads=${R1_READS}, " \
        "R2_reads=${R2_READS}, " \
        "total_sequences=${TOTAL_SEQUENCES} " \
        "(paired-end read pairs)"


    ###########################################################################
    # RETURN TAB-SEPARATED RESULT
    #
    # Output:
    #
    #     SAMPLE_ID    R1_READS    R2_READS    TOTAL_SEQUENCES
    ###########################################################################

    printf "%s\t%s\t%s\t%s\n" \
        "${SAMPLE_ID}" \
        "${R1_READS}" \
        "${R2_READS}" \
        "${TOTAL_SEQUENCES}"

}


###############################################################################
# GENERATE PROJECT-LEVEL SEQUENCING DEPTH TABLE
#
# This function scans the host-depleted R1 files to identify samples and
# calculates sequencing depth from the actual FASTQ files used as RGI inputs.
#
# Output:
#
#     ${SEQUENCING_DEPTH_FILE}
#
# Columns:
#
#     SampleID
#     R1_reads
#     R2_reads
#     total_sequences
#
# Definition:
#
#     total_sequences = R1_reads = R2_reads
#
###############################################################################

generate_sequencing_depth() {

    print_step_header \
        "Module 04 - Generate Sequencing Depth Summary"


    ###########################################################################
    # VALIDATE INPUT DIRECTORY
    ###########################################################################

    check_directory_exists "${HOST_DEPLETED_READS_DIR}"


    ###########################################################################
    # CREATE OUTPUT DIRECTORY
    #
    # SEQUENCING_DEPTH_FILE is located inside READ_RESISTOME_DIR according
    # to config.sh.
    ###########################################################################

    create_directory "${READ_RESISTOME_DIR}"


    ###########################################################################
    # TEMPORARY FILE
    #
    # The temporary file is written first so that a failed calculation does
    # not leave behind a partially written sequencing-depth.csv.
    ###########################################################################

    local TEMP_FILE="${SEQUENCING_DEPTH_FILE}.tmp"

    rm -f "${TEMP_FILE}"


    ###########################################################################
    # WRITE HEADER
    ###########################################################################

    printf "SampleID,R1_reads,R2_reads,total_sequences\n" \
        > "${TEMP_FILE}"


    ###########################################################################
    # IDENTIFY HOST-DEPLETED R1 FILES
    ###########################################################################

    local R1_FILE
    local SAMPLE_ID
    local RESULT
    local SAMPLE_COUNT=0


    shopt -s nullglob

    local R1_FILES=(
        "${HOST_DEPLETED_READS_DIR}"/*_R1_depleted.fastq.gz
    )

    shopt -u nullglob


    if [[ "${#R1_FILES[@]}" -eq 0 ]]; then

        rm -f "${TEMP_FILE}"

        log_message \
            "ERROR" \
            "No host-depleted R1 FASTQ files found in " \
            "${HOST_DEPLETED_READS_DIR}."

        return 1

    fi


    ###########################################################################
    # PROCESS EACH SAMPLE
    ###########################################################################

    for R1_FILE in "${R1_FILES[@]}"
    do

        SAMPLE_ID="$(basename "${R1_FILE}")"

        SAMPLE_ID="${SAMPLE_ID%_R1_depleted.fastq.gz}"


        log_message \
            "INFO" \
            "Generating sequencing depth for ${SAMPLE_ID}..."


        RESULT="$(calculate_sample_depth "${SAMPLE_ID}")"


        printf "%s\n" "${RESULT}" | \
            awk -F '\t' \
                'NF == 4 {
                    print $1 "," $2 "," $3 "," $4
                }' \
            >> "${TEMP_FILE}"


        SAMPLE_COUNT=$((SAMPLE_COUNT + 1))

    done


    ###########################################################################
    # SORT OUTPUT BY SAMPLE ID
    #
    # Keep the header as the first line.
    ###########################################################################

    {
        head -n 1 "${TEMP_FILE}"

        tail -n +2 "${TEMP_FILE}" | sort -t ',' -k1,1

    } > "${SEQUENCING_DEPTH_FILE}"


    ###########################################################################
    # REMOVE TEMPORARY FILE
    ###########################################################################

    rm -f "${TEMP_FILE}"


    ###########################################################################
    # VALIDATE OUTPUT
    ###########################################################################

    check_output "${SEQUENCING_DEPTH_FILE}"


    ###########################################################################
    # REPORT
    ###########################################################################

    log_message \
        "SUCCESS" \
        "Sequencing-depth summary generated for ${SAMPLE_COUNT} samples."

    log_message \
        "INFO" \
        "Sequencing-depth definition: R1 reads = R2 reads = total_sequences."

    log_message \
        "INFO" \
        "Output: ${SEQUENCING_DEPTH_FILE}"

}


###############################################################################
# RUN READ-BASED RESISTOME ANALYSIS
###############################################################################

run_sample() {

    local SAMPLE_ID="$1"

    print_step_header \
        "Module 04 - Read-based Resistome (${SAMPLE_ID})"


    ###########################################################################
    # INPUT FILES
    ###########################################################################

    local R1="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R1_depleted.fastq.gz"

    local R2="${HOST_DEPLETED_READS_DIR}/${SAMPLE_ID}_R2_depleted.fastq.gz"

    check_paired_reads "${R1}" "${R2}"


    ###########################################################################
    # OUTPUT DIRECTORIES
    ###########################################################################

    create_directory "${RGI_READS_DIR}"


    ###########################################################################
    # SOFTWARE
    ###########################################################################

    activate_env "${RGI_ENV}"

    check_tool rgi


    ###########################################################################
    # OUTPUT FILES
    ###########################################################################

    local OUTPUT_PREFIX="${RGI_READS_DIR}/${SAMPLE_ID}"

    local GENE_MAPPING_OUTPUT="${OUTPUT_PREFIX}.gene_mapping_data.txt"


    ###########################################################################
    # SKIP COMPLETED SAMPLE
    ###########################################################################

    if skip_step \
        "Read-based resistome analysis (${SAMPLE_ID})" \
        "${GENE_MAPPING_OUTPUT}"
    then

        return 0

    fi


    ###########################################################################
    # RGI BWT
    ###########################################################################

    log_message \
        "INFO" \
        "Running RGI BWT..."


    run_command \
        rgi bwt \
        --read_one "${R1}" \
        --read_two "${R2}" \
        --output_file "${OUTPUT_PREFIX}" \
        --local \
        -n "${THREADS}"


    ###########################################################################
    # VALIDATE OUTPUT
    ###########################################################################

    check_output "${GENE_MAPPING_OUTPUT}"


    ###########################################################################
    # COMPLETION
    ###########################################################################

    log_message \
        "SUCCESS" \
        "Module 04 completed for sample ${SAMPLE_ID}."

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


###############################################################################
# SUMMARY MODE
###############################################################################

if [[ "$1" == "summary" ]]; then

    generate_sequencing_depth

    exit 0

fi


###############################################################################
# SAMPLE MODE
###############################################################################

run_sample "$1"