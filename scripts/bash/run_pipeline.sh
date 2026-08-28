#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
#
# Master Pipeline Launcher
#
# Purpose:
#   Single entry point for the complete metagenomic AMR workflow.
#
# Architecture:
#
#   run_pipeline.sh
#       |
#       +-- config/config.sh
#       |
#       +-- Phase I: Bash bioinformatics modules
#       |      +-- 01_quality_control.sh
#       |      +-- 02_host_depletion.sh
#       |      +-- 03_taxonomic_profiling.sh
#       |      +-- 04_read_based_resistome.sh
#       |      +-- 05_metagenome_assembly.sh
#       |      +-- 06_assembly_quality.sh
#       |      +-- 07_contig_annotation.sh
#       |
#       +-- Phase II: R downstream analysis
#              +-- 08_launch_R_analysis.sh
#                       |
#                       +-- run_analysis.R
#                              +-- taxonomy_analysis.R
#                              +-- resistome_analysis.R
#                              +-- plasmid_analysis.R
#
# Execution modes:
#
#   full
#       Phase I + Phase II
#
#   bioinformatics
#       Phase I only
#
#   analysis
#       Phase II only
#
###############################################################################

set -euo pipefail


###############################################################################
# 1. DETERMINE SCRIPT AND PROJECT LOCATIONS
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(dirname "$(dirname "${SCRIPT_DIR}")")"


###############################################################################
# 2. LOAD CONFIGURATION AND UTILITY FUNCTIONS
###############################################################################

CONFIG_FILE="${PROJECT_ROOT}/config/config.sh"

if [[ ! -f "${CONFIG_FILE}" ]]; then

    echo "[ERROR] Configuration file not found:"
    echo "        ${CONFIG_FILE}"
    exit 1

fi


source "${CONFIG_FILE}"


if [[ -z "${SCRIPT_DIR_BASH:-}" ]]; then

    echo "[ERROR] SCRIPT_DIR_BASH is not defined by config.sh."
    exit 1

fi


if [[ ! -f "${SCRIPT_DIR_BASH}/utils.sh" ]]; then

    echo "[ERROR] Utility script not found:"
    echo "        ${SCRIPT_DIR_BASH}/utils.sh"
    exit 1

fi


source "${SCRIPT_DIR_BASH}/utils.sh"


###############################################################################
# 3. INITIALIZE PIPELINE
###############################################################################

PIPELINE_START_TIME="$(start_timer)"

initialize_logger

pipeline_banner

initialize_directories


log_message "INFO" "Pipeline version : ${PIPELINE_VERSION}"
log_message "INFO" "Execution mode   : ${PIPELINE_MODE}"
log_message "INFO" "Project root     : ${PROJECT_ROOT}"


###############################################################################
# 4. BASIC CONFIGURATION VALIDATION
###############################################################################

print_step_header "Pipeline validation"


bash "${SCRIPT_DIR_BASH}/validate_installation.sh"


print_step_footer "Pipeline validation completed"


###############################################################################
# 5. SAMPLE DISCOVERY
###############################################################################
#
# Sample discovery is required for Phase I because the Bash modules operate
# on individual samples.
#
# Phase II does not require discovery here because the R analysis launcher
# operates on the upstream result directories configured through config.sh.
#
###############################################################################

SAMPLE_IDS=()


if [[ "${PIPELINE_MODE}" != "analysis" ]]; then

    print_step_header "Sample discovery"


    if [[ ! -d "${RAW_READS_DIR}" ]]; then

        log_message \
            "ERROR" \
            "Raw reads directory not found: ${RAW_READS_DIR}"

        exit 1

    fi


    mapfile -t SAMPLE_IDS < <(

        find "${RAW_READS_DIR}" \
            -type f \
            -name "*_R1*.fastq.gz" \
            -exec basename {} \; |

        sed 's/_R1.*//' |

        sort -u

    )


    if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then

        log_message \
            "ERROR" \
            "No paired-end FASTQ files were found in ${RAW_READS_DIR}."

        exit 1

    fi


    log_message \
        "INFO" \
        "Discovered ${#SAMPLE_IDS[@]} sequencing sample(s)."


    printf '%s\n' "${SAMPLE_IDS[@]}" | while read -r SAMPLE_ID
    do
        log_message "INFO" "  Sample: ${SAMPLE_ID}"
    done


    print_step_footer "Sample discovery completed"

fi


###############################################################################
# 6. PHASE I — BIOINFORMATICS
###############################################################################

if [[ "${PIPELINE_MODE}" == "full" ||
      "${PIPELINE_MODE}" == "bioinformatics" ]]; then


    print_step_header "Phase I — Bioinformatics"


    ###########################################################################
    # Module 01 — Quality Control
    ###########################################################################

    if [[ "${RUN_QC}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 01 — Quality Control"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/01_quality_control.sh" \
                "${SAMPLE_ID}"

        done


        log_message \
            "INFO" \
            "Generating project-level Quality Control summary."


        bash \
            "${BASH_MODULE_DIR}/01_quality_control.sh" \
            summary

    else

        log_message \
            "INFO" \
            "Module 01 — Quality Control skipped."

    fi


    ###########################################################################
    # Module 02 — Host Depletion
    ###########################################################################

    if [[ "${RUN_HOST_DEPLETION}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 02 — Host Depletion"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/02_host_depletion.sh" \
                "${SAMPLE_ID}"

        done


        log_message \
            "INFO" \
            "Generating project-level Host Depletion summary."


        bash \
            "${BASH_MODULE_DIR}/02_host_depletion.sh" \
            summary

    else

        log_message \
            "INFO" \
            "Module 02 — Host Depletion skipped."

    fi


    ###########################################################################
    # Module 03 — Taxonomic Profiling
    ###########################################################################

    if [[ "${RUN_TAXONOMY}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 03 — Taxonomic Profiling"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/03_taxonomic_profiling.sh" \
                "${SAMPLE_ID}"

        done

    else

        log_message \
            "INFO" \
            "Module 03 — Taxonomic Profiling skipped."

    fi


    ###########################################################################
    # Module 04 — Read-based Resistome
    ###########################################################################

    if [[ "${RUN_READ_RESISTOME}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 04 — Read-based Resistome"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/04_read_based_resistome.sh" \
                "${SAMPLE_ID}"

        done

         log_message \
            "INFO" \
            "Generating sequencing-depth summary..."


        bash \
            "${BASH_MODULE_DIR}/04_read_based_resistome.sh" \
            summary

    else

        log_message \
            "INFO" \
            "Module 04 — Read-based Resistome skipped."

    fi


    ###########################################################################
    # Module 05 — Metagenome Assembly
    ###########################################################################

    if [[ "${RUN_ASSEMBLY}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 05 — Metagenome Assembly"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/05_metagenome_assembly.sh" \
                "${SAMPLE_ID}"

        done

    else

        log_message \
            "INFO" \
            "Module 05 — Metagenome Assembly skipped."

    fi


    ###########################################################################
    # Module 06 — Assembly Quality
    ###########################################################################

    if [[ "${RUN_ASSEMBLY_QC}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 06 — Assembly Quality"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/06_assembly_quality.sh" \
                "${SAMPLE_ID}"

        done

    else

        log_message \
            "INFO" \
            "Module 06 — Assembly Quality skipped."

    fi


    ###########################################################################
    # Module 07 — Contig Annotation
    ###########################################################################

    if [[ "${RUN_CONTIG_ANNOTATION}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 07 — Contig Annotation"


        for SAMPLE_ID in "${SAMPLE_IDS[@]}"
        do

            bash \
                "${BASH_MODULE_DIR}/07_contig_annotation.sh" \
                "${SAMPLE_ID}"

        done

    else

        log_message \
            "INFO" \
            "Module 07 — Contig Annotation skipped."

    fi


    print_step_footer "Phase I — Bioinformatics completed"

fi


###############################################################################
# 7. PHASE II — DOWNSTREAM ANALYSIS
###############################################################################
#
# Module 08 is deliberately treated as a single launcher.
#
# run_pipeline.sh does NOT directly call:
#
#   run_analysis.R
#   taxonomy_analysis.R
#   resistome_analysis.R
#   plasmid_analysis.R
#
# Those relationships are owned by:
#
#   08_launch_R_analysis.sh
#           |
#           +-- run_analysis.R
#                  |
#                  +-- taxonomy_analysis.R
#                  +-- resistome_analysis.R
#                  +-- plasmid_analysis.R
#
###############################################################################

if [[ "${PIPELINE_MODE}" == "full" ||
      "${PIPELINE_MODE}" == "analysis" ]]; then


    print_step_header "Phase II — Downstream Analysis"


    if [[ "${RUN_R_ANALYSIS}" == true ]]; then

        log_message \
            "INFO" \
            "Running Module 08 — Downstream Analysis"


        bash \
            "${BASH_MODULE_DIR}/08_launch_R_analysis.sh"


    else

        log_message \
            "INFO" \
            "Module 08 — Downstream Analysis skipped."

    fi


    print_step_footer "Phase II — Downstream Analysis completed"

fi


###############################################################################
# 8. PIPELINE SUMMARY
###############################################################################

print_step_header "Pipeline Summary"


log_message \
    "INFO" \
    "Execution mode : ${PIPELINE_MODE}"


if [[ "${PIPELINE_MODE}" != "analysis" ]]; then

    log_message \
        "INFO" \
        "Samples processed : ${#SAMPLE_IDS[@]}"

fi


echo


log_message "INFO" "Pipeline components executed:"


###############################################################################
# Phase I summary
###############################################################################

if [[ "${PIPELINE_MODE}" == "full" ||
      "${PIPELINE_MODE}" == "bioinformatics" ]]; then


    if [[ "${RUN_QC}" == true ]]; then
        log_message "INFO" "  ✓ Module 01 — Quality Control"
    fi


    if [[ "${RUN_HOST_DEPLETION}" == true ]]; then
        log_message "INFO" "  ✓ Module 02 — Host Depletion"
    fi


    if [[ "${RUN_TAXONOMY}" == true ]]; then
        log_message "INFO" "  ✓ Module 03 — Taxonomic Profiling"
    fi


    if [[ "${RUN_READ_RESISTOME}" == true ]]; then
        log_message "INFO" "  ✓ Module 04 — Read-based Resistome"
    fi


    if [[ "${RUN_ASSEMBLY}" == true ]]; then
        log_message "INFO" "  ✓ Module 05 — Metagenome Assembly"
    fi


    if [[ "${RUN_ASSEMBLY_QC}" == true ]]; then
        log_message "INFO" "  ✓ Module 06 — Assembly Quality"
    fi


    if [[ "${RUN_CONTIG_ANNOTATION}" == true ]]; then
        log_message "INFO" "  ✓ Module 07 — Contig Annotation"
    fi

fi


###############################################################################
# Phase II summary
###############################################################################

if [[ "${PIPELINE_MODE}" == "full" ||
      "${PIPELINE_MODE}" == "analysis" ]]; then


    if [[ "${RUN_R_ANALYSIS}" == true ]]; then

        log_message \
            "INFO" \
            "  ✓ Module 08 — Downstream Analysis"

    fi

fi


###############################################################################
# 9. COMPLETION
###############################################################################

TOTAL_RUNTIME="$(stop_timer "${PIPELINE_START_TIME}")"


echo


log_message \
    "INFO" \
    "Pipeline completed successfully."


log_message \
    "INFO" \
    "Total runtime : ${TOTAL_RUNTIME}"


log_message \
    "INFO" \
    "Results directory : ${RESULTS_DIR}"


print_step_footer "Pipeline Finished"


###############################################################################
# END OF run_pipeline.sh
###############################################################################
