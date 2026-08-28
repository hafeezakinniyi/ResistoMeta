#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Module 08 - Launch R Analysis
#
# Purpose:
#   Transitions the pipeline from the Bash-based bioinformatics workflow to
#   the downstream R analytical workflow.
#
#   This module:
#
#       1. Validates the R analytical entry point.
#       2. Validates required upstream results for enabled R modules.
#       3. Validates sequencing-depth information when required.
#       4. Creates downstream analysis directories.
#       5. Activates the configured R environment.
#       6. Validates Rscript availability.
#       7. Exports configuration required by the R workflow.
#       8. Launches run_analysis.R.
#
#   Module 08 does NOT perform statistical analyses directly.
#
# Important:
#
#   config.sh remains the single source of truth for pipeline configuration.
#
#   Module 08 validates, exports, and passes configuration to the R workflow.
#   It does not redefine analytical parameters.
#
###############################################################################

set -euo pipefail


###############################################################################
# 1. DETERMINE PROJECT ROOT
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="$(dirname "$(dirname "$(dirname "${SCRIPT_DIR}")")")"


###############################################################################
# 2. LOAD CONFIGURATION
###############################################################################

source "${PROJECT_ROOT}/config/config.sh"

source "${SCRIPT_DIR_BASH}/utils.sh"


###############################################################################
# 3. MODULE HEADER
###############################################################################

print_step_header "Module 08 - Launch R Analysis"


###############################################################################
# 4. CHECK MASTER R ANALYSIS SWITCH
###############################################################################

if [[ "${RUN_R_ANALYSIS}" != true ]]; then

    log_message "INFO" \
        "Downstream R analysis is disabled in configuration."

    log_message "INFO" \
        "Module 08 skipped."

    exit 0

fi


###############################################################################
# 5. VALIDATE R ANALYTICAL WORKFLOW
###############################################################################

log_message "INFO" \
    "Validating R analytical workflow..."


###############################################################################
# Master R entry point
###############################################################################

check_file_exists \
    "${RUN_ANALYSIS_SCRIPT}"


###############################################################################
# Taxonomy module
###############################################################################

if [[ "${RUN_TAXONOMY_ANALYSIS}" == true ]]; then

    check_file_exists \
        "${TAXONOMY_ANALYSIS_SCRIPT}"

fi


###############################################################################
# Resistome module
###############################################################################

if [[ "${RUN_RESISTOME_ANALYSIS}" == true ]]; then

    check_file_exists \
        "${RESISTOME_ANALYSIS_SCRIPT}"

fi


###############################################################################
# Plasmid module
###############################################################################

if [[ "${RUN_PLASMID_ANALYSIS}" == true ]]; then

    check_file_exists \
        "${PLASMID_ANALYSIS_SCRIPT}"

fi


###############################################################################
# 6. VALIDATE REQUIRED UPSTREAM WORKFLOW DEPENDENCIES
###############################################################################
#
# Taxonomy:
#       Module 03
#
# Resistome:
#       Module 04 - RGI read-based resistome
#       Module 07 - ABRicate-CARD / contig annotation
#
# Plasmid:
#       Module 07 - ABRicate-CARD / PlasmidFinder
#
###############################################################################


###############################################################################
# 6A. TAXONOMY ANALYSIS
###############################################################################

if [[ "${RUN_TAXONOMY_ANALYSIS}" == true ]]; then

    if [[ "${RUN_TAXONOMY}" != true ]]; then

        log_error \
            "Taxonomy analysis is enabled but RUN_TAXONOMY is disabled."

        exit 1

    fi


    log_message "INFO" \
        "Checking taxonomic profiling results..."

    if [[ "${CHECK_INPUTS}" == true ]]; then

        check_directory_exists \
            "${TAXONOMY_DIR}"

    fi

fi


###############################################################################
# 6B. RESISTOME ANALYSIS
###############################################################################

if [[ "${RUN_RESISTOME_ANALYSIS}" == true ]]; then

    if [[ "${RUN_READ_RESISTOME}" != true ]]; then

        log_error \
            "Resistome analysis is enabled but RUN_READ_RESISTOME is disabled."

        exit 1

    fi


    log_message "INFO" \
        "Checking read-based resistome results..."

    if [[ "${CHECK_INPUTS}" == true ]]; then

        check_directory_exists \
            "${READ_RESISTOME_DIR}"


        log_message "INFO" \
            "Checking sequencing-depth information required for RPKM..."

        check_file_exists \
            "${SEQUENCING_DEPTH_FILE}"

    fi

fi


###############################################################################
# 6C. PLASMID ANALYSIS
###############################################################################

if [[ "${RUN_PLASMID_ANALYSIS}" == true ]]; then

    if [[ "${RUN_CONTIG_ANNOTATION}" != true ]]; then

        log_error \
            "Plasmid analysis is enabled but RUN_CONTIG_ANNOTATION is disabled."

        exit 1

    fi


    log_message "INFO" \
        "Checking contig annotation results..."

    if [[ "${CHECK_INPUTS}" == true ]]; then

        check_directory_exists \
            "${CONTIG_ANNOTATION_DIR}"


        #######################################################################
        # PlasmidFinder results
        #######################################################################

        PLASMIDFINDER_DIR="${CONTIG_ANNOTATION_DIR}/plasmidfinder"

        check_directory_exists \
            "${PLASMIDFINDER_DIR}"


        #######################################################################
        # ABRicate-CARD results
        #######################################################################

        ABRICATE_CARD_DIR="${CONTIG_ANNOTATION_DIR}/abricate_card"

        check_directory_exists \
            "${ABRICATE_CARD_DIR}"

    fi

fi


###############################################################################
# 7. CREATE DOWNSTREAM ANALYSIS DIRECTORIES
###############################################################################
#
# Directory definitions are provided by config.sh.
#
# Module 08 creates them but does not define their locations.
#
###############################################################################

log_message "INFO" \
    "Preparing downstream analysis directories..."


create_directory \
    "${ANALYSIS_DIR}"


create_directory \
    "${TAXONOMY_ANALYSIS_DIR}"


create_directory \
    "${RESISTOME_ANALYSIS_DIR}"


create_directory \
    "${PLASMID_ANALYSIS_DIR}"


create_directory \
    "${FIGURE_DIR}"


create_directory \
    "${OBJECT_DIR}"


###############################################################################
# 8. ACTIVATE R ENVIRONMENT
###############################################################################

log_message "INFO" \
    "Activating R environment: ${R_ENV}"


activate_env \
    "${R_ENV}"


###############################################################################
# 9. VALIDATE RSCRIPT
###############################################################################

RSCRIPT="$(command -v Rscript || true)"


if [[ -z "${RSCRIPT}" ]]; then

    log_error \
        "Rscript could not be found in the active R environment."

    exit 1

fi


###############################################################################
# 10. DISPLAY R VERSION
###############################################################################

log_message "INFO" \
    "Using Rscript: ${RSCRIPT}"


log_message \
    "INFO" \
    "R version: $(${RSCRIPT} --version 2>&1)"


###############################################################################
# 11. EXPORT R WORKFLOW CONFIGURATION
###############################################################################
#
# config.sh remains the single source of truth.
#
# Module 08 exports the configuration required by run_analysis.R and the
# downstream R modules.
#
###############################################################################


###############################################################################
# 11A. PROJECT PATHS
###############################################################################

export METAGENOMIC_AMR_PROJECT_ROOT="${PROJECT_ROOT}"

export METAGENOMIC_AMR_RESULTS_DIR="${RESULTS_DIR}"

export METAGENOMIC_AMR_ANALYSIS_DIR="${ANALYSIS_DIR}"


###############################################################################
# 11B. UPSTREAM RESULT DIRECTORIES
###############################################################################

export METAGENOMIC_AMR_TAXONOMY_DIR="${TAXONOMY_DIR}"

export METAGENOMIC_AMR_READ_RESISTOME_DIR="${READ_RESISTOME_DIR}"

export METAGENOMIC_AMR_CONTIG_ANNOTATION_DIR="${CONTIG_ANNOTATION_DIR}"


###############################################################################
# 11C. DOWNSTREAM ANALYSIS DIRECTORIES
###############################################################################

export METAGENOMIC_AMR_TAXONOMY_ANALYSIS_DIR="${TAXONOMY_ANALYSIS_DIR}"

export METAGENOMIC_AMR_RESISTOME_ANALYSIS_DIR="${RESISTOME_ANALYSIS_DIR}"

export METAGENOMIC_AMR_PLASMID_ANALYSIS_DIR="${PLASMID_ANALYSIS_DIR}"

export METAGENOMIC_AMR_ANALYSIS_FIGURE_DIR="${FIGURE_DIR}"

export METAGENOMIC_AMR_ANALYSIS_OBJECT_DIR="${OBJECT_DIR}"


###############################################################################
# 11D. R SCRIPT PATHS
###############################################################################

export METAGENOMIC_AMR_R_SCRIPT_DIR="${SCRIPT_DIR_R}"

export METAGENOMIC_AMR_R_MODULE_DIR="${R_MODULE_DIR}"

export METAGENOMIC_AMR_RUN_ANALYSIS_SCRIPT="${RUN_ANALYSIS_SCRIPT}"

export METAGENOMIC_AMR_TAXONOMY_ANALYSIS_SCRIPT="${TAXONOMY_ANALYSIS_SCRIPT}"

export METAGENOMIC_AMR_RESISTOME_ANALYSIS_SCRIPT="${RESISTOME_ANALYSIS_SCRIPT}"

export METAGENOMIC_AMR_PLASMID_ANALYSIS_SCRIPT="${PLASMID_ANALYSIS_SCRIPT}"


###############################################################################
# 11E. UPSTREAM MODULE STATUS
###############################################################################

export METAGENOMIC_AMR_RUN_TAXONOMY="${RUN_TAXONOMY}"

export METAGENOMIC_AMR_RUN_READ_RESISTOME="${RUN_READ_RESISTOME}"

export METAGENOMIC_AMR_RUN_CONTIG_ANNOTATION="${RUN_CONTIG_ANNOTATION}"


###############################################################################
# 11F. DOWNSTREAM R MODULE SWITCHES
###############################################################################

export METAGENOMIC_AMR_RUN_R_ANALYSIS="${RUN_R_ANALYSIS}"

export METAGENOMIC_AMR_RUN_TAXONOMY_ANALYSIS="${RUN_TAXONOMY_ANALYSIS}"

export METAGENOMIC_AMR_RUN_RESISTOME_ANALYSIS="${RUN_RESISTOME_ANALYSIS}"

export METAGENOMIC_AMR_RUN_PLASMID_ANALYSIS="${RUN_PLASMID_ANALYSIS}"

###############################################################################
# 11G. SEQUENCING-DEPTH INPUT
###############################################################################
#
# Used by resistome_analysis.R for:
#
#       - sequencing-depth normalization
#       - RPKM calculation
#       - resistome sample-universe definition
#
###############################################################################

export METAGENOMIC_AMR_SEQUENCING_DEPTH_FILE="${SEQUENCING_DEPTH_FILE}"


###############################################################################
# 11H. TAXONOMIC ANALYSIS PARAMETERS
###############################################################################

export METAGENOMIC_AMR_TOP_GENERA="${TOP_GENERA}"

export METAGENOMIC_AMR_TOP_PHYLA="${TOP_PHYLA}"

export METAGENOMIC_AMR_TOP_SPECIES="${TOP_SPECIES}"

export METAGENOMIC_AMR_TAXONOMIC_LEVELS="${TAXONOMIC_LEVELS}"

export METAGENOMIC_AMR_TAXONOMY_ABUNDANCE_SCALE="${TAXONOMY_ABUNDANCE_SCALE}"


###############################################################################
# 11I. RESISTOME ANALYSIS PARAMETERS
###############################################################################

export METAGENOMIC_AMR_RESISTOME_ABUNDANCE_METRIC="${RESISTOME_ABUNDANCE_METRIC}"

export METAGENOMIC_AMR_RGI_MIN_MAPPED_READS="${RGI_MIN_MAPPED_READS}"

export METAGENOMIC_AMR_RGI_MIN_COVERAGE="${RGI_MIN_COVERAGE}"

export METAGENOMIC_AMR_CALCULATE_ARG_RICHNESS="${CALCULATE_ARG_RICHNESS}"

export METAGENOMIC_AMR_CALCULATE_ARG_PREVALENCE="${CALCULATE_ARG_PREVALENCE}"

export METAGENOMIC_AMR_CALCULATE_ARG_CLASS_ABUNDANCE="${CALCULATE_ARG_CLASS_ABUNDANCE}"

export METAGENOMIC_AMR_CALCULATE_ARG_CLASS_RICHNESS="${CALCULATE_ARG_CLASS_RICHNESS}"


###############################################################################
# 11J. ABRICATE PARAMETERS
###############################################################################
#
# Module 07 is the authoritative stage for primary ABRicate hit filtering.
#
# These values are exported to R for configuration/provenance.
#
# Downstream R modules do NOT independently redefine or relax the upstream
# CARD/PlasmidFinder identity and coverage filtering.
#
###############################################################################

export METAGENOMIC_AMR_ABRICATE_CARD_DB="${ABRICATE_CARD_DB}"

export METAGENOMIC_AMR_ABRICATE_PLASMID_DB="${ABRICATE_PLASMID_DB}"

export METAGENOMIC_AMR_ABRICATE_MIN_IDENTITY="${ABRICATE_MIN_IDENTITY}"

export METAGENOMIC_AMR_ABRICATE_MIN_COVERAGE="${ABRICATE_MIN_COVERAGE}"

export METAGENOMIC_AMR_CALCULATE_ABRICATE_RICHNESS="${CALCULATE_ABRICATE_RICHNESS}"

export METAGENOMIC_AMR_CALCULATE_ABRICATE_PREVALENCE="${CALCULATE_ABRICATE_PREVALENCE}"


###############################################################################
# 11K. PLASMID ANALYSIS PARAMETERS
###############################################################################

export METAGENOMIC_AMR_CALCULATE_PLASMID_RICHNESS="${CALCULATE_PLASMID_RICHNESS}"

export METAGENOMIC_AMR_CALCULATE_PLASMID_PREVALENCE="${CALCULATE_PLASMID_PREVALENCE}"

export METAGENOMIC_AMR_CALCULATE_SAMPLE_CO_OCCURRENCE="${CALCULATE_SAMPLE_CO_OCCURRENCE}"

export METAGENOMIC_AMR_CALCULATE_CONTIG_COLOCALIZATION="${CALCULATE_CONTIG_COLOCALIZATION}"


###############################################################################
# 11L. FIGURE PARAMETERS
###############################################################################

export METAGENOMIC_AMR_GENERATE_FIGURES="${GENERATE_FIGURES}"

export METAGENOMIC_AMR_SAVE_PNG="${SAVE_PNG}"

export METAGENOMIC_AMR_SAVE_PDF="${SAVE_PDF}"

export METAGENOMIC_AMR_FIGURE_DPI="${FIGURE_DPI}"

export METAGENOMIC_AMR_FIGURE_WIDTH="${FIGURE_WIDTH}"

export METAGENOMIC_AMR_FIGURE_HEIGHT="${FIGURE_HEIGHT}"


###############################################################################
# 12. DISPLAY FINAL CONFIGURATION SUMMARY
###############################################################################

log_message "INFO" \
    "R analysis configuration exported successfully."


log_message "INFO" \
    "Project root       : ${PROJECT_ROOT}"


log_message "INFO" \
    "Taxonomy results   : ${TAXONOMY_DIR}"


log_message "INFO" \
    "Read resistome     : ${READ_RESISTOME_DIR}"


log_message "INFO" \
    "Contig annotation  : ${CONTIG_ANNOTATION_DIR}"


log_message "INFO" \
    "Analysis output    : ${ANALYSIS_DIR}"


log_message "INFO" \
    "Figure output      : ${FIGURE_DIR}"


log_message "INFO" \
    "R object output    : ${OBJECT_DIR}"


log_message "INFO" \
    "Sequencing depth   : ${SEQUENCING_DEPTH_FILE}"


log_message "INFO" \
    "ABRicate threshold : ${ABRICATE_MIN_IDENTITY}% identity / ${ABRICATE_MIN_COVERAGE}% coverage"


log_message "INFO" \
    "RGI threshold      : ${RGI_MIN_MAPPED_READS} mapped reads / ${RGI_MIN_COVERAGE}% coverage"


log_message "INFO" \
    "Resistome metric   : ${RESISTOME_ABUNDANCE_METRIC}"


log_message "INFO" \
    "R modules          : taxonomy=${RUN_TAXONOMY_ANALYSIS}, resistome=${RUN_RESISTOME_ANALYSIS}, plasmid=${RUN_PLASMID_ANALYSIS}"


###############################################################################
# 13. LAUNCH R ANALYSIS
###############################################################################

log_message "INFO" \
    "Launching downstream R analysis..."


run_command \
    "${RSCRIPT}" \
    "${RUN_ANALYSIS_SCRIPT}"


###############################################################################
# 14. COMPLETION
###############################################################################

print_step_footer \
    "Module 08 - Launch R Analysis"