#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Installation Validation
#
# Purpose:
#   Validates that the pipeline installation is complete before execution.
#
#   Checks:
#     - project architecture
#     - pipeline scripts
#     - Bash modules
#     - R workflow
#     - Conda environments
#     - required software
#     - databases
#     - sequencing input structure
#
# Usage:
#     bash validate_installation.sh
###############################################################################

set -euo pipefail


###############################################################################
# 1. DETERMINE PROJECT ROOT
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "${SCRIPT_DIR}")")"


###############################################################################
# 2. LOAD CONFIGURATION
###############################################################################

source "${PROJECT_ROOT}/config/config.sh"


###############################################################################
# 3. LOAD SHARED UTILITIES
###############################################################################

source "${SCRIPT_DIR_BASH}/utils.sh"


###############################################################################
# 4. INITIALIZE CONDA
###############################################################################

initialize_conda


###############################################################################
# 5. VALIDATE PROJECT ARCHITECTURE
###############################################################################

print_step_header "Installation Validation"

print_step_header "Validating Project Architecture"


###############################################################################
# Core project directories
###############################################################################

check_directory_exists "${PROJECT_ROOT}"

check_directory_exists "${CONFIG_DIR}"

check_directory_exists "${DATABASE_DIR}"

check_directory_exists "${DOCS_DIR}"

check_directory_exists "${LOG_DIR}"

check_directory_exists "${RAW_READS_DIR}"

check_directory_exists "${SCRIPT_DIR_BASH}"

check_directory_exists "${SCRIPT_DIR_R}"

check_directory_exists "${BASH_MODULE_DIR}"

check_directory_exists "${R_MODULE_DIR}"

###############################################################################
# Results directories
###############################################################################

check_directory_exists "${RESULTS_DIR}"

###############################################################################
# Pipeline entry scripts
###############################################################################

check_file_exists "${CONFIG_DIR}/config.sh"

check_file_exists "${RUN_PIPELINE_SCRIPT}"

check_file_exists "${VALIDATE_INSTALLATION_SCRIPT}"

check_file_exists "${RUN_ANALYSIS_SCRIPT}"


###############################################################################
# Bash modules
###############################################################################

BASH_MODULES=(

    "01_quality_control.sh"
    "02_host_depletion.sh"
    "03_taxonomic_profiling.sh"
    "04_read_based_resistome.sh"
    "05_metagenome_assembly.sh"
    "06_assembly_quality.sh"
    "07_contig_annotation.sh"
    "08_launch_R_analysis.sh"

)

for MODULE in "${BASH_MODULES[@]}"
do

    check_file_exists "${BASH_MODULE_DIR}/${MODULE}"

done


###############################################################################
# R modules
###############################################################################

R_MODULES=(

    "taxonomy_analysis.R"
    "resistome_analysis.R"
    "plasmid_analysis.R"

)

for MODULE in "${R_MODULES[@]}"
do

    check_file_exists "${R_MODULE_DIR}/${MODULE}"

done


print_step_footer "Project Architecture" "PASSED"


###############################################################################
# 6. VALIDATE CONDA ENVIRONMENTS
###############################################################################

print_step_header "Validating Conda Environments"


check_conda


check_conda_environment "${MAIN_ENV}"

check_conda_environment "${RGI_ENV}"

check_conda_environment "${R_ENV}"


print_step_footer "Conda Environments" "PASSED"


###############################################################################
# 7. VALIDATE BIOINFORMATICS SOFTWARE
###############################################################################

print_step_header "Validating Bioinformatics Software"


###############################################################################
# Main bioinformatics environment
###############################################################################

activate_env "${MAIN_ENV}"


MAIN_TOOLS=(

    "fastqc"
    "fastp"
    "multiqc"
    "bowtie2"
    "kraken2"
    "bracken"
    "megahit"
    "metaquast.py"
    "abricate"

)


for TOOL in "${MAIN_TOOLS[@]}"
do

    check_tool "${TOOL}"

done


###############################################################################
# RGI environment
###############################################################################

activate_env "${RGI_ENV}"

check_tool rgi


print_step_footer "Bioinformatics Software" "PASSED"


###############################################################################
# 8. VALIDATE DATABASES
###############################################################################

print_step_header "Validating Databases"


###############################################################################
# Host reference
###############################################################################

check_bowtie2_index "${HOST_REFERENCE_INDEX}"


###############################################################################
# Kraken2 database
###############################################################################

check_directory_exists "${KRAKEN_DATABASE}"


###############################################################################
# ABRicate databases
###############################################################################

activate_env "${MAIN_ENV}"

check_tool abricate


if ! abricate --list \
    | awk '{print $1}' \
    | grep -qx "${ABRICATE_CARD_DB}"; then

    log_error \
        "ABRicate database not available: ${ABRICATE_CARD_DB}"

    exit 1

fi


if ! abricate --list \
    | awk '{print $1}' \
    | grep -qx "${ABRICATE_PLASMID_DB}"; then

    log_error \
        "ABRicate database not available: ${ABRICATE_PLASMID_DB}"

    exit 1

fi


print_step_footer "Databases" "PASSED"


###############################################################################
# 9. VALIDATE R WORKFLOW
###############################################################################

print_step_header "Validating R Workflow"


activate_env "${R_ENV}"


###############################################################################
# Rscript
###############################################################################

check_tool Rscript


###############################################################################
# R entry point
###############################################################################

check_file_exists "${RUN_ANALYSIS_SCRIPT}"


###############################################################################
# R analysis modules
###############################################################################

for MODULE in "${R_MODULES[@]}"
do

    check_file_exists "${R_MODULE_DIR}/${MODULE}"

done


print_step_footer "R Workflow" "PASSED"


###############################################################################
# 10. VALIDATE SEQUENCING DATA
###############################################################################

print_step_header "Discovering Sequencing Samples"


mapfile -t SAMPLE_IDS < <(

    find "${RAW_READS_DIR}" \
        -type f \
        -name "*_R1_*.fastq.gz" \
        -exec basename {} \; |
    sed 's/_R1_.*//' |
    sort -u

)


if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then

    log_warning "No sequencing samples were detected."

else

    log_message \
        "INFO" \
        "Detected ${#SAMPLE_IDS[@]} sample(s)."

fi


print_step_footer \
    "Sequencing Sample Discovery" \
    "COMPLETED"


###############################################################################
# 11. INSTALLATION SUMMARY
###############################################################################

print_step_header "Installation Summary"


log_message "SUCCESS" "Project architecture validated."

log_message "SUCCESS" "Conda environments validated."

log_message "SUCCESS" "Bioinformatics software validated."

log_message "SUCCESS" "Databases validated."

log_message "SUCCESS" "R workflow validated."


if [[ ${#SAMPLE_IDS[@]} -gt 0 ]]; then

    log_message \
        "SUCCESS" \
        "Sequencing data detected (${#SAMPLE_IDS[@]} sample(s))."

else

    log_warning "No sequencing data detected."

fi


echo

log_message "INFO" "Pipeline Version        : ${PIPELINE_VERSION}"

log_message "INFO" "Pipeline Mode           : ${PIPELINE_MODE}"

log_message "INFO" "Bioinformatics Threads  : ${THREADS}"

log_message "INFO" "R Threads               : ${R_THREADS}"


###############################################################################
# 12. COMPLETION
###############################################################################

print_step_header "Installation Validation Completed"


log_message \
    "SUCCESS" \
    "Installation validation completed successfully."

echo

log_message \
    "INFO" \
    "Launch the complete pipeline using:"

echo

echo "    bash ${RUN_PIPELINE_SCRIPT}"

echo
