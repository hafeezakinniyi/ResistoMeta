#!/usr/bin/env bash

###############################################################################
# RESISTOMETA
# Shared Bash Utility Functions
#
# Purpose:
#   Reusable helper functions supporting the Bash phase of the
#   Metagenomic AMR Pipeline.
#
# Responsibilities:
#   - logging
#   - directory management
#   - input/output validation
#   - software validation
#   - Conda environment handling
#   - command execution
#   - timing
#   - pipeline status reporting
###############################################################################

set -euo pipefail


###############################################################################
# 1. TIMESTAMP
###############################################################################

#------------------------------------------------------------------------------
# Generate formatted timestamp
#------------------------------------------------------------------------------

timestamp() {

    date "${DATE_FORMAT}"

}


###############################################################################
# 2. LOGGING
###############################################################################

#------------------------------------------------------------------------------
# Initialize pipeline logger
#------------------------------------------------------------------------------

initialize_logger() {

    create_directory "${LOG_DIR}"

    LOG_FILE="${LOG_DIR}/pipeline.log"

    exec > >(tee -a "${LOG_FILE}") 2>&1

    log_message "INFO" "Pipeline logging initialized."

}


#------------------------------------------------------------------------------
# Generic log message
#------------------------------------------------------------------------------

log_message() {

    local LEVEL="$1"
    shift

    echo "$(timestamp) [${LEVEL}] $*"

}


#------------------------------------------------------------------------------
# Warning message
#------------------------------------------------------------------------------

log_warning() {

    log_message "WARNING" "$@"

}


#------------------------------------------------------------------------------
# Error message
#------------------------------------------------------------------------------

log_error() {

    local MESSAGE="$*"

    echo "$(timestamp) [ERROR] ${MESSAGE}" >&2

}


###############################################################################
# 3. PIPELINE PRINTING
###############################################################################

#------------------------------------------------------------------------------
# Print pipeline banner
#------------------------------------------------------------------------------

pipeline_banner() {

    echo
    echo "=============================================================================="
    echo "${PIPELINE_NAME}"
    echo "=============================================================================="

    log_message "INFO" "Version            : ${PIPELINE_VERSION}"
    log_message "INFO" "Author             : ${PIPELINE_AUTHOR}"
    log_message "INFO" "Execution mode     : ${PIPELINE_MODE}"
    log_message "INFO" "Project root       : ${PROJECT_ROOT}"
    log_message "INFO" "Start time         : $(timestamp)"
    log_message "INFO" "CPU threads        : ${THREADS}"
    log_message "INFO" "R threads          : ${R_THREADS}"

    echo

}


#------------------------------------------------------------------------------
# Print section header
#------------------------------------------------------------------------------

print_step_header() {

    echo
    echo "=============================================================================="
    echo "$1"
    echo "=============================================================================="

}


#------------------------------------------------------------------------------
# Print section footer
#------------------------------------------------------------------------------

print_step_footer() {

    local STEP_NAME="$1"
    local STATUS="${2:-SUCCESS}"

    log_message "${STATUS}" "${STEP_NAME}."

}


###############################################################################
# 4. DIRECTORY UTILITIES
###############################################################################

#------------------------------------------------------------------------------
# Create directory if it does not exist
#------------------------------------------------------------------------------

create_directory() {

    local DIRECTORY="$1"

    [[ -d "${DIRECTORY}" ]] || mkdir -p "${DIRECTORY}"

}


#------------------------------------------------------------------------------
# Initialize pipeline directory structure
#------------------------------------------------------------------------------

initialize_directories() {

    local DIRECTORIES=(

        "${RESULTS_DIR}"

        # Bioinformatics outputs
        "${QC_DIR}"
        "${HOST_DEPLETION_DIR}"
        "${TAXONOMY_DIR}"
        "${READ_RESISTOME_DIR}"
        "${ASSEMBLY_DIR}"
        "${ASSEMBLY_QC_DIR}"
        "${CONTIG_ANNOTATION_DIR}"

        # Downstream R outputs
        "${ANALYSIS_DIR}"
        "${TAXONOMY_ANALYSIS_DIR}"
        "${RESISTOME_ANALYSIS_DIR}"
        "${PLASMID_ANALYSIS_DIR}"

        # Logs
        "${LOG_DIR}"

    )

    for DIR in "${DIRECTORIES[@]}"
    do
        create_directory "${DIR}"
    done

    log_message "INFO" "Pipeline directory structure initialized."

}


###############################################################################
# 5. FILE AND DIRECTORY VALIDATION
###############################################################################

#------------------------------------------------------------------------------
# Check directory exists
#------------------------------------------------------------------------------

check_directory_exists() {

    local DIRECTORY="$1"

    if [[ ! -d "${DIRECTORY}" ]]; then

        log_error "Directory not found: ${DIRECTORY}"

        exit 1

    fi

}


#------------------------------------------------------------------------------
# Check file exists
#------------------------------------------------------------------------------

check_file_exists() {

    local FILE="$1"

    if [[ ! -f "${FILE}" ]]; then

        log_error "File not found: ${FILE}"

        exit 1

    fi

}


#------------------------------------------------------------------------------
# Check multiple required inputs
#------------------------------------------------------------------------------

check_required_inputs() {

    for ITEM in "$@"
    do

        if [[ ! -e "${ITEM}" ]]; then

            log_error "Required input missing: ${ITEM}"

            exit 1

        fi

    done

}


#------------------------------------------------------------------------------
# Validate expected non-empty file output
#------------------------------------------------------------------------------

check_output() {

    local FILE="$1"

    if [[ ! -s "${FILE}" ]]; then

        log_error "Expected output not generated or empty: ${FILE}"

        exit 1

    fi

}


#------------------------------------------------------------------------------
# Validate paired-end reads
#------------------------------------------------------------------------------

check_paired_reads() {

    check_file_exists "$1"
    check_file_exists "$2"

}


###############################################################################
# 6. BOWTIE2 VALIDATION
###############################################################################

#------------------------------------------------------------------------------
# Check Bowtie2 index
#------------------------------------------------------------------------------

check_bowtie2_index() {

    local PREFIX="$1"

    if [[ -z "${PREFIX}" ]]; then

        log_error "Bowtie2 index prefix is empty."

        exit 1

    fi

    local REQUIRED=(
        "${PREFIX}.1.bt2"
        "${PREFIX}.2.bt2"
        "${PREFIX}.3.bt2"
        "${PREFIX}.4.bt2"
        "${PREFIX}.rev.1.bt2"
        "${PREFIX}.rev.2.bt2"
    )

    for FILE in "${REQUIRED[@]}"
    do

        if [[ ! -f "${FILE}" ]]; then

            log_error "Missing Bowtie2 index file: ${FILE}"

            exit 1

        fi

    done

    log_message "SUCCESS" "Bowtie2 index validated: ${PREFIX}"

}


###############################################################################
# 7. SOFTWARE VALIDATION
###############################################################################

#------------------------------------------------------------------------------
# Check software is installed
#------------------------------------------------------------------------------

check_tool() {

    local TOOL="$1"

    if ! command -v "${TOOL}" >/dev/null 2>&1; then

        log_error "Required software not found: ${TOOL}"

        exit 1

    fi

}


#------------------------------------------------------------------------------
# Check Conda is available
#------------------------------------------------------------------------------

check_conda() {

    if ! command -v conda >/dev/null 2>&1; then

        log_error "Conda executable not found."

        exit 1

    fi

}


#------------------------------------------------------------------------------
# Check Conda environment exists
#------------------------------------------------------------------------------

check_conda_environment() {

    local ENVIRONMENT="$1"

    check_conda

    if ! conda env list | awk '{print $1}' | grep -qx "${ENVIRONMENT}"; then

        log_error "Conda environment not found: ${ENVIRONMENT}"

        exit 1

    fi

}


###############################################################################
# 8. CONDA UTILITIES
###############################################################################

#------------------------------------------------------------------------------
# Initialize Conda shell functions
#------------------------------------------------------------------------------

initialize_conda() {

    local CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"

    if [[ ! -f "${CONDA_SH}" ]]; then

        log_error "Conda initialization script not found: ${CONDA_SH}"

        exit 1

    fi

    # shellcheck disable=SC1090
    source "${CONDA_SH}"

}


#------------------------------------------------------------------------------
# Activate Conda environment
#------------------------------------------------------------------------------

activate_env() {

    local ENVIRONMENT="$1"

    # Initialize Conda shell integration before using conda activate
    initialize_conda

    check_conda_environment "${ENVIRONMENT}"

    if [[ "${CONDA_DEFAULT_ENV:-}" != "${ENVIRONMENT}" ]]; then

        log_message "INFO" "Activating Conda environment: ${ENVIRONMENT}"

        conda activate "${ENVIRONMENT}"

    fi

}


#------------------------------------------------------------------------------
# Deactivate current Conda environment
#------------------------------------------------------------------------------

deactivate_env() {

    if [[ -n "${CONDA_DEFAULT_ENV:-}" ]] &&
       [[ "${CONDA_DEFAULT_ENV}" != "base" ]]; then

        conda deactivate

    fi

}


###############################################################################
# 9. COMMAND EXECUTION
###############################################################################

#------------------------------------------------------------------------------
# Execute command with logging and explicit error handling
#------------------------------------------------------------------------------

run_command() {

    log_message "INFO" "Running: $*"

    if "$@"; then

        return 0

    else

        local STATUS=$?

        log_error "Command failed (exit code: ${STATUS})"

        exit "${STATUS}"

    fi

}


#------------------------------------------------------------------------------
# Execute command and validate expected output
#------------------------------------------------------------------------------

run_command_checked() {

    local OUTPUT_FILE="$1"

    shift

    run_command "$@"

    check_output "${OUTPUT_FILE}"

}


###############################################################################
# 10. TIMING
###############################################################################

#------------------------------------------------------------------------------
# Start timer
#------------------------------------------------------------------------------

start_timer() {

    date +%s

}


#------------------------------------------------------------------------------
# Stop timer and return elapsed time
#------------------------------------------------------------------------------

stop_timer() {

    local START_TIME="$1"

    local END_TIME
    END_TIME=$(date +%s)

    local ELAPSED=$((END_TIME - START_TIME))

    printf "%02d:%02d:%02d" \
        $((ELAPSED / 3600)) \
        $(((ELAPSED % 3600) / 60)) \
        $((ELAPSED % 60))

}


###############################################################################
# 11. PIPELINE HELPERS
###############################################################################

#------------------------------------------------------------------------------
# Determine whether a pipeline step has already completed
#------------------------------------------------------------------------------

skip_step() {

    local STEP_NAME="$1"

    shift

    local COMPLETE=true

    for FILE in "$@"
    do

        if [[ ! -s "${FILE}" ]]; then

            COMPLETE=false

            break

        fi

    done

    if ${COMPLETE}; then

        log_message "INFO" \
            "Skipping ${STEP_NAME} (expected outputs already exist)."

        return 0

    fi

    return 1

}


#------------------------------------------------------------------------------
# Print pipeline summary
#------------------------------------------------------------------------------

pipeline_summary() {

    print_step_header "Pipeline Summary"

    log_message "INFO" "Pipeline           : ${PIPELINE_NAME}"
    log_message "INFO" "Version            : ${PIPELINE_VERSION}"
    log_message "INFO" "Execution mode     : ${PIPELINE_MODE}"
    log_message "INFO" "Finish time        : $(timestamp)"
    log_message "INFO" "Results directory  : ${RESULTS_DIR}"

}


#------------------------------------------------------------------------------
# Print pipeline completion message
#------------------------------------------------------------------------------

pipeline_completed() {

    local RUNTIME="$1"

    print_step_header "Pipeline Completed"

    log_message "SUCCESS" "Pipeline finished successfully."

    log_message "INFO" "Total runtime      : ${RUNTIME}"

}
