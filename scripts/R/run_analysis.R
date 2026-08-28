#!/usr/bin/env Rscript

###############################################################################
# RESISTOMETA
# Downstream R Analysis Orchestrator
#
# File:
#   run_analysis.R
#
# Purpose:
#   Coordinate downstream R analyses for:
#
#       1. Taxonomic / microbiome characterization
#       2. Read-based RGI-CARD resistome characterization
#       3. ABRicate-CARD / PlasmidFinder plasmid analysis
#
# Configuration architecture:
#
#       config/config.sh
#              |
#              v
#       08_launch_R_analysis.sh
#              |
#              v
#       run_analysis.R
#              |
#              +--> taxonomy_analysis.R
#              +--> resistome_analysis.R
#              +--> plasmid_analysis.R
#
# IMPORTANT:
#
#   config.sh remains the SINGLE SOURCE OF TRUTH.
#
#   Module 08 exports the configuration required by this workflow.
#
#   This orchestrator:
#
#       - validates exported configuration
#       - validates required upstream inputs
#       - creates analysis directories
#       - validates enabled analysis modules
#       - records the active configuration
#       - launches enabled analysis modules
#       - reports completion
#
#   This orchestrator does NOT:
#
#       - perform statistical analyses
#       - calculate abundance metrics
#       - perform normalization
#       - define analytical thresholds
#       - redefine pipeline switches
#
# Analytical logic belongs to the individual R modules.
#
###############################################################################


###############################################################################
# 1. STARTUP
###############################################################################

options(
    stringsAsFactors = FALSE,
    warn = 1
)


analysis_start_time <- Sys.time()


cat("\n")
cat("============================================================\n")
cat(" RESISTOMETA\n")
cat(" DOWNSTREAM R ANALYSIS\n")
cat("============================================================\n")
cat("\n")


cat(
    "Analysis started: ",
    format(analysis_start_time),
    "\n\n",
    sep = ""
)


###############################################################################
# 2. HELPER FUNCTIONS
###############################################################################

get_env <- function(
    name,
    required = TRUE,
    default = NULL
) {

    value <- Sys.getenv(
        name,
        unset = ""
    )

    value <- trimws(value)


    if (
        identical(value, "") &&
        !is.null(default)
    ) {

        return(default)

    }


    if (
        required &&
        identical(value, "")
    ) {

        stop(
            sprintf(
                "Required environment variable '%s' is not defined.",
                name
            )
        )

    }


    value

}


as_bool <- function(x, name = "Boolean configuration value") {

    if (
        length(x) != 1 ||
        is.na(x) ||
        identical(trimws(x), "")
    ) {

        stop(
            sprintf(
                "%s is missing.",
                name
            )
        )

    }


    value <- tolower(
        trimws(x)
    )


    valid_values <- c(
        "true",
        "false",
        "t",
        "f",
        "1",
        "0",
        "yes",
        "no",
        "y",
        "n"
    )


    if (!value %in% valid_values) {

        stop(
            sprintf(
                "Invalid Boolean configuration value for %s: '%s'.",
                name,
                x
            )
        )

    }


    value %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
    )

}


as_numeric <- function(x, name) {

    value <- suppressWarnings(
        as.numeric(x)
    )


    if (
        length(value) != 1 ||
        is.na(value) ||
        !is.finite(value)
    ) {

        stop(
            sprintf(
                "Invalid numeric configuration value for %s: '%s'.",
                name,
                x
            )
        )

    }


    value

}


as_integer <- function(x, name) {

    value <- suppressWarnings(
        as.numeric(x)
    )


    if (
        length(value) != 1 ||
        is.na(value) ||
        !is.finite(value) ||
        value != as.integer(value)
    ) {

        stop(
            sprintf(
                "Invalid integer configuration value for %s: '%s'.",
                name,
                x
            )
        )

    }


    as.integer(value)

}


log_info <- function(...) {

    cat(
        "[INFO] ",
        paste0(..., collapse = ""),
        "\n",
        sep = ""
    )

}


log_success <- function(...) {

    cat(
        "[SUCCESS] ",
        paste0(..., collapse = ""),
        "\n",
        sep = ""
    )

}


log_warning <- function(...) {

    cat(
        "[WARNING] ",
        paste0(..., collapse = ""),
        "\n",
        sep = ""
    )

}


###############################################################################
# 3. READ EXPORTED CONFIGURATION
###############################################################################

log_info(
    "Reading pipeline configuration exported by Module 08..."
)


###############################################################################
# Project paths
###############################################################################

PROJECT_ROOT <- get_env(
    "METAGENOMIC_AMR_PROJECT_ROOT"
)


ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_DIR"
)


###############################################################################
# Upstream result directories
###############################################################################

TAXONOMY_DIR <- get_env(
    "METAGENOMIC_AMR_TAXONOMY_DIR"
)


READ_RESISTOME_DIR <- get_env(
    "METAGENOMIC_AMR_READ_RESISTOME_DIR"
)


CONTIG_ANNOTATION_DIR <- get_env(
    "METAGENOMIC_AMR_CONTIG_ANNOTATION_DIR"
)


###############################################################################
# Downstream analysis directories
###############################################################################

TAXONOMY_ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_TAXONOMY_ANALYSIS_DIR"
)


RESISTOME_ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_RESISTOME_ANALYSIS_DIR"
)


PLASMID_ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_PLASMID_ANALYSIS_DIR"
)


FIGURE_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_FIGURE_DIR"
)


OBJECT_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_OBJECT_DIR"
)


###############################################################################
# R script/module directories
###############################################################################

R_MODULE_DIR <- get_env(
    "METAGENOMIC_AMR_R_MODULE_DIR"
)


###############################################################################
# Sequencing-depth input
#
# Required only when read-based resistome analysis is enabled.
###############################################################################

SEQUENCING_DEPTH_FILE <- get_env(
    "METAGENOMIC_AMR_SEQUENCING_DEPTH_FILE",
    required = FALSE,
    default = ""
)


###############################################################################
# 4. READ DOWNSTREAM R ANALYSIS SWITCHES
###############################################################################

RUN_R_ANALYSIS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_R_ANALYSIS"
    ),
    "METAGENOMIC_AMR_RUN_R_ANALYSIS"
)


RUN_TAXONOMY_ANALYSIS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_TAXONOMY_ANALYSIS"
    ),
    "METAGENOMIC_AMR_RUN_TAXONOMY_ANALYSIS"
)


RUN_RESISTOME_ANALYSIS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_RESISTOME_ANALYSIS"
    ),
    "METAGENOMIC_AMR_RUN_RESISTOME_ANALYSIS"
)


RUN_PLASMID_ANALYSIS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_PLASMID_ANALYSIS"
    ),
    "METAGENOMIC_AMR_RUN_PLASMID_ANALYSIS"
)


###############################################################################
# 5. READ UPSTREAM PRODUCTION SWITCHES
###############################################################################
#
# These switches describe whether the upstream data-producing modules were
# enabled.
#
# They are used for dependency validation and provenance.
#
###############################################################################

RUN_TAXONOMY <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_TAXONOMY"
    ),
    "METAGENOMIC_AMR_RUN_TAXONOMY"
)


RUN_READ_RESISTOME <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_READ_RESISTOME"
    ),
    "METAGENOMIC_AMR_RUN_READ_RESISTOME"
)


RUN_CONTIG_ANNOTATION <- as_bool(
    get_env(
        "METAGENOMIC_AMR_RUN_CONTIG_ANNOTATION"
    ),
    "METAGENOMIC_AMR_RUN_CONTIG_ANNOTATION"
)


###############################################################################
# 6. READ ANALYTICAL PARAMETERS
###############################################################################


###############################################################################
# Taxonomy parameters
###############################################################################

TOP_GENERA <- as_integer(
    get_env(
        "METAGENOMIC_AMR_TOP_GENERA"
    ),
    "METAGENOMIC_AMR_TOP_GENERA"
)


TOP_PHYLA <- as_integer(
    get_env(
        "METAGENOMIC_AMR_TOP_PHYLA"
    ),
    "METAGENOMIC_AMR_TOP_PHYLA"
)


TOP_SPECIES <- as_integer(
    get_env(
        "METAGENOMIC_AMR_TOP_SPECIES"
    ),
    "METAGENOMIC_AMR_TOP_SPECIES"
)


TAXONOMIC_LEVELS <- get_env(
    "METAGENOMIC_AMR_TAXONOMIC_LEVELS"
)


TAXONOMY_ABUNDANCE_SCALE <- get_env(
    "METAGENOMIC_AMR_TAXONOMY_ABUNDANCE_SCALE"
)


###############################################################################
# RGI resistome parameters
###############################################################################

RGI_MIN_MAPPED_READS <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_RGI_MIN_MAPPED_READS"
    ),
    "METAGENOMIC_AMR_RGI_MIN_MAPPED_READS"
)


RGI_MIN_COVERAGE <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_RGI_MIN_COVERAGE"
    ),
    "METAGENOMIC_AMR_RGI_MIN_COVERAGE"
)


RESISTOME_ABUNDANCE_METRIC <- get_env(
    "METAGENOMIC_AMR_RESISTOME_ABUNDANCE_METRIC"
)


###############################################################################
# ABRicate parameters
###############################################################################

ABRICATE_CARD_DB <- get_env(
    "METAGENOMIC_AMR_ABRICATE_CARD_DB"
)


ABRICATE_PLASMID_DB <- get_env(
    "METAGENOMIC_AMR_ABRICATE_PLASMID_DB"
)


ABRICATE_MIN_IDENTITY <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_ABRICATE_MIN_IDENTITY"
    ),
    "METAGENOMIC_AMR_ABRICATE_MIN_IDENTITY"
)


ABRICATE_MIN_COVERAGE <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_ABRICATE_MIN_COVERAGE"
    ),
    "METAGENOMIC_AMR_ABRICATE_MIN_COVERAGE"
)


###############################################################################
# Resistome calculation switches
###############################################################################

CALCULATE_ARG_RICHNESS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ARG_RICHNESS"
    ),
    "METAGENOMIC_AMR_CALCULATE_ARG_RICHNESS"
)


CALCULATE_ARG_PREVALENCE <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ARG_PREVALENCE"
    ),
    "METAGENOMIC_AMR_CALCULATE_ARG_PREVALENCE"
)


CALCULATE_ARG_CLASS_ABUNDANCE <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ARG_CLASS_ABUNDANCE"
    ),
    "METAGENOMIC_AMR_CALCULATE_ARG_CLASS_ABUNDANCE"
)


CALCULATE_ARG_CLASS_RICHNESS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ARG_CLASS_RICHNESS"
    ),
    "METAGENOMIC_AMR_CALCULATE_ARG_CLASS_RICHNESS"
)


CALCULATE_ABRICATE_RICHNESS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ABRICATE_RICHNESS"
    ),
    "METAGENOMIC_AMR_CALCULATE_ABRICATE_RICHNESS"
)


CALCULATE_ABRICATE_PREVALENCE <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_ABRICATE_PREVALENCE"
    ),
    "METAGENOMIC_AMR_CALCULATE_ABRICATE_PREVALENCE"
)


###############################################################################
# Plasmid calculation switches
###############################################################################

CALCULATE_PLASMID_RICHNESS <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_PLASMID_RICHNESS"
    ),
    "METAGENOMIC_AMR_CALCULATE_PLASMID_RICHNESS"
)


CALCULATE_PLASMID_PREVALENCE <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_PLASMID_PREVALENCE"
    ),
    "METAGENOMIC_AMR_CALCULATE_PLASMID_PREVALENCE"
)


CALCULATE_SAMPLE_CO_OCCURRENCE <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_SAMPLE_CO_OCCURRENCE"
    ),
    "METAGENOMIC_AMR_CALCULATE_SAMPLE_CO_OCCURRENCE"
)


CALCULATE_CONTIG_COLOCALIZATION <- as_bool(
    get_env(
        "METAGENOMIC_AMR_CALCULATE_CONTIG_COLOCALIZATION"
    ),
    "METAGENOMIC_AMR_CALCULATE_CONTIG_COLOCALIZATION"
)


###############################################################################
# 7. FIGURE PARAMETERS
###############################################################################
#
# These values are passed through from config.sh for downstream modules.
# run_analysis.R does not perform figure generation itself.
#
###############################################################################

GENERATE_FIGURES <- as_bool(
    get_env(
        "METAGENOMIC_AMR_GENERATE_FIGURES"
    ),
    "METAGENOMIC_AMR_GENERATE_FIGURES"
)


SAVE_PNG <- as_bool(
    get_env(
        "METAGENOMIC_AMR_SAVE_PNG"
    ),
    "METAGENOMIC_AMR_SAVE_PNG"
)


SAVE_PDF <- as_bool(
    get_env(
        "METAGENOMIC_AMR_SAVE_PDF"
    ),
    "METAGENOMIC_AMR_SAVE_PDF"
)


FIGURE_DPI <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_FIGURE_DPI"
    ),
    "METAGENOMIC_AMR_FIGURE_DPI"
)


FIGURE_WIDTH <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_FIGURE_WIDTH"
    ),
    "METAGENOMIC_AMR_FIGURE_WIDTH"
)


FIGURE_HEIGHT <- as_numeric(
    get_env(
        "METAGENOMIC_AMR_FIGURE_HEIGHT"
    ),
    "METAGENOMIC_AMR_FIGURE_HEIGHT"
)


###############################################################################
# 8. MASTER SWITCH
###############################################################################

if (!RUN_R_ANALYSIS) {

    log_info(
        "Downstream R analysis is disabled by configuration."
    )


    log_success(
        "Downstream R analysis skipped."
    )


    quit(
        save = "no",
        status = 0
    )

}


###############################################################################
# 9. DISPLAY ACTIVE CONFIGURATION
###############################################################################

log_info(
    "------------------------------------------------------------"
)

log_info(
    "Active downstream R analysis configuration"
)

log_info(
    "------------------------------------------------------------"
)


###############################################################################
# Downstream analysis switches
###############################################################################

log_info(
    "Taxonomy analysis: ",
    RUN_TAXONOMY_ANALYSIS
)


log_info(
    "Resistome analysis: ",
    RUN_RESISTOME_ANALYSIS
)


log_info(
    "Plasmid analysis: ",
    RUN_PLASMID_ANALYSIS
)


###############################################################################
# Upstream dependency switches
###############################################################################

log_info(
    "Upstream taxonomy available: ",
    RUN_TAXONOMY
)


log_info(
    "Upstream read resistome available: ",
    RUN_READ_RESISTOME
)


log_info(
    "Upstream contig annotation available: ",
    RUN_CONTIG_ANNOTATION
)


###############################################################################
# Analytical parameters
###############################################################################

log_info(
    "Top genera: ",
    TOP_GENERA
)


log_info(
    "Top phyla: ",
    TOP_PHYLA
)


log_info(
    "Top species: ",
    TOP_SPECIES
)


log_info(
    "Taxonomic levels: ",
    TAXONOMIC_LEVELS
)


log_info(
    "Taxonomy abundance scale: ",
    TAXONOMY_ABUNDANCE_SCALE
)


log_info(
    "RGI thresholds: ",
    RGI_MIN_MAPPED_READS,
    " mapped reads / ",
    RGI_MIN_COVERAGE,
    "% coverage"
)


log_info(
    "ABRicate thresholds: ",
    ABRICATE_MIN_IDENTITY,
    "% identity / ",
    ABRICATE_MIN_COVERAGE,
    "% coverage"
)


log_info(
    "Resistome abundance metric: ",
    RESISTOME_ABUNDANCE_METRIC
)


###############################################################################
# Figure parameters
###############################################################################

log_info(
    "Figure generation: ",
    GENERATE_FIGURES
)


log_info(
    "PNG output: ",
    SAVE_PNG
)


log_info(
    "PDF output: ",
    SAVE_PDF
)


log_info(
    "Figure dimensions: ",
    FIGURE_WIDTH,
    " x ",
    FIGURE_HEIGHT,
    " inches at ",
    FIGURE_DPI,
    " DPI"
)


###############################################################################
# 10. VALIDATE ANALYSIS DEPENDENCIES
###############################################################################

log_info(
    "Validating downstream analysis dependencies..."
)


###############################################################################
# Taxonomy
###############################################################################

if (RUN_TAXONOMY_ANALYSIS) {

    log_info(
        "Validating taxonomy analysis dependency..."
    )


    if (!RUN_TAXONOMY) {

        stop(
            paste0(
                "Taxonomy analysis is enabled, but upstream taxonomy ",
                "production (RUN_TAXONOMY) is disabled."
            )
        )

    }


    if (!dir.exists(TAXONOMY_DIR)) {

        stop(
            paste0(
                "Taxonomy results directory not found: ",
                TAXONOMY_DIR
            )
        )

    }

}


###############################################################################
# Read-based resistome
###############################################################################

if (RUN_RESISTOME_ANALYSIS) {

    log_info(
        "Validating read-based resistome analysis dependency..."
    )


    if (!RUN_READ_RESISTOME) {

        stop(
            paste0(
                "Resistome analysis is enabled, but upstream read-based ",
                "resistome production (RUN_READ_RESISTOME) is disabled."
            )
        )

    }


    if (!dir.exists(READ_RESISTOME_DIR)) {

        stop(
            paste0(
                "Read-based resistome directory not found: ",
                READ_RESISTOME_DIR
            )
        )

    }


    if (identical(SEQUENCING_DEPTH_FILE, "")) {

        stop(
            paste0(
                "Read-based resistome analysis requires a sequencing-depth ",
                "file, but METAGENOMIC_AMR_SEQUENCING_DEPTH_FILE is not defined."
            )
        )

    }


    if (!file.exists(SEQUENCING_DEPTH_FILE)) {

        stop(
            paste0(
                "Sequencing-depth file not found: ",
                SEQUENCING_DEPTH_FILE
            )
        )

    }

}


###############################################################################
# Plasmid analysis
###############################################################################

if (RUN_PLASMID_ANALYSIS) {

    log_info(
        "Validating plasmid analysis dependency..."
    )


    if (!RUN_CONTIG_ANNOTATION) {

        stop(
            paste0(
                "Plasmid analysis is enabled, but upstream contig annotation ",
                "production (RUN_CONTIG_ANNOTATION) is disabled."
            )
        )

    }


    if (!dir.exists(CONTIG_ANNOTATION_DIR)) {

        stop(
            paste0(
                "Contig annotation directory not found: ",
                CONTIG_ANNOTATION_DIR
            )
        )

    }

}


log_success(
    "All enabled downstream analysis dependencies validated."
)


###############################################################################
# 11. CREATE OUTPUT DIRECTORIES
###############################################################################

log_info(
    "Creating downstream analysis directories..."
)


analysis_directories <- unique(
    c(
        ANALYSIS_DIR,
        TAXONOMY_ANALYSIS_DIR,
        RESISTOME_ANALYSIS_DIR,
        PLASMID_ANALYSIS_DIR,
        FIGURE_DIR,
        OBJECT_DIR
    )
)


for (directory in analysis_directories) {

    if (!dir.exists(directory)) {

        dir.create(
            directory,
            recursive = TRUE,
            showWarnings = FALSE
        )

    }


    if (!dir.exists(directory)) {

        stop(
            paste0(
                "Failed to create analysis directory: ",
                directory
            )
        )

    }

}


###############################################################################
# 12. DEFINE R ANALYSIS MODULES
###############################################################################

TAXONOMY_MODULE <- file.path(
    R_MODULE_DIR,
    "taxonomy_analysis.R"
)


RESISTOME_MODULE <- file.path(
    R_MODULE_DIR,
    "resistome_analysis.R"
)


PLASMID_MODULE <- file.path(
    R_MODULE_DIR,
    "plasmid_analysis.R"
)


###############################################################################
# 13. VALIDATE ENABLED ANALYSIS MODULES
###############################################################################

log_info(
    "Validating enabled R analysis modules..."
)


enabled_modules <- character()


if (RUN_TAXONOMY_ANALYSIS) {

    enabled_modules <- c(
        enabled_modules,
        TAXONOMY_MODULE
    )

}


if (RUN_RESISTOME_ANALYSIS) {

    enabled_modules <- c(
        enabled_modules,
        RESISTOME_MODULE
    )

}


if (RUN_PLASMID_ANALYSIS) {

    enabled_modules <- c(
        enabled_modules,
        PLASMID_MODULE
    )

}


missing_modules <- enabled_modules[
    !file.exists(enabled_modules)
]


if (length(missing_modules) > 0) {

    stop(
        paste0(
            "The following enabled R analysis modules are missing:\n",
            paste(
                missing_modules,
                collapse = "\n"
            )
        )
    )

}


log_success(
    "All enabled R analysis modules are available."
)


###############################################################################
# 14. WRITE ACTIVE CONFIGURATION SUMMARY
###############################################################################

analysis_config_file <- file.path(
    ANALYSIS_DIR,
    "analysis_configuration.txt"
)


log_info(
    "Writing downstream analysis configuration..."
)


sink(
    analysis_config_file
)


cat(
    "RESISTOMETA\n"
)

cat(
    "Downstream R Analysis Configuration\n"
)

cat(
    "========================================\n\n"
)


cat(
    "Analysis started: ",
    format(analysis_start_time),
    "\n\n",
    sep = ""
)


###############################################################################
# Paths
###############################################################################

cat(
    "PATHS\n"
)

cat(
    "-----\n"
)


cat(
    "Project root: ",
    PROJECT_ROOT,
    "\n",
    sep = ""
)


cat(
    "Analysis directory: ",
    ANALYSIS_DIR,
    "\n",
    sep = ""
)


cat(
    "Taxonomy results: ",
    TAXONOMY_DIR,
    "\n",
    sep = ""
)


cat(
    "Read-based resistome: ",
    READ_RESISTOME_DIR,
    "\n",
    sep = ""
)


cat(
    "Contig annotation: ",
    CONTIG_ANNOTATION_DIR,
    "\n",
    sep = ""
)


cat(
    "Taxonomy analysis output: ",
    TAXONOMY_ANALYSIS_DIR,
    "\n",
    sep = ""
)


cat(
    "Resistome analysis output: ",
    RESISTOME_ANALYSIS_DIR,
    "\n",
    sep = ""
)


cat(
    "Plasmid analysis output: ",
    PLASMID_ANALYSIS_DIR,
    "\n",
    sep = ""
)


cat(
    "Figure output: ",
    FIGURE_DIR,
    "\n",
    sep = ""
)


cat(
    "R object output: ",
    OBJECT_DIR,
    "\n",
    sep = ""
)


cat(
    "Sequencing-depth file: ",
    ifelse(
        identical(SEQUENCING_DEPTH_FILE, ""),
        "not required / not configured",
        SEQUENCING_DEPTH_FILE
    ),
    "\n\n",
    sep = ""
)


###############################################################################
# Analysis switches
###############################################################################

cat(
    "DOWNSTREAM ANALYSIS SWITCHES\n"
)

cat(
    "----------------------------\n"
)


cat(
    "R analysis: ",
    RUN_R_ANALYSIS,
    "\n",
    sep = ""
)


cat(
    "Taxonomy analysis: ",
    RUN_TAXONOMY_ANALYSIS,
    "\n",
    sep = ""
)


cat(
    "Resistome analysis: ",
    RUN_RESISTOME_ANALYSIS,
    "\n",
    sep = ""
)


cat(
    "Plasmid analysis: ",
    RUN_PLASMID_ANALYSIS,
    "\n\n",
    sep = ""
)


###############################################################################
# Upstream dependencies
###############################################################################

cat(
    "UPSTREAM DEPENDENCIES\n"
)

cat(
    "---------------------\n"
)


cat(
    "Taxonomy production: ",
    RUN_TAXONOMY,
    "\n",
    sep = ""
)


cat(
    "Read-based resistome production: ",
    RUN_READ_RESISTOME,
    "\n",
    sep = ""
)


cat(
    "Contig annotation production: ",
    RUN_CONTIG_ANNOTATION,
    "\n\n",
    sep = ""
)


###############################################################################
# Analytical parameters
###############################################################################

cat(
    "ANALYTICAL PARAMETERS\n"
)

cat(
    "---------------------\n"
)


cat(
    "Top genera: ",
    TOP_GENERA,
    "\n",
    sep = ""
)


cat(
    "Top phyla: ",
    TOP_PHYLA,
    "\n",
    sep = ""
)


cat(
    "Top species: ",
    TOP_SPECIES,
    "\n",
    sep = ""
)


cat(
    "Taxonomic levels: ",
    TAXONOMIC_LEVELS,
    "\n",
    sep = ""
)


cat(
    "Taxonomy abundance scale: ",
    TAXONOMY_ABUNDANCE_SCALE,
    "\n",
    sep = ""
)


cat(
    "RGI minimum mapped reads: ",
    RGI_MIN_MAPPED_READS,
    "\n",
    sep = ""
)


cat(
    "RGI minimum coverage: ",
    RGI_MIN_COVERAGE,
    "%\n",
    sep = ""
)


cat(
    "Resistome abundance metric: ",
    RESISTOME_ABUNDANCE_METRIC,
    "\n",
    sep = ""
)


cat(
    "ABRicate CARD database: ",
    ABRICATE_CARD_DB,
    "\n",
    sep = ""
)


cat(
    "ABRicate PlasmidFinder database: ",
    ABRICATE_PLASMID_DB,
    "\n",
    sep = ""
)


cat(
    "ABRicate minimum identity: ",
    ABRICATE_MIN_IDENTITY,
    "%\n",
    sep = ""
)


cat(
    "ABRicate minimum coverage: ",
    ABRICATE_MIN_COVERAGE,
    "%\n\n",
    sep = ""
)


###############################################################################
# Figure parameters
###############################################################################

cat(
    "FIGURE PARAMETERS\n"
)

cat(
    "-----------------\n"
)


cat(
    "Generate figures: ",
    GENERATE_FIGURES,
    "\n",
    sep = ""
)


cat(
    "Save PNG: ",
    SAVE_PNG,
    "\n",
    sep = ""
)


cat(
    "Save PDF: ",
    SAVE_PDF,
    "\n",
    sep = ""
)


cat(
    "Figure DPI: ",
    FIGURE_DPI,
    "\n",
    sep = ""
)


cat(
    "Figure width: ",
    FIGURE_WIDTH,
    "\n",
    sep = ""
)


cat(
    "Figure height: ",
    FIGURE_HEIGHT,
    "\n\n",
    sep = ""
)


###############################################################################
# Calculation switches
###############################################################################

cat(
    "CALCULATION SWITCHES\n"
)

cat(
    "--------------------\n"
)


cat(
    "ARG richness: ",
    CALCULATE_ARG_RICHNESS,
    "\n",
    sep = ""
)


cat(
    "ARG prevalence: ",
    CALCULATE_ARG_PREVALENCE,
    "\n",
    sep = ""
)


cat(
    "ARG class abundance: ",
    CALCULATE_ARG_CLASS_ABUNDANCE,
    "\n",
    sep = ""
)


cat(
    "ARG class richness: ",
    CALCULATE_ARG_CLASS_RICHNESS,
    "\n",
    sep = ""
)


cat(
    "ABRicate richness: ",
    CALCULATE_ABRICATE_RICHNESS,
    "\n",
    sep = ""
)


cat(
    "ABRicate prevalence: ",
    CALCULATE_ABRICATE_PREVALENCE,
    "\n",
    sep = ""
)


cat(
    "Plasmid richness: ",
    CALCULATE_PLASMID_RICHNESS,
    "\n",
    sep = ""
)


cat(
    "Plasmid prevalence: ",
    CALCULATE_PLASMID_PREVALENCE,
    "\n",
    sep = ""
)


cat(
    "Sample-level ARG-plasmid co-occurrence: ",
    CALCULATE_SAMPLE_CO_OCCURRENCE,
    "\n",
    sep = ""
)


cat(
    "Contig-level ARG-plasmid co-localization: ",
    CALCULATE_CONTIG_COLOCALIZATION,
    "\n",
    sep = ""
)


sink()


log_success(
    "Analysis configuration summary written: ",
    analysis_config_file
)


###############################################################################
# 15. RUN TAXONOMIC / MICROBIOME ANALYSIS
###############################################################################

if (RUN_TAXONOMY_ANALYSIS) {

    log_info(
        "------------------------------------------------------------"
    )


    log_info(
        "Running taxonomic / microbiome analysis..."
    )


    log_info(
        "------------------------------------------------------------"
    )


    tryCatch(

        {

            source(
                TAXONOMY_MODULE,
                local = FALSE
            )

        },

        error = function(e) {

            stop(
                paste0(
                    "Taxonomic analysis failed: ",
                    conditionMessage(e)
                )
            )

        }

    )


    log_success(
        "Taxonomic analysis completed successfully."
    )

} else {

    log_info(
        "Taxonomic analysis disabled."
    )

}


###############################################################################
# 16. RUN READ-BASED RESISTOME ANALYSIS
###############################################################################

if (RUN_RESISTOME_ANALYSIS) {

    log_info(
        "------------------------------------------------------------"
    )


    log_info(
        "Running read-based RGI-CARD resistome analysis..."
    )


    log_info(
        "------------------------------------------------------------"
    )


    tryCatch(

        {

            source(
                RESISTOME_MODULE,
                local = FALSE
            )

        },

        error = function(e) {

            stop(
                paste0(
                    "Read-based resistome analysis failed: ",
                    conditionMessage(e)
                )
            )

        }

    )


    log_success(
        "Read-based resistome analysis completed successfully."
    )

} else {

    log_info(
        "Read-based resistome analysis disabled."
    )

}


###############################################################################
# 17. RUN PLASMID / ARG LOCALIZATION ANALYSIS
###############################################################################

if (RUN_PLASMID_ANALYSIS) {

    log_info(
        "------------------------------------------------------------"
    )


    log_info(
        "Running ABRicate-CARD / PlasmidFinder analysis..."
    )


    log_info(
        "------------------------------------------------------------"
    )


    tryCatch(

        {

            source(
                PLASMID_MODULE,
                local = FALSE
            )

        },

        error = function(e) {

            stop(
                paste0(
                    "Plasmid / ARG localization analysis failed: ",
                    conditionMessage(e)
                )
            )

        }

    )


    log_success(
        "Plasmid / ARG localization analysis completed successfully."
    )

} else {

    log_info(
        "Plasmid analysis disabled."
    )

}


###############################################################################
# 18. COMPLETION SUMMARY
###############################################################################

analysis_end_time <- Sys.time()


analysis_duration <- difftime(
    analysis_end_time,
    analysis_start_time,
    units = "mins"
)


cat("\n")
cat("============================================================\n")
cat(" DOWNSTREAM R ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("\n")


cat(
    "Started:  ",
    format(analysis_start_time),
    "\n",
    sep = ""
)


cat(
    "Finished: ",
    format(analysis_end_time),
    "\n",
    sep = ""
)


cat(
    "Duration: ",
    round(
        as.numeric(analysis_duration),
        2
    ),
    " minutes\n",
    sep = ""
)


cat(
    "\nAnalysis directory:\n",
    ANALYSIS_DIR,
    "\n\n",
    sep = ""
)


log_success(
    "All requested downstream analyses completed successfully."
)