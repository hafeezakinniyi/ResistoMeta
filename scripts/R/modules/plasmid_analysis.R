#!/usr/bin/env Rscript

###############################################################################
# RESISTOMETA
# Module: Plasmid Analysis
#
# File:
#   plasmid_analysis.R
#
# Purpose:
#   Perform generic downstream characterization of plasmid replicons detected
#   from assembled metagenomic contigs using ABRicate with the PlasmidFinder
#   database, together with ARG–plasmid association analysis using
#   ABRicate-CARD results.
#
# Scope:
#
#   1. Import ABRicate-PlasmidFinder results
#   2. Characterize plasmid replicon profiles
#   3. Deduplicate plasmid replicons within samples
#   4. Calculate plasmid richness
#   5. Calculate plasmid prevalence
#   6. Import ABRicate-CARD ARG results
#   7. Characterize ARG profiles
#   8. Calculate ARG–plasmid sample-level co-occurrence
#   9. Calculate ARG–plasmid contig-level co-localization
#  10. Generate integrated plasmid/ARG summaries
#  11. Save structured analysis objects
#
# This module intentionally does NOT perform:
#
#   - Metadata-dependent analysis
#   - Geographic comparisons
#   - Site comparisons
#   - Host-group comparisons
#   - Clinical/epidemiological comparisons
#   - Treatment-associated analysis
#   - Statistical association with external metadata
#   - Microbiome-resistome correlation
#   - Differential abundance analysis
#
# Sample identity is derived exclusively from upstream result filenames.
#
# IMPORTANT INTERPRETATION:
#
#   Sample-level ARG–plasmid co-occurrence:
#       ARG and plasmid replicon detected in the same biological sample.
#
#       This does NOT demonstrate physical linkage.
#
#   Contig-level ARG–plasmid co-localization:
#       ARG and plasmid replicon detected on the same assembled contig
#       within the same sample.
#
#       This provides evidence of sequence-level co-localization but does
#       not by itself establish complete plasmid structure, circularity,
#       mobility, or biological transfer.
#
# Upstream thresholds:
#
#   ABRicate identity and coverage thresholds are applied upstream.
#   This module analyzes the resulting ABRicate outputs directly.
#
###############################################################################


###############################################################################
# 1. LOAD REQUIRED PACKAGES
###############################################################################

required_packages <- c(
    "dplyr",
    "tidyr",
    "stringr",
    "readr",
    "purrr",
    "tibble"
)


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )
]


if (length(missing_packages) > 0) {

    stop(
        paste0(
            "The following R packages are required but not installed:\n",
            paste(
                missing_packages,
                collapse = ", "
            )
        )
    )

}


suppressPackageStartupMessages({

    library(dplyr)
    library(tidyr)
    library(stringr)
    library(readr)
    library(purrr)
    library(tibble)

})


###############################################################################
# 2. HELPER FUNCTIONS
###############################################################################

# IMPORTANT:
#
# This module is executed by run_analysis.R after Module 08 has exported the
# configuration from config.sh into the environment.
#
# Therefore:
#   - no analytical defaults are defined here;
#   - no configuration is sourced directly from config.sh;
#   - missing/invalid exported configuration causes an explicit failure.
#
# config.sh -> 08_launch_R_analysis.sh -> run_analysis.R -> this module
# remains the authoritative configuration flow.

get_env <- function(
    name,
    required = TRUE,
    default = NULL
) {

    value <- Sys.getenv(
        name,
        unset = ""
    )


    if (identical(value, "")) {

        if (required) {

            stop(
                paste0(
                    "Required exported environment variable is not set: ",
                    name
                )
            )

        }

        return(default)

    }


    value

}


safe_max <- function(
    x
) {

    x <- x[!is.na(x)]

    if (length(x) == 0) {
        return(NA_real_)
    }

    max(x)

}


as_bool <- function(
    x,
    name
) {

    if (
        length(x) == 0 ||
        is.na(x) ||
        identical(trimws(x), "")
    ) {

        stop(
            paste0(
                "Boolean configuration value is missing: ",
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
            paste0(
                "Invalid Boolean configuration value for ",
                name,
                ": ",
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


as_numeric <- function(
    x,
    name
) {

    value <- suppressWarnings(
        as.numeric(x)
    )


    if (
        length(value) != 1 ||
        is.na(value) ||
        !is.finite(value)
    ) {

        stop(
            paste0(
                "Invalid numeric configuration value for ",
                name,
                ": ",
                x
            )
        )

    }


    value

}


log_plasmid <- function(...) {

    message(
        paste0(
            "[PLASMID] ",
            paste0(..., collapse = "")
        )
    )

}


###############################################################################
# 3. READ EXPORTED R WORKFLOW CONFIGURATION
###############################################################################

message("============================================================")
message("PLASMID ANALYSIS")
message("============================================================")


###############################################################################
# Input paths exported by Module 08
###############################################################################

CONTIG_ANNOTATION_DIR <- get_env(
    "METAGENOMIC_AMR_CONTIG_ANNOTATION_DIR"
)


###############################################################################
# Output paths exported by Module 08
###############################################################################

PLASMID_ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_PLASMID_ANALYSIS_DIR"
)


OBJECT_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_OBJECT_DIR"
)


###############################################################################
# Plasmid analysis switches exported by Module 08
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
# ABRicate thresholds exported by Module 08
#
# These values are recorded for provenance only.
# Primary ABRicate filtering is performed upstream in Module 07.
###############################################################################

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
# 4. VALIDATE MODULE OUTPUT DIRECTORIES
###############################################################################
#
# Module 08 creates the downstream analysis directories before launching
# run_analysis.R. run_analysis.R also validates/creates its shared output
# directories. This module therefore validates the exported paths rather than
# redefining or creating them.
#
###############################################################################

if (!dir.exists(PLASMID_ANALYSIS_DIR)) {

    stop(
        paste0(
            "Plasmid analysis output directory not found: ",
            PLASMID_ANALYSIS_DIR
        )
    )

}


if (!dir.exists(OBJECT_DIR)) {

    stop(
        paste0(
            "R analysis object directory not found: ",
            OBJECT_DIR
        )
    )

}


###############################################################################
# 5. LOG ACTIVE MODULE CONFIGURATION
###############################################################################

log_plasmid(
    "Contig annotation directory: ",
    CONTIG_ANNOTATION_DIR
)


log_plasmid(
    "Output directory: ",
    PLASMID_ANALYSIS_DIR
)


log_plasmid(
    "R object directory: ",
    OBJECT_DIR
)


log_plasmid(
    "ABRicate upstream threshold: ",
    ABRICATE_MIN_IDENTITY,
    "% identity / ",
    ABRICATE_MIN_COVERAGE,
    "% coverage"
)


log_plasmid(
    "Plasmid richness: ",
    CALCULATE_PLASMID_RICHNESS
)


log_plasmid(
    "Plasmid prevalence: ",
    CALCULATE_PLASMID_PREVALENCE
)


log_plasmid(
    "Sample-level ARG-plasmid co-occurrence: ",
    CALCULATE_SAMPLE_CO_OCCURRENCE
)


log_plasmid(
    "Contig-level ARG-plasmid co-localization: ",
    CALCULATE_CONTIG_COLOCALIZATION
)


###############################################################################
# 6. VALIDATE INPUT DIRECTORY
###############################################################################

if (!dir.exists(CONTIG_ANNOTATION_DIR)) {

    stop(
        paste0(
            "Contig annotation directory not found:\n",
            CONTIG_ANNOTATION_DIR
        )
    )

}


###############################################################################
# 7. DEFINE UPSTREAM DIRECTORIES
###############################################################################

ABRICATE_CARD_DIR <- file.path(
    CONTIG_ANNOTATION_DIR,
    "abricate_card"
)


ABRICATE_PLASMIDFINDER_DIR <- file.path(
    CONTIG_ANNOTATION_DIR,
    "plasmidfinder"
)


if (!dir.exists(ABRICATE_PLASMIDFINDER_DIR)) {

    stop(
        paste0(
            "ABRicate-PlasmidFinder directory not found:\n",
            ABRICATE_PLASMIDFINDER_DIR
        )
    )

}


if (!dir.exists(ABRICATE_CARD_DIR)) {

    stop(
        paste0(
            "ABRicate-CARD directory not found:\n",
            ABRICATE_CARD_DIR
        )
    )

}


###############################################################################
# 8. LOCATE PLASMIDFINDER FILES
###############################################################################

message("")
log_plasmid(
    "Searching for ABRicate-PlasmidFinder results..."
)


plasmidfinder_files <- list.files(

    ABRICATE_PLASMIDFINDER_DIR,

    pattern = "_plasmidfinder\\.tsv$",

    full.names = TRUE,

    recursive = FALSE

)


if (length(plasmidfinder_files) == 0) {

    stop(
        paste0(
            "No ABRicate-PlasmidFinder files were found in:\n",
            ABRICATE_PLASMIDFINDER_DIR
        )
    )

}


log_plasmid(
    "PlasmidFinder files detected: ",
    length(plasmidfinder_files)
)


###############################################################################
# 9. DERIVE SAMPLE IDS FROM PLASMIDFINDER FILENAMES
###############################################################################

plasmid_sample_ids <- str_remove(

    basename(plasmidfinder_files),

    "_plasmidfinder\\.tsv$"

)


if (
    any(
        is.na(plasmid_sample_ids) |
        plasmid_sample_ids == ""
    )
) {

    stop(
        "One or more PlasmidFinder files produced an invalid SampleID."
    )

}


if (
    anyDuplicated(plasmid_sample_ids) > 0
) {

    duplicated_ids <- unique(
        plasmid_sample_ids[
            duplicated(plasmid_sample_ids)
        ]
    )


    stop(
        paste0(
            "Multiple PlasmidFinder files detected for the same SampleID:\n",
            paste(
                duplicated_ids,
                collapse = ", "
            )
        )
    )

}


all_plasmid_samples <- plasmid_sample_ids


n_plasmid_samples <- length(
    all_plasmid_samples
)


log_plasmid(
    "Samples represented by PlasmidFinder results: ",
    n_plasmid_samples
)


###############################################################################
# 10. READ ABRICATE-PLASMIDFINDER
###############################################################################

read_abricate_plasmidfinder <- function(file) {

    sample_id <- str_remove(

        basename(file),

        "_plasmidfinder\\.tsv$"

    )


    log_plasmid(
        "Reading PlasmidFinder: ",
        basename(file)
    )


    df <- read_tsv(

        file,

        show_col_types = FALSE,

        progress = FALSE,

        col_types = cols(

            START = col_double(),

            END = col_double(),

            `%COVERAGE` = col_double(),

            `%IDENTITY` = col_double(),

            .default = col_character()

        )

    )


    required_columns <- c(

        "#FILE",
        "SEQUENCE",
        "START",
        "END",
        "STRAND",
        "GENE",
        "%COVERAGE",
        "COVERAGE_MAP",
        "GAPS",
        "%IDENTITY",
        "DATABASE",
        "ACCESSION",
        "PRODUCT",
        "RESISTANCE"

    )


    missing_columns <- setdiff(

        required_columns,

        names(df)

    )


    if (length(missing_columns) > 0) {

        stop(
            paste0(
                "PlasmidFinder file ",
                basename(file),
                " is missing required columns:\n",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )

    }


    if (nrow(df) == 0) {

        return(

            tibble(

                SampleID = character(),

                Contig = character(),

                START = numeric(),

                END = numeric(),

                STRAND = character(),

                Plasmid = character(),

                Coverage = numeric(),

                Identity = numeric(),

                DATABASE = character(),

                ACCESSION = character(),

                PRODUCT = character()

            )

        )

    }


    df %>%

        transmute(

            SampleID =
                sample_id,

            Contig =
                str_trim(
                    as.character(SEQUENCE)
                ),

            START =
                as.numeric(START),

            END =
                as.numeric(END),

            STRAND =
                as.character(STRAND),

            Plasmid =
                str_trim(
                    as.character(GENE)
                ),

            Coverage =
                as.numeric(`%COVERAGE`),

            Identity =
                as.numeric(`%IDENTITY`),

            DATABASE =
                as.character(DATABASE),

            ACCESSION =
                as.character(ACCESSION),

            PRODUCT =
                as.character(PRODUCT)

        ) %>%

        filter(

            !is.na(Contig),
            Contig != "",

            !is.na(Plasmid),
            Plasmid != ""

        )

}


###############################################################################
# 11. IMPORT PLASMIDFINDER RESULTS
###############################################################################

plasmidfinder_all <- map_dfr(

    plasmidfinder_files,

    read_abricate_plasmidfinder

)


log_plasmid(
    "PlasmidFinder hits detected: ",
    nrow(plasmidfinder_all)
)


###############################################################################
# 12. SAVE PARSED PLASMIDFINDER RESULTS
###############################################################################

write_csv(

    plasmidfinder_all,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "plasmidfinder_parsed_hits.csv"

    )

)


###############################################################################
# 13. WITHIN-SAMPLE PLASMID DEDUPLICATION
###############################################################################
#
# A plasmid replicon may occur on multiple contigs.
#
# For sample-level plasmid characterization, each replicon is counted once
# per sample.
#
###############################################################################

plasmid_sample <- plasmidfinder_all %>%

    group_by(

        SampleID,

        Plasmid

    ) %>%

    summarise(

        n_contigs =
            n_distinct(
                Contig
            ),

        max_identity =
            safe_max(
                Identity
            ),

        max_coverage =
            safe_max(
                Coverage
            ),

        .groups = "drop"

    )


###############################################################################
# 14. SAVE SAMPLE × PLASMID TABLE
###############################################################################

write_csv(

    plasmid_sample,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "plasmid_sample_table.csv"

    )

)


###############################################################################
# 15. PLASMID RICHNESS
###############################################################################
#
# Every sample represented by a PlasmidFinder file is retained.
#
# Samples without plasmid hits receive richness = 0.
#
###############################################################################

if (CALCULATE_PLASMID_RICHNESS) {

    message("")
    log_plasmid(
        "Calculating plasmid richness..."
    )


    plasmid_richness <- tibble(

        SampleID =
            all_plasmid_samples

    ) %>%

        left_join(

            plasmid_sample %>%

                group_by(
                    SampleID
                ) %>%

                summarise(

                    Plasmid_Richness =
                        n_distinct(
                            Plasmid
                        ),

                    .groups = "drop"

                ),

            by = "SampleID"

        ) %>%

        mutate(

            Plasmid_Richness =
                replace_na(
                    Plasmid_Richness,
                    0
                )

        )


    write_csv(

        plasmid_richness,

        file.path(

            PLASMID_ANALYSIS_DIR,

            "plasmid_richness.csv"

        )

    )

} else {

    plasmid_richness <- tibble()

}


###############################################################################
# 16. PLASMID PREVALENCE
###############################################################################
#
# Denominator:
#
#   Number of samples represented by PlasmidFinder result files.
#
###############################################################################

if (CALCULATE_PLASMID_PREVALENCE) {

    message("")
    log_plasmid(
        "Calculating plasmid prevalence..."
    )


    plasmid_prevalence <- plasmid_sample %>%

        group_by(
            Plasmid
        ) %>%

        summarise(

            n_samples =
                n_distinct(
                    SampleID
                ),

            .groups = "drop"

        ) %>%

        mutate(

            prevalence =
                n_samples /
                n_plasmid_samples,

            prevalence_percent =
                prevalence *
                100

        ) %>%

        arrange(

            desc(prevalence),

            desc(n_samples),

            Plasmid

        )


    write_csv(

        plasmid_prevalence,

        file.path(

            PLASMID_ANALYSIS_DIR,

            "plasmid_prevalence.csv"

        )

    )

} else {

    plasmid_prevalence <- tibble()

}


###############################################################################
# 17. PLASMID PROFILE SUMMARY
###############################################################################
#
# Provides an overall characterization of detected plasmid replicons.
#
###############################################################################

plasmid_profile <- plasmid_sample %>%

    group_by(

        Plasmid

    ) %>%

    summarise(

        n_samples =
            n_distinct(
                SampleID
            ),

        n_contigs =
            sum(
                n_contigs,
                na.rm = TRUE
            ),

        max_identity =
            safe_max(
                max_identity
            ),

        max_coverage =
            safe_max(
                max_coverage
            ),

        .groups = "drop"

    ) %>%

    mutate(

        prevalence =
            n_samples /
            n_plasmid_samples,

        prevalence_percent =
            prevalence *
            100

    ) %>%

    arrange(

        desc(n_samples),

        desc(n_contigs),

        Plasmid

    )


write_csv(

    plasmid_profile,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "plasmid_profile.csv"

    )

)


###############################################################################
# 18. LOCATE ABRICATE-CARD FILES
###############################################################################

message("")
log_plasmid(
    "Searching for ABRicate-CARD results..."
)


abricate_card_files <- list.files(

    ABRICATE_CARD_DIR,

    pattern = "_card\\.tsv$",

    full.names = TRUE,

    recursive = FALSE

)


if (length(abricate_card_files) == 0) {

    stop(
        paste0(
            "No ABRicate-CARD files were found in:\n",
            ABRICATE_CARD_DIR
        )
    )

}


log_plasmid(
    "ABRicate-CARD files detected: ",
    length(abricate_card_files)
)


###############################################################################
# 19. READ ABRICATE-CARD
###############################################################################

read_abricate_card <- function(file) {

    sample_id <- str_remove(

        basename(file),

        "_card\\.tsv$"

    )


    log_plasmid(
        "Reading ABRicate-CARD: ",
        basename(file)
    )


    df <- read_tsv(

        file,

        show_col_types = FALSE,

        progress = FALSE,

        col_types = cols(

            START = col_double(),

            END = col_double(),

            `%COVERAGE` = col_double(),

            `%IDENTITY` = col_double(),

            .default = col_character()

        )

    )


    required_columns <- c(

        "#FILE",
        "SEQUENCE",
        "START",
        "END",
        "STRAND",
        "GENE",
        "%COVERAGE",
        "COVERAGE_MAP",
        "GAPS",
        "%IDENTITY",
        "DATABASE",
        "ACCESSION",
        "PRODUCT",
        "RESISTANCE"

    )


    missing_columns <- setdiff(

        required_columns,

        names(df)

    )


    if (length(missing_columns) > 0) {

        stop(
            paste0(
                "ABRicate-CARD file ",
                basename(file),
                " is missing required columns:\n",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )

    }


    if (nrow(df) == 0) {

        return(

            tibble(

                SampleID = character(),

                Contig = character(),

                START = numeric(),

                END = numeric(),

                STRAND = character(),

                ARG = character(),

                Coverage = numeric(),

                Identity = numeric(),

                DATABASE = character(),

                ACCESSION = character(),

                PRODUCT = character(),

                RESISTANCE = character()

            )

        )

    }


    df %>%

        transmute(

            SampleID =
                sample_id,

            Contig =
                str_trim(
                    as.character(SEQUENCE)
                ),

            START =
                as.numeric(START),

            END =
                as.numeric(END),

            STRAND =
                as.character(STRAND),

            ARG =
                str_trim(
                    as.character(GENE)
                ),

            Coverage =
                as.numeric(`%COVERAGE`),

            Identity =
                as.numeric(`%IDENTITY`),

            DATABASE =
                as.character(DATABASE),

            ACCESSION =
                as.character(ACCESSION),

            PRODUCT =
                as.character(PRODUCT),

            RESISTANCE =
                as.character(RESISTANCE)

        ) %>%

        filter(

            !is.na(Contig),
            Contig != "",

            !is.na(ARG),
            ARG != ""

        )

}


###############################################################################
# 20. IMPORT ABRICATE-CARD RESULTS
###############################################################################

abricate_card_all <- map_dfr(

    abricate_card_files,

    read_abricate_card

)


log_plasmid(
    "ABRicate-CARD ARG hits detected: ",
    nrow(abricate_card_all)
)


###############################################################################
# 21. DEFINE CARD SAMPLE COHORT
###############################################################################

card_sample_ids <- str_remove(

    basename(abricate_card_files),

    "_card\\.tsv$"

)


n_card_samples <- length(
    card_sample_ids
)


log_plasmid(
    "Samples represented by CARD results: ",
    n_card_samples
)


###############################################################################
# 22. SAVE PARSED CARD RESULTS
###############################################################################

write_csv(

    abricate_card_all,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "abricate_card_parsed_hits.csv"

    )

)


###############################################################################
# 23. DEDUPLICATE ARG–CONTIG OBSERVATIONS
###############################################################################
#
# Multiple overlapping ABRicate hits for the same ARG on the same contig
# are collapsed to one observation.
#
###############################################################################

arg_contig <- abricate_card_all %>%

    group_by(

        SampleID,

        Contig,

        ARG

    ) %>%

    summarise(

        max_identity =
            safe_max(
                Identity
            ),

        max_coverage =
            safe_max(
                Coverage
            ),

        .groups = "drop"

    )


write_csv(

    arg_contig,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "abricate_arg_contig_table.csv"

    )

)


###############################################################################
# 24. ARG PROFILE SUMMARY
###############################################################################

arg_profile <- abricate_card_all %>%

    distinct(

        SampleID,

        ARG

    ) %>%

    group_by(

        ARG

    ) %>%

    summarise(

        n_samples =
            n_distinct(
                SampleID
            ),

        .groups = "drop"

    ) %>%

    mutate(

        prevalence =
            n_samples /
            n_card_samples,

        prevalence_percent =
            prevalence *
            100

    ) %>%

    arrange(

        desc(n_samples),

        ARG

    )


write_csv(

    arg_profile,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "arg_profile.csv"

    )

)


###############################################################################
# 25. ARG–PLASMID SAMPLE-LEVEL CO-OCCURRENCE
###############################################################################
#
# Definition:
#
#   ARG and plasmid replicon are detected in the same biological sample.
#
# This does NOT demonstrate physical linkage.
#
# Denominator:
#
#   Number of samples represented by the PlasmidFinder result files.
#
###############################################################################

if (CALCULATE_SAMPLE_CO_OCCURRENCE) {

    message("")
    log_plasmid(
        "Calculating ARG–plasmid sample-level co-occurrence..."
    )


    arg_sample <- arg_contig %>%

        distinct(

            SampleID,

            ARG

        )


    plasmid_sample_pairs <- plasmid_sample %>%

        select(

            SampleID,

            Plasmid

        ) %>%

        distinct()


    arg_plasmid_cooccurrence <- arg_sample %>%

        inner_join(

            plasmid_sample_pairs,

            by = "SampleID"

        ) %>%

        group_by(

            ARG,

            Plasmid

        ) %>%

        summarise(

            n_samples =
                n_distinct(
                    SampleID
                ),

            prevalence =
                n_samples /
                n_plasmid_samples,

            prevalence_percent =
                prevalence *
                100,

            .groups = "drop"

        ) %>%

        arrange(

            desc(n_samples),

            ARG,

            Plasmid

        )


    write_csv(

        arg_plasmid_cooccurrence,

        file.path(

            PLASMID_ANALYSIS_DIR,

            "arg_plasmid_sample_cooccurrence.csv"

        )

    )

} else {

    arg_plasmid_cooccurrence <- tibble()

}


###############################################################################
# 26. ARG–PLASMID CONTIG-LEVEL CO-LOCALIZATION
###############################################################################
#
# Definition:
#
#   ARG and plasmid replicon occur on the same assembled contig
#   within the same sample.
#
# This is sequence-level co-localization.
#
###############################################################################

if (CALCULATE_CONTIG_COLOCALIZATION) {

    message("")
    log_plasmid(
        "Calculating ARG–plasmid contig-level co-localization..."
    )


    plasmid_contig <- plasmidfinder_all %>%

        select(

            SampleID,

            Contig,

            Plasmid

        ) %>%

        distinct()


    arg_plasmid_contig_pairs <- arg_contig %>%

        inner_join(

            plasmid_contig,

            by = c(
                "SampleID",
                "Contig"
            )

        ) %>%

        distinct(

            SampleID,

            Contig,

            ARG,

            Plasmid

        )


    arg_plasmid_colocalization <- arg_plasmid_contig_pairs %>%

        group_by(

            ARG,

            Plasmid

        ) %>%

        summarise(

            n_samples =
                n_distinct(
                    SampleID
                ),

            n_contigs =
                n_distinct(
                    paste(
                        SampleID,
                        Contig,
                        sep = "::"
                    )
                ),

            .groups = "drop"

        ) %>%

        arrange(

            desc(n_samples),

            desc(n_contigs),

            ARG,

            Plasmid

        )


    write_csv(

        arg_plasmid_colocalization,

        file.path(

            PLASMID_ANALYSIS_DIR,

            "arg_plasmid_contig_colocalization.csv"

        )

    )


    ###########################################################################
    # Detailed same-contig observations
    ###########################################################################

    write_csv(

        arg_plasmid_contig_pairs,

        file.path(

            PLASMID_ANALYSIS_DIR,

            "arg_plasmid_contig_colocalization_details.csv"

        )

    )

} else {

    arg_plasmid_colocalization <- tibble()

    arg_plasmid_contig_pairs <- tibble()

}


###############################################################################
# 27. SAMPLE-LEVEL PLASMID SUMMARY
###############################################################################
#
# This summary is based exclusively on PlasmidFinder result files.
#
###############################################################################

plasmid_summary <- plasmid_richness


if (nrow(plasmid_summary) > 0) {

    plasmid_summary <- plasmid_summary %>%

        arrange(
            SampleID
        )

}


write_csv(

    plasmid_summary,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "plasmid_sample_summary.csv"

    )

)


###############################################################################
# 28. OVERALL PLASMID SUMMARY
###############################################################################

overall_plasmid_summary <- tibble(

    n_samples =
        n_plasmid_samples,

    n_samples_with_plasmid_hits =
        n_distinct(
            plasmid_sample$SampleID
        ),

    n_samples_without_plasmid_hits =
        n_plasmid_samples -
        n_distinct(
            plasmid_sample$SampleID
        ),

    n_unique_plasmid_replicons =
        n_distinct(
            plasmidfinder_all$Plasmid
        ),

    n_plasmid_hit_records =
        nrow(
            plasmidfinder_all
        ),

    n_sample_plasmid_observations =
        nrow(
            plasmid_sample
        )

)


write_csv(

    overall_plasmid_summary,

    file.path(

        PLASMID_ANALYSIS_DIR,

        "plasmid_overall_summary.csv"

    )

)


###############################################################################
# 29. SAVE STRUCTURED R OBJECT
###############################################################################

plasmid_analysis <- list(

    configuration = list(

        plasmidfinder_directory =
            ABRICATE_PLASMIDFINDER_DIR,

        card_directory =
            ABRICATE_CARD_DIR,

        output_directory =
            PLASMID_ANALYSIS_DIR,

        object_directory =
            OBJECT_DIR,

        abricate_min_identity =
            ABRICATE_MIN_IDENTITY,

        abricate_min_coverage =
            ABRICATE_MIN_COVERAGE,

        calculate_plasmid_richness =
            CALCULATE_PLASMID_RICHNESS,

        calculate_plasmid_prevalence =
            CALCULATE_PLASMID_PREVALENCE,

        calculate_sample_cooccurrence =
            CALCULATE_SAMPLE_CO_OCCURRENCE,

        calculate_contig_colocalization =
            CALCULATE_CONTIG_COLOCALIZATION

    ),

    sample_information = list(

        plasmidfinder_samples =
            all_plasmid_samples,

        n_plasmidfinder_samples =
            n_plasmid_samples,

        card_samples =
            card_sample_ids,

        n_card_samples =
            n_card_samples

    ),

    plasmid = list(

        raw_hits =
            plasmidfinder_all,

        sample_table =
            plasmid_sample,

        richness =
            plasmid_richness,

        prevalence =
            plasmid_prevalence,

        profile =
            plasmid_profile,

        sample_summary =
            plasmid_summary,

        overall_summary =
            overall_plasmid_summary

    ),

    resistance = list(

        raw_hits =
            abricate_card_all,

        arg_contig =
            arg_contig,

        profile =
            arg_profile

    ),

    associations = list(

        sample_cooccurrence =
            arg_plasmid_cooccurrence,

        contig_colocalization =
            arg_plasmid_colocalization,

        contig_details =
            arg_plasmid_contig_pairs

    ),

    analysis = list(

        generated =
            Sys.time(),

        description =
            paste(
                "Generic plasmid replicon characterization and",
                "ARG-plasmid association analysis."
            )

    )

)


saveRDS(

    plasmid_analysis,

    file.path(

        OBJECT_DIR,

        "plasmid_analysis.rds"

    )

)


###############################################################################
# 30. WRITE ANALYSIS SUMMARY
###############################################################################

summary_file <- file.path(

    PLASMID_ANALYSIS_DIR,

    "plasmid_analysis_summary.txt"

)


sink(summary_file)


cat(
    "METAGENOMIC AMR PIPELINE\n"
)

cat(
    "Plasmid Analysis Summary\n"
)

cat(
    "========================================\n\n"
)


cat(
    "Analysis date: ",
    format(Sys.time()),
    "\n\n",
    sep = ""
)


cat(
    "PlasmidFinder samples: ",
    n_plasmid_samples,
    "\n",
    sep = ""
)


cat(
    "Samples with plasmid hits: ",
    n_distinct(
        plasmid_sample$SampleID
    ),
    "\n",
    sep = ""
)


cat(
    "Samples without plasmid hits: ",
    n_plasmid_samples -
        n_distinct(
            plasmid_sample$SampleID
        ),
    "\n",
    sep = ""
)


cat(
    "ABRicate identity threshold: ",
    ABRICATE_MIN_IDENTITY,
    "%\n",
    sep = ""
)


cat(
    "ABRicate coverage threshold: ",
    ABRICATE_MIN_COVERAGE,
    "%\n\n",
    sep = ""
)


cat(
    "PlasmidFinder hit records: ",
    nrow(plasmidfinder_all),
    "\n",
    sep = ""
)


cat(
    "Unique plasmid replicons: ",
    n_distinct(
        plasmidfinder_all$Plasmid
    ),
    "\n",
    sep = ""
)


cat(
    "Sample × plasmid observations: ",
    nrow(plasmid_sample),
    "\n",
    sep = ""
)


cat(
    "CARD samples: ",
    n_card_samples,
    "\n",
    sep = ""
)


cat(
    "ABRicate-CARD ARG records: ",
    nrow(abricate_card_all),
    "\n",
    sep = ""
)


cat(
    "Unique ABRicate-CARD ARGs: ",
    n_distinct(
        abricate_card_all$ARG
    ),
    "\n",
    sep = ""
)


cat(
    "ARG–plasmid sample-level pairs: ",
    nrow(arg_plasmid_cooccurrence),
    "\n",
    sep = ""
)


cat(
    "ARG–plasmid contig-level pairs: ",
    nrow(arg_plasmid_colocalization),
    "\n",
    sep = ""
)


cat(
    "\nOutput directory:\n",
    PLASMID_ANALYSIS_DIR,
    "\n",
    sep = ""
)


sink()


###############################################################################
# 31. COMPLETION SUMMARY
###############################################################################

message("")
message("============================================================")
message("PLASMID ANALYSIS COMPLETED")
message("============================================================")


log_plasmid(
    "PlasmidFinder samples: ",
    n_plasmid_samples
)


log_plasmid(
    "PlasmidFinder hit records: ",
    nrow(plasmidfinder_all)
)


log_plasmid(
    "Unique plasmid replicons: ",
    n_distinct(
        plasmidfinder_all$Plasmid
    )
)


log_plasmid(
    "Samples with plasmid hits: ",
    n_distinct(
        plasmid_sample$SampleID
    )
)


log_plasmid(
    "Samples without plasmid hits: ",
    n_plasmid_samples -
        n_distinct(
            plasmid_sample$SampleID
        )
)


log_plasmid(
    "ABRicate-CARD ARG records: ",
    nrow(abricate_card_all)
)


log_plasmid(
    "Unique ABRicate-CARD ARGs: ",
    n_distinct(
        abricate_card_all$ARG
    )
)


log_plasmid(
    "ARG–plasmid sample co-occurrence pairs: ",
    nrow(arg_plasmid_cooccurrence)
)


log_plasmid(
    "ARG–plasmid contig co-localization pairs: ",
    nrow(arg_plasmid_colocalization)
)


log_plasmid(
    "Results written to: ",
    PLASMID_ANALYSIS_DIR
)


log_plasmid(
    "R object written to: ",
    file.path(
        OBJECT_DIR,
        "plasmid_analysis.rds"
    )
)


message("============================================================")
