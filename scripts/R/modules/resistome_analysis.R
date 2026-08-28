###############################################################################
# RESISTOMETA
# Module: Resistome Analysis
#
# Purpose:
#   Performs downstream analysis of antimicrobial resistance genes detected
#   using:
#
#       1. RGI-CARD read-based resistance profiling
#       2. ABRicate-CARD contig-based resistance profiling
#
# SAMPLE UNIVERSE:
#   The sequencing-depth table is the sole source defining the samples
#   included in downstream resistome analysis.
#
#   Required sequencing-depth columns:
#
#       SampleID
#       total_sequences
#
# PREVALENCE:
#   Prevalence is calculated using the number of validated samples in the
#   sequencing-depth table as the denominator.
#
#       prevalence =
#           number of samples with ARG
#           --------------------------
#           number of validated samples
#
# RGI-CARD:
#       - high-confidence ARG detection
#       - ARG richness
#       - ARG abundance (RPKM)
#       - ARG prevalence
#       - drug-class ARG counts
#       - drug-class abundance
#       - drug-class prevalence
#
# ABRicate-CARD:
#       - high-confidence contig-level ARG detection
#       - within-sample gene deduplication
#       - ARG richness
#       - ARG prevalence
#
# IMPORTANT:
#
#   ARG richness and ARG abundance are distinct measures.
#
#   Richness:
#       Number of unique ARGs detected.
#
#   Abundance:
#       Normalized RGI mapped-read signal expressed as RPKM.
#
#   Prevalence:
#       Proportion of validated sequenced samples in which an ARG is detected.
#
#   Sequencing depth:
#       Used both to define the sample universe and to normalize RGI
#       read-based abundance.
#
# CONFIGURATION ARCHITECTURE:
#
#   config.sh
#       |
#       v
#   08_launch_R_analysis.sh
#       |
#       v
#   resistome_analysis.R
#
#   This module consumes exported configuration. It does not redefine
#   analytical thresholds or calculation switches.
#
###############################################################################


###############################################################################
# 1. LOAD REQUIRED PACKAGES
###############################################################################

suppressPackageStartupMessages({

    library(dplyr)
    library(tidyr)
    library(stringr)
    library(readr)
    library(purrr)
    library(tibble)

})


###############################################################################
# 2. ENVIRONMENT HELPERS
###############################################################################

get_env <- function(name, required = TRUE, default = NULL) {

    value <- Sys.getenv(
        name,
        unset = ""
    )

    if (identical(value, "")) {

        if (required) {

            stop(
                paste0(
                    "Required environment variable is not set: ",
                    name
                )
            )

        }

        return(default)

    }

    value

}


as_bool <- function(x, name) {

    if (
        length(x) != 1 ||
        is.na(x) ||
        identical(trimws(x), "")
    ) {

        stop(
            paste0(
                "Missing Boolean configuration value: ",
                name
            )
        )

    }

    value <- tolower(trimws(x))

    valid_values <- c(
        "true", "false",
        "t", "f",
        "1", "0",
        "yes", "no",
        "y", "n"
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


###############################################################################
# 3. READ ENVIRONMENT CONFIGURATION
###############################################################################


###############################################################################
# Sequencing-depth input
###############################################################################

SEQUENCING_DEPTH_FILE <- get_env(
    "METAGENOMIC_AMR_SEQUENCING_DEPTH_FILE"
)


###############################################################################
# Upstream result directories
###############################################################################

READ_RESISTOME_DIR <- get_env(
    "METAGENOMIC_AMR_READ_RESISTOME_DIR"
)


CONTIG_ANNOTATION_DIR <- get_env(
    "METAGENOMIC_AMR_CONTIG_ANNOTATION_DIR"
)


###############################################################################
# Output directories
###############################################################################

RESISTOME_ANALYSIS_DIR <- get_env(
    "METAGENOMIC_AMR_RESISTOME_ANALYSIS_DIR"
)


FIGURE_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_FIGURE_DIR"
)


OBJECT_DIR <- get_env(
    "METAGENOMIC_AMR_ANALYSIS_OBJECT_DIR"
)


###############################################################################
# Analysis parameters
###############################################################################

RGI_MIN_MAPPED_READS <- as.numeric(
    get_env(
        "METAGENOMIC_AMR_RGI_MIN_MAPPED_READS"
    )
)


RGI_MIN_COVERAGE <- as.numeric(
    get_env(
        "METAGENOMIC_AMR_RGI_MIN_COVERAGE"
    )
)


###############################################################################
# Resistome abundance metric
#
# This value is defined by config.sh and exported by Module 08.
###############################################################################

RESISTOME_ABUNDANCE_METRIC <- get_env(
    "METAGENOMIC_AMR_RESISTOME_ABUNDANCE_METRIC"
)


###############################################################################
# ABRicate-CARD configuration
#
# Database identity and thresholds are configuration values supplied by
# config.sh and exported by Module 08.
###############################################################################

ABRICATE_CARD_DB <- get_env(
    "METAGENOMIC_AMR_ABRICATE_CARD_DB"
)


ABRICATE_MIN_IDENTITY <- as.numeric(
    get_env(
        "METAGENOMIC_AMR_ABRICATE_MIN_IDENTITY"
    )
)


ABRICATE_MIN_COVERAGE <- as.numeric(
    get_env(
        "METAGENOMIC_AMR_ABRICATE_MIN_COVERAGE"
    )
)


###############################################################################
# Analysis switches
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
# 4. VALIDATE ANALYSIS PARAMETERS
###############################################################################

if (
    !identical(
        toupper(RESISTOME_ABUNDANCE_METRIC),
        "RPKM"
    )
) {

    stop(
        paste0(
            "Unsupported resistome abundance metric: ",
            RESISTOME_ABUNDANCE_METRIC,
            ". This module currently requires RPKM."
        )
    )

}


if (
    is.na(RGI_MIN_MAPPED_READS) ||
    !is.finite(RGI_MIN_MAPPED_READS) ||
    RGI_MIN_MAPPED_READS < 0
) {

    stop(
        "RGI minimum mapped reads must be a finite non-negative number."
    )

}


if (
    is.na(RGI_MIN_COVERAGE) ||
    !is.finite(RGI_MIN_COVERAGE) ||
    RGI_MIN_COVERAGE < 0 ||
    RGI_MIN_COVERAGE > 100
) {

    stop(
        "RGI minimum coverage must be a finite number between 0 and 100."
    )

}


if (
    is.na(ABRICATE_MIN_IDENTITY) ||
    !is.finite(ABRICATE_MIN_IDENTITY) ||
    ABRICATE_MIN_IDENTITY < 0 ||
    ABRICATE_MIN_IDENTITY > 100
) {

    stop(
        "ABRicate minimum identity must be a finite number between 0 and 100."
    )

}


if (
    is.na(ABRICATE_MIN_COVERAGE) ||
    !is.finite(ABRICATE_MIN_COVERAGE) ||
    ABRICATE_MIN_COVERAGE < 0 ||
    ABRICATE_MIN_COVERAGE > 100
) {

    stop(
        "ABRicate minimum coverage must be a finite number between 0 and 100."
    )

}


###############################################################################
# 5. CREATE OUTPUT DIRECTORIES
###############################################################################

dir.create(
    RESISTOME_ANALYSIS_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    FIGURE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    OBJECT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


###############################################################################
# 6. LOGGING
###############################################################################

message("============================================================")
message("RESISTOME ANALYSIS")
message("============================================================")

message(
    "Read-based resistome directory: ",
    READ_RESISTOME_DIR
)

message(
    "Contig annotation directory: ",
    CONTIG_ANNOTATION_DIR
)

message(
    "Sequencing-depth file: ",
    SEQUENCING_DEPTH_FILE
)

message(
    "Output directory: ",
    RESISTOME_ANALYSIS_DIR
)

message(
    "RGI threshold: ",
    RGI_MIN_MAPPED_READS,
    " mapped reads / ",
    RGI_MIN_COVERAGE,
    "% coverage"
)

message(
    "ABRicate-CARD database: ",
    ABRICATE_CARD_DB
)

message(
    "ABRicate threshold: ",
    ABRICATE_MIN_IDENTITY,
    "% identity / ",
    ABRICATE_MIN_COVERAGE,
    "% coverage"
)

message(
    "Abundance metric: ",
    RESISTOME_ABUNDANCE_METRIC
)


###############################################################################
# 7. VALIDATE INPUTS
###############################################################################

message("")
message("Validating resistome analysis inputs...")


if (!file.exists(SEQUENCING_DEPTH_FILE)) {

    stop(
        paste0(
            "Sequencing-depth file not found:\n",
            SEQUENCING_DEPTH_FILE,
            "\n\n",
            "The sequencing-depth table is required because it defines ",
            "the analysis sample universe and is required for RPKM."
        )
    )

}


if (!dir.exists(READ_RESISTOME_DIR)) {

    stop(
        paste0(
            "Read-based resistome directory not found:\n",
            READ_RESISTOME_DIR
        )
    )

}


if (!dir.exists(CONTIG_ANNOTATION_DIR)) {

    warning(
        paste0(
            "Contig annotation directory not found:\n",
            CONTIG_ANNOTATION_DIR,
            "\n",
            "ABRicate-CARD analyses will be unavailable."
        )
    )

}


###############################################################################
# 8. READ AND VALIDATE SEQUENCING DEPTH
###############################################################################

message("")
message("Reading sequencing-depth information...")


sequencing_depth <- read_csv(
    SEQUENCING_DEPTH_FILE,
    show_col_types = FALSE
)


required_depth_columns <- c(
    "SampleID",
    "R1_reads",
    "R2_reads",
    "total_sequences"
)


missing_depth_columns <- setdiff(
    required_depth_columns,
    names(sequencing_depth)
)


if (length(missing_depth_columns) > 0) {

    stop(
        paste0(
            "Sequencing-depth file is missing required columns:\n",
            paste(
                missing_depth_columns,
                collapse = ", "
            )
        )
    )

}


sequencing_depth <- sequencing_depth %>%

    select(
        SampleID,
        R1_reads,
        R2_reads,
        total_sequences
    ) %>%

    mutate(

        SampleID =
            as.character(SampleID),

        R1_reads =
            as.numeric(R1_reads),

        R2_reads =
            as.numeric(R2_reads),

        total_sequences =
            as.numeric(total_sequences)

    )


###############################################################################
# 8A. VALIDATE SAMPLE IDs
###############################################################################

if (
    any(
        is.na(sequencing_depth$SampleID) |
        sequencing_depth$SampleID == ""
    )
) {

    stop(
        "Sequencing-depth table contains missing or empty SampleID values."
    )

}


duplicate_depth_samples <- sequencing_depth %>%

    count(SampleID) %>%

    filter(n > 1)


if (nrow(duplicate_depth_samples) > 0) {

    stop(
        paste0(
            "Sequencing-depth table contains duplicate SampleID values:\n",
            paste(
                duplicate_depth_samples$SampleID,
                collapse = ", "
            )
        )
    )

}


###############################################################################
# 8B. VALIDATE SEQUENCING DEPTH
###############################################################################

invalid_depth <- sequencing_depth %>%

    filter(
        is.na(total_sequences) |
        total_sequences <= 0
    )


if (nrow(invalid_depth) > 0) {

    stop(
        paste0(
            "Invalid sequencing depth detected for sample(s):\n",
            paste(
                invalid_depth$SampleID,
                collapse = ", "
            ),
            "\n\n",
            "'total_sequences' must be a positive numeric value."
        )
    )

}


###############################################################################
# 8C. DEFINE ANALYSIS SAMPLE UNIVERSE
###############################################################################
#
# The sequencing-depth table defines the complete sample universe.
#
###############################################################################

all_samples <- sequencing_depth$SampleID


n_total_samples <- length(all_samples)


if (n_total_samples == 0) {

    stop(
        "No valid samples were found in the sequencing-depth table."
    )

}


message(
    "Validated samples in analysis universe: ",
    n_total_samples
)

message(
    "Sequencing depth source: Module 04 host-depleted RGI input reads"
)

###############################################################################
# 9. LOCATE RGI OUTPUT FILES
###############################################################################

message("")
message("Searching for RGI gene-mapping files...")


rgi_files <- list.files(
    READ_RESISTOME_DIR,
    pattern = "gene_mapping_data.*\\.txt$",
    recursive = TRUE,
    full.names = TRUE
)


if (length(rgi_files) == 0) {

    stop(
        paste0(
            "No RGI gene_mapping_data files were found in:\n",
            READ_RESISTOME_DIR
        )
    )

}


message(
    "RGI files found: ",
    length(rgi_files)
)


###############################################################################
# 10. IMPORT RGI RESULTS
###############################################################################

read_rgi_file <- function(file) {

    message(
        "Reading RGI: ",
        basename(file)
    )

    df <- read_tsv(
        file,
        show_col_types = FALSE,
        progress = FALSE
    )


    required_columns <- c(
        "ARO Term",
        "All Mapped Reads",
        "Average Percent Coverage",
        "Reference Length",
        "AMR Gene Family",
        "Drug Class",
        "Resistance Mechanism"
    )


    missing_columns <- setdiff(
        required_columns,
        names(df)
    )


    if (length(missing_columns) > 0) {

        stop(
            paste0(
                "Missing required RGI columns in ",
                basename(file),
                ": ",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )

    }


 ###########################################################################
    # Extract SampleID from RGI filename
    #
    # Expected filename:
    #
    #   SAMPLE_ID.gene_mapping_data.txt
 ###########################################################################

    sample_id <- str_remove(
        basename(file),
        "\\.gene_mapping_data.*$"
    )


    if (
        is.na(sample_id) ||
        sample_id == ""
    ) {

        stop(
            paste0(
                "Could not extract SampleID from RGI filename: ",
                basename(file)
            )
        )

    }


    df %>%

        mutate(
            SampleID = sample_id
        )

}

    rgi_all <- map_dfr(
	rgi_files,
	read_rgi_file
    )

    message(
        "RGI samples detected: ",
        paste(
            sort(unique(rgi_all$SampleID)),
            collapse = ", "
        )
    )

###############################################################################
# 11. CLEAN RGI DATA
###############################################################################

rgi_all <- rgi_all %>%

    mutate(

        `All Mapped Reads` =
            as.numeric(`All Mapped Reads`),

        `Average Percent Coverage` =
            as.numeric(`Average Percent Coverage`),

        `Reference Length` =
            as.numeric(`Reference Length`)

    )


invalid_rgi_numeric <- rgi_all %>%

    filter(

        is.na(`All Mapped Reads`) |
        !is.finite(`All Mapped Reads`) |
        `All Mapped Reads` < 0 |

        is.na(`Average Percent Coverage`) |
        !is.finite(`Average Percent Coverage`) |
        `Average Percent Coverage` < 0 |
        `Average Percent Coverage` > 100 |

        (
            !is.na(`Reference Length`) &
            (
                !is.finite(`Reference Length`) |
                `Reference Length` <= 0
            )
        )

    )


if (nrow(invalid_rgi_numeric) > 0) {

    stop(
        paste0(
            "Invalid numeric values detected in RGI gene-mapping results. ",
            "Mapped reads and coverage must be valid non-negative values; ",
            "coverage must be <=100; reference length must be positive when present."
        )
    )

}


###############################################################################
# 12. VALIDATE RGI SAMPLE IDs
###############################################################################

rgi_samples <- unique(
    rgi_all$SampleID
)


###############################################################################
# 12A. RGI samples not represented in sequencing-depth table
###############################################################################

rgi_without_depth <- setdiff(
    rgi_samples,
    all_samples
)


if (length(rgi_without_depth) > 0) {

    stop(
        paste0(
            "RGI results contain samples not present in the ",
            "sequencing-depth table:\n",
            paste(
                rgi_without_depth,
                collapse = ", "
            ),
            "\n\n",
            "Unexpected sample IDs are not permitted."
        )
    )

}


###############################################################################
# 12B. Sequencing-depth samples without RGI results
###############################################################################

depth_without_rgi <- setdiff(
    all_samples,
    rgi_samples
)


if (length(depth_without_rgi) > 0) {

    warning(
        paste0(
            "The following sequencing-depth samples have no RGI results:\n",
            paste(
                depth_without_rgi,
                collapse = ", "
            ),
            "\n\n",
            "These samples remain part of the analysis universe and will ",
            "receive zero RGI ARG richness/abundance where applicable."
        )
    )

}


###############################################################################
# 12C. Report matched samples
###############################################################################

matched_rgi_samples <- intersect(
    rgi_samples,
    all_samples
)


message(
    "RGI samples matched to sequencing-depth table: ",
    length(matched_rgi_samples),
    " / ",
    n_total_samples
)


###############################################################################
# 13. FILTER HIGH-CONFIDENCE RGI HITS
###############################################################################

message("")
message("Applying RGI high-confidence filters...")


rgi_filtered <- rgi_all %>%

    filter(

        `All Mapped Reads` >=
            RGI_MIN_MAPPED_READS,

        `Average Percent Coverage` >=
            RGI_MIN_COVERAGE

    )


message(
    "RGI records before filtering: ",
    nrow(rgi_all)
)

message(
    "RGI records after filtering: ",
    nrow(rgi_filtered)
)


###############################################################################
# 14. SAVE FILTERED RGI RESULTS
###############################################################################

write_csv(
    rgi_filtered,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "rgi_high_confidence_hits.csv"
    )
)


###############################################################################
# 15. CREATE SAMPLE × ARG TABLE
###############################################################################

rgi_gene_sample <- rgi_filtered %>%

    group_by(
        SampleID,
        `ARO Term`
    ) %>%

    summarise(

        total_mapped_reads =
            sum(
                `All Mapped Reads`,
                na.rm = TRUE
            ),

        mean_coverage =
            mean(
                `Average Percent Coverage`,
                na.rm = TRUE
            ),

        max_coverage =
            max(
                `Average Percent Coverage`,
                na.rm = TRUE
            ),

        gene_length =
            first(
                `Reference Length`
            ),

        AMR_Gene_Family =
            first(
                `AMR Gene Family`
            ),

        Drug_Class =
            first(
                `Drug Class`
            ),

        Resistance_Mechanism =
            first(
                `Resistance Mechanism`
            ),

        .groups = "drop"

    )


###############################################################################
# 16. HANDLE INVALID GENE LENGTHS
###############################################################################

if (
    any(
        is.na(rgi_gene_sample$gene_length) |
        rgi_gene_sample$gene_length <= 0
    )
) {

    warning(
        paste0(
            "Some RGI hits have missing or invalid reference lengths. ",
            "RPKM will be unavailable for those records."
        )
    )

}


###############################################################################
# 17. JOIN SEQUENCING DEPTH
###############################################################################

message("")
message("Joining sequencing depth to RGI results...")


rgi_gene_sample <- rgi_gene_sample %>%

    left_join(
        sequencing_depth,
        by = "SampleID"
    )


###############################################################################
# 18. VALIDATE RGI SAMPLE DEPTH
###############################################################################

missing_rgi_depth <- unique(
    rgi_gene_sample$SampleID[
        is.na(rgi_gene_sample$total_sequences)
    ]
)


if (length(missing_rgi_depth) > 0) {

    stop(
        paste0(
            "RGI results contain samples without sequencing-depth information:\n",
            paste(
                missing_rgi_depth,
                collapse = ", "
            )
        )
    )

}


###############################################################################
# 19. CALCULATE RPKM
###############################################################################
#
# RPKM =
#
#       mapped reads × 1,000,000,000
#       --------------------------------
#       gene length (bp) × total reads
#
###############################################################################

message("")
message("Calculating RGI abundance as RPKM...")


rgi_gene_sample <- rgi_gene_sample %>%

    mutate(

        gene_length_kb =
            gene_length / 1000,

        total_sequences_million =
            total_sequences / 1e6,

        RPKM =
            if_else(

                !is.na(gene_length_kb) &
                gene_length_kb > 0 &
                !is.na(total_sequences_million) &
                total_sequences_million > 0,

                total_mapped_reads /
                (
                    gene_length_kb *
                    total_sequences_million
                ),

                NA_real_

            )

    )


###############################################################################
# 20. SAVE SAMPLE × ARG TABLE
###############################################################################

write_csv(
    rgi_gene_sample,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "rgi_arg_sample_table.csv"
    )
)


###############################################################################
# 21. ARG RICHNESS PER SAMPLE
###############################################################################

if (CALCULATE_ARG_RICHNESS) {

    message("")
    message("Calculating RGI ARG richness...")


    rgi_richness <- rgi_gene_sample %>%

        group_by(SampleID) %>%

        summarise(

            ARG_Richness =
                n_distinct(`ARO Term`),

            .groups = "drop"

        ) %>%

        right_join(
            sequencing_depth %>%
                select(SampleID),
            by = "SampleID"
        ) %>%

        mutate(

            ARG_Richness =
                replace_na(
                    ARG_Richness,
                    0
                )

        )


    write_csv(
        rgi_richness,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "rgi_arg_richness.csv"
        )
    )

}


###############################################################################
# 22. SAMPLE-LEVEL RGI RESISTOME ABUNDANCE
###############################################################################
#
# ARG_RPKM:
#   Sum of gene-level RPKM values detected within a sample.
#
# This represents aggregate normalized resistome signal.
#
###############################################################################

rgi_resistome_load <- rgi_gene_sample %>%

    group_by(SampleID) %>%

    summarise(

        total_mapped_reads =
            sum(
                total_mapped_reads,
                na.rm = TRUE
            ),

        ARG_RPKM =
            sum(
                RPKM,
                na.rm = TRUE
            ),

        n_detected_ARGs =
            n_distinct(`ARO Term`),

        .groups = "drop"

    ) %>%

    right_join(
        sequencing_depth,
        by = "SampleID"
    ) %>%

    mutate(

        total_mapped_reads =
            replace_na(
                total_mapped_reads,
                0
            ),

        ARG_RPKM =
            replace_na(
                ARG_RPKM,
                0
            ),

        n_detected_ARGs =
            replace_na(
                n_detected_ARGs,
                0
            )

    ) %>%

    select(
        SampleID,
        total_sequences,
        total_mapped_reads,
        n_detected_ARGs,
        ARG_RPKM
    )


write_csv(
    rgi_resistome_load,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "rgi_resistome_load.csv"
    )
)


###############################################################################
# 23. ARG PREVALENCE
###############################################################################

if (CALCULATE_ARG_PREVALENCE) {

    message("")
    message("Calculating ARG prevalence...")


    rgi_prevalence <- rgi_gene_sample %>%

        group_by(`ARO Term`) %>%

        summarise(

            n_samples =
                n_distinct(SampleID),

            .groups = "drop"

        ) %>%

        mutate(

            prevalence =
                n_samples /
                n_total_samples,

            prevalence_percent =
                prevalence * 100

        ) %>%

        arrange(
            desc(prevalence),
            desc(n_samples)
        )


    write_csv(
        rgi_prevalence,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "rgi_arg_prevalence.csv"
        )
    )

}


###############################################################################
# 24. DRUG-CLASS STANDARDIZATION
###############################################################################

normalize_drug_class <- function(x) {

    x <- as.character(x)

    x <- str_squish(x)

    x <- str_replace_all(
        x,
        "\\s*;\\s*",
        ";"
    )

    x

}


###############################################################################
# 25. CREATE ARG × DRUG-CLASS TABLE
###############################################################################

rgi_gene_classes <- rgi_gene_sample %>%

    select(
        SampleID,
        `ARO Term`,
        Drug_Class
    ) %>%

    distinct() %>%

    separate_rows(
        Drug_Class,
        sep = ";"
    ) %>%

    mutate(

        Drug_Class =
            normalize_drug_class(
                Drug_Class
            )

    ) %>%

    filter(

        !is.na(Drug_Class),

        Drug_Class != ""

    )


###############################################################################
# 26. SAVE DRUG-CLASS MAPPING
###############################################################################

write_csv(
    rgi_gene_classes,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "rgi_arg_drug_class_mapping.csv"
    )
)


###############################################################################
# 27. DRUG-CLASS ARG COUNTS
###############################################################################

if (CALCULATE_ARG_CLASS_RICHNESS) {

    message("")
    message("Calculating drug-class ARG counts...")


    rgi_class_richness <- rgi_gene_classes %>%

        distinct(
            `ARO Term`,
            Drug_Class
        ) %>%

        count(
            Drug_Class,
            name = "n_genes"
        ) %>%

        arrange(
            desc(n_genes)
        )


    write_csv(
        rgi_class_richness,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "rgi_drug_class_arg_counts.csv"
        )
    )

}


###############################################################################
# 28. DRUG-CLASS ABUNDANCE
###############################################################################
#
# An ARG may belong to multiple CARD drug classes.
#
# Therefore, class-specific RPKM values may overlap and should not be
# interpreted as mutually exclusive components of a 100% composition.
#
###############################################################################

if (CALCULATE_ARG_CLASS_ABUNDANCE) {

    message("")
    message("Calculating drug-class abundance...")


    rgi_class_abundance <- rgi_gene_sample %>%

        select(
            SampleID,
            `ARO Term`,
            RPKM
        ) %>%

        inner_join(
            rgi_gene_classes,
            by = c(
                "SampleID",
                "ARO Term"
            )
        ) %>%

        group_by(
            SampleID,
            Drug_Class
        ) %>%

        summarise(

            RPKM =
                sum(
                    RPKM,
                    na.rm = TRUE
                ),

            n_ARGs =
                n_distinct(
                    `ARO Term`
                ),

            .groups = "drop"

        )


    write_csv(
        rgi_class_abundance,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "rgi_drug_class_abundance.csv"
        )
    )


    ###########################################################################
    # Global drug-class summary
    ###########################################################################

    rgi_class_global <- rgi_class_abundance %>%

        group_by(
            Drug_Class
        ) %>%

        summarise(

            mean_RPKM =
                mean(
                    RPKM,
                    na.rm = TRUE
                ),

            median_RPKM =
                median(
                    RPKM,
                    na.rm = TRUE
                ),

            total_RPKM =
                sum(
                    RPKM,
                    na.rm = TRUE
                ),

            samples_with_class =
                sum(
                    RPKM > 0
                ),

            prevalence =
                samples_with_class /
                n_total_samples,

            prevalence_percent =
                prevalence * 100,

            .groups = "drop"

        ) %>%

        arrange(
            desc(mean_RPKM)
        )


    write_csv(
        rgi_class_global,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "rgi_drug_class_global_summary.csv"
        )
    )

}


###############################################################################
# 28B. DRUG-CLASS PREVALENCE
###############################################################################
#
# Prevalence uses the complete sequencing-depth sample universe as the
# denominator. An ARG may map to multiple drug classes.
#
###############################################################################

rgi_class_prevalence <- rgi_gene_classes %>%

    distinct(
        SampleID,
        Drug_Class
    ) %>%

    count(
        Drug_Class,
        name = "n_samples"
    ) %>%

    mutate(

        prevalence =
            n_samples /
            n_total_samples,

        prevalence_percent =
            prevalence * 100

    ) %>%

    arrange(
        desc(prevalence),
        desc(n_samples),
        Drug_Class
    )


write_csv(
    rgi_class_prevalence,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "rgi_drug_class_prevalence.csv"
    )
)


###############################################################################
# 29. LOCATE ABRICATE-CARD FILES
###############################################################################

message("")
message("Searching for ABRicate-CARD results...")


abricate_card_dir <- file.path(
    CONTIG_ANNOTATION_DIR,
    "abricate_card"
)


if (dir.exists(abricate_card_dir)) {

    abricate_card_files <- list.files(

        abricate_card_dir,

        pattern = "\\.tsv$",

        recursive = TRUE,

        full.names = TRUE

    )

} else {

    abricate_card_files <- character(0)

}


if (length(abricate_card_files) == 0) {

    warning(
        paste0(
            "No ABRicate-CARD TSV files found in:\n",
            abricate_card_dir
        )
    )

} else {

    message(
        "ABRicate-CARD files found: ",
        length(abricate_card_files)
    )

}


###############################################################################
# 30. READ ABRICATE-CARD
###############################################################################

read_abricate_card <- function(file) {

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
        "SEQUENCE",
        "GENE",
        "%COVERAGE",
        "%IDENTITY"
    )


    missing_columns <- setdiff(
        required_columns,
        names(df)
    )


    if (length(missing_columns) > 0) {

        stop(
            paste0(
                "Missing required ABRicate-CARD columns in ",
                basename(file),
                ": ",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )

    }


    sample_id <- tools::file_path_sans_ext(
        basename(file)
    )


    df %>%

        mutate(

            SampleID =
                str_remove(
                    sample_id,
                    "_card$"
                ),

            `%COVERAGE` =
                as.numeric(`%COVERAGE`),

            `%IDENTITY` =
                as.numeric(`%IDENTITY`)

        )

}


if (length(abricate_card_files) > 0) {

    abricate_card_all <- map_dfr(
        abricate_card_files,
        read_abricate_card
    )

} else {

    abricate_card_all <- tibble()

}


###############################################################################
# 31. VALIDATE ABRICATE SAMPLE IDs
###############################################################################

if (nrow(abricate_card_all) > 0) {

    abricate_samples <- unique(
        abricate_card_all$SampleID
    )


    abricate_without_depth <- setdiff(
        abricate_samples,
        all_samples
    )


    if (length(abricate_without_depth) > 0) {

        stop(
            paste0(
                "ABRicate-CARD results contain samples not present in the ",
                "sequencing-depth table:\n",
                paste(
                    abricate_without_depth,
                    collapse = ", "
                ),
                "\n\n",
                "Unexpected sample IDs are not permitted."
            )
        )

    }

}


###############################################################################
# 32. APPLY CONFIGURED ABRICATE THRESHOLD
###############################################################################

if (nrow(abricate_card_all) > 0) {

    abricate_card_filtered <- abricate_card_all %>%

        filter(

            `%IDENTITY` >=
                ABRICATE_MIN_IDENTITY,

            `%COVERAGE` >=
                ABRICATE_MIN_COVERAGE

        )

} else {

    abricate_card_filtered <- tibble()

}


###############################################################################
# 33. SAVE FILTERED ABRICATE RESULTS
###############################################################################

if (nrow(abricate_card_filtered) > 0) {

    write_csv(
        abricate_card_filtered,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "abricate_card_high_confidence_hits.csv"
        )
    )

}


###############################################################################
# 34. DEDUPLICATE ABRICATE-CARD WITHIN SAMPLE
###############################################################################

if (nrow(abricate_card_filtered) > 0) {

    abricate_gene_sample <- abricate_card_filtered %>%

        group_by(
            SampleID,
            GENE
        ) %>%

        summarise(

            n_contigs =
                n_distinct(SEQUENCE),

            max_identity =
                max(
                    `%IDENTITY`,
                    na.rm = TRUE
                ),

            max_coverage =
                max(
                    `%COVERAGE`,
                    na.rm = TRUE
                ),

            .groups = "drop"

        )

} else {

    abricate_gene_sample <- tibble()

}


###############################################################################
# 35. ABRICATE ARG RICHNESS
###############################################################################

if (
    CALCULATE_ABRICATE_RICHNESS &&
    nrow(abricate_gene_sample) > 0
) {

    message("")
    message("Calculating ABRicate-CARD ARG richness...")


    abricate_richness <- abricate_gene_sample %>%

        group_by(SampleID) %>%

        summarise(

            ABRicate_ARG_Richness =
                n_distinct(GENE),

            .groups = "drop"

        ) %>%

        right_join(
            sequencing_depth %>%
                select(SampleID),
            by = "SampleID"
        ) %>%

        mutate(

            ABRicate_ARG_Richness =
                replace_na(
                    ABRicate_ARG_Richness,
                    0
                )

        )


    write_csv(
        abricate_richness,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "abricate_arg_richness.csv"
        )
    )

}


###############################################################################
# 36. ABRICATE ARG PREVALENCE
###############################################################################

if (
    CALCULATE_ABRICATE_PREVALENCE &&
    nrow(abricate_gene_sample) > 0
) {

    message("")
    message("Calculating ABRicate-CARD ARG prevalence...")


    abricate_prevalence <- abricate_gene_sample %>%

        group_by(GENE) %>%

        summarise(

            n_samples =
                n_distinct(SampleID),

            .groups = "drop"

        ) %>%

        mutate(

            prevalence =
                n_samples /
                n_total_samples,

            prevalence_percent =
                prevalence * 100

        ) %>%

        arrange(
            desc(prevalence),
            desc(n_samples)
        )


    write_csv(
        abricate_prevalence,
        file.path(
            RESISTOME_ANALYSIS_DIR,
            "abricate_arg_prevalence.csv"
        )
    )

}


###############################################################################
# 37. COMBINED RESISTOME SUMMARY
###############################################################################

message("")
message("Creating integrated resistome summary...")


###############################################################################
# Start with sequencing-depth sample universe
###############################################################################

resistome_summary <- sequencing_depth %>%

    select(
        SampleID,
        total_sequences
    )


###############################################################################
# Add RGI richness
###############################################################################

if (exists("rgi_richness")) {

    resistome_summary <- resistome_summary %>%

        left_join(
            rgi_richness,
            by = "SampleID"
        )

}


###############################################################################
# Add RGI resistome abundance
###############################################################################

if (exists("rgi_resistome_load")) {

    resistome_summary <- resistome_summary %>%

        left_join(
            rgi_resistome_load %>%
                select(
                    SampleID,
                    total_mapped_reads,
                    n_detected_ARGs,
                    ARG_RPKM
                ),
            by = "SampleID"
        )

}


###############################################################################
# Add ABRicate richness
###############################################################################

if (exists("abricate_richness")) {

    resistome_summary <- resistome_summary %>%

        left_join(
            abricate_richness,
            by = "SampleID"
        )

}


###############################################################################
# Replace missing values with zero
###############################################################################

# Ensure the integrated sample table has stable columns even when an
# individual calculation switch is disabled.
if (!"ARG_Richness" %in% names(resistome_summary)) {

    resistome_summary$ARG_Richness <- 0

}


if (!"ABRicate_ARG_Richness" %in% names(resistome_summary)) {

    resistome_summary$ABRicate_ARG_Richness <- 0

}


if (!"total_mapped_reads" %in% names(resistome_summary)) {

    resistome_summary$total_mapped_reads <- 0

}


if (!"n_detected_ARGs" %in% names(resistome_summary)) {

    resistome_summary$n_detected_ARGs <- 0

}


if (!"ARG_RPKM" %in% names(resistome_summary)) {

    resistome_summary$ARG_RPKM <- 0

}


resistome_summary <- resistome_summary %>%

    mutate(

        ARG_Richness =
            replace_na(
                ARG_Richness,
                0
            ),

        ABRicate_ARG_Richness =
            replace_na(
                ABRicate_ARG_Richness,
                0
            ),

        total_mapped_reads =
            replace_na(
                total_mapped_reads,
                0
            ),

        n_detected_ARGs =
            replace_na(
                n_detected_ARGs,
                0
            ),

        ARG_RPKM =
            replace_na(
                ARG_RPKM,
                0
            )

    ) %>%

    select(
        SampleID,
        total_sequences,
        ARG_Richness,
        ABRicate_ARG_Richness,
        total_mapped_reads,
        n_detected_ARGs,
        ARG_RPKM
    )


write_csv(
    resistome_summary,
    file.path(
        RESISTOME_ANALYSIS_DIR,
        "resistome_sample_summary.csv"
    )
)


###############################################################################
# 38. SAVE R OBJECT
###############################################################################

saveRDS(

    list(

        sequencing_depth =
            sequencing_depth,

        rgi_all =
            rgi_all,

        rgi_filtered =
            rgi_filtered,

        rgi_gene_sample =
            rgi_gene_sample,

        rgi_gene_classes =
            rgi_gene_classes,

        rgi_richness =
            if (exists("rgi_richness"))
                rgi_richness
            else
                NULL,

        rgi_resistome_load =
            rgi_resistome_load,

        rgi_class_richness =
            if (exists("rgi_class_richness"))
                rgi_class_richness
            else
                NULL,

        rgi_class_abundance =
            if (exists("rgi_class_abundance"))
                rgi_class_abundance
            else
                NULL,

        rgi_class_global =
            if (exists("rgi_class_global"))
                rgi_class_global
            else
                NULL,

        rgi_class_prevalence =
            rgi_class_prevalence,

        resistome_summary =
            resistome_summary,

        configuration =
            list(
                abundance_metric =
                    RESISTOME_ABUNDANCE_METRIC,

                rgi_min_mapped_reads =
                    RGI_MIN_MAPPED_READS,

                rgi_min_coverage =
                    RGI_MIN_COVERAGE,

                abricate_card_db =
                    ABRICATE_CARD_DB,

                abricate_min_identity =
                    ABRICATE_MIN_IDENTITY,

                abricate_min_coverage =
                    ABRICATE_MIN_COVERAGE
            ),

        abricate_card_all =
            abricate_card_all,

        abricate_card_filtered =
            abricate_card_filtered,

        abricate_gene_sample =
            abricate_gene_sample

    ),

    file.path(
        OBJECT_DIR,
        "resistome_analysis.rds"
    )

)


###############################################################################
# 39. COMPLETION
###############################################################################

message("")
message("============================================================")
message("RESISTOME ANALYSIS COMPLETED")
message("============================================================")


message(
    "Validated analysis samples: ",
    n_total_samples
)


message(
    "RGI high-confidence ARG records: ",
    nrow(rgi_filtered)
)


message(
    "RGI unique ARGs: ",
    n_distinct(
        rgi_filtered$`ARO Term`
    )
)


if (nrow(abricate_card_filtered) > 0) {

    message(
        "ABRicate-CARD high-confidence hits: ",
        nrow(abricate_card_filtered)
    )

    message(
        "ABRicate-CARD unique ARGs: ",
        n_distinct(
            abricate_card_filtered$GENE
        )
    )

}


message(
    "Prevalence denominator: ",
    n_total_samples,
    " validated sequenced samples"
)


message(
    "ABRicate threshold: ",
    ABRICATE_MIN_IDENTITY,
    "% identity / ",
    ABRICATE_MIN_COVERAGE,
    "% coverage"
)


message(
    "Results written to: ",
    RESISTOME_ANALYSIS_DIR
)


message("============================================================")
