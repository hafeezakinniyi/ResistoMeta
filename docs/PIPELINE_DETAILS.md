# ResistoMeta Pipeline Details

This document provides the technical reference for the ResistoMeta workflow. The main [`README.md`](../README.md) contains the concise user-facing introduction.

## 1. Execution model

ResistoMeta uses a centralized configuration and a two-phase workflow:

``` text
config/config.sh
      │
      ▼
validate_installation.sh
      │
      ▼
run_pipeline.sh
      │
      ├── Phase I: Bioinformatics
      │     ├── Modules 01–07
      │     └── Module 06 optional
      │
      └── Phase II: Downstream R analysis
            └── Module 08
                  └── run_analysis.R
                        ├── taxonomy_analysis.R
                        ├── resistome_analysis.R
                        └── plasmid_analysis.R
```

`config/config.sh` is the single configuration source. Module 08 exports the relevant configuration to the R workflow rather than redefining analytical parameters.

## 2. Bioinformatics modules

### Module 01 — Quality Control

`01_quality_control.sh` performs FastQC on raw reads, adapter/quality preprocessing with fastp, FastQC on trimmed reads, and MultiQC reporting.

Usage:

``` bash
bash scripts/bash/modules/01_quality_control.sh SAMPLE_ID
```

Generate project-level summaries:

``` bash
bash scripts/bash/modules/01_quality_control.sh summary
```

### Module 02 — Host Depletion

`02_host_depletion.sh` removes host-derived reads from quality-controlled paired-end reads using Bowtie2 and the configured host-reference index. Host-depleted reads are the principal input for taxonomic profiling, read-based AMR detection, and assembly.

Usage:

``` bash
bash scripts/bash/modules/02_host_depletion.sh SAMPLE_ID
```

Project-level summary:

``` bash
bash scripts/bash/modules/02_host_depletion.sh summary
```

### Module 03 — Taxonomic Profiling

`03_taxonomic_profiling.sh` classifies host-depleted reads with Kraken2 and estimates abundance with Bracken. Profiles are generated at phylum, genus, and species levels and are consumed by the downstream R taxonomy analysis.

Usage:

``` bash
bash scripts/bash/modules/03_taxonomic_profiling.sh SAMPLE_ID
```

### Module 04 — Read-based Resistome

`04_read_based_resistome.sh` detects ARGs directly from host-depleted paired-end reads using the CARD RGI BWT workflow. It also validates paired-end read counts and generates the project-level sequencing-depth table used for downstream resistome normalization.

Usage:

``` bash
bash scripts/bash/modules/04_read_based_resistome.sh SAMPLE_ID
```

Generate the project-level sequencing-depth table:

``` bash
bash scripts/bash/modules/04_read_based_resistome.sh summary
```

For paired-end data with `N` R1 reads and `N` R2 reads, the pipeline defines:

``` text
R1_reads        = N
R2_reads        = N
total_sequences = N
```

R1 and R2 are paired mates representing the same sequencing fragments, so they are not counted as `2N` independent sequences.

### Module 05 — Metagenome Assembly

`05_metagenome_assembly.sh` performs de novo assembly of host-depleted paired-end reads with MEGAHIT. The principal assembly is `final.contigs.fa`.

Usage:

``` bash
bash scripts/bash/modules/05_metagenome_assembly.sh SAMPLE_ID
```

### Module 06 — Assembly Quality

`06_assembly_quality.sh` evaluates MEGAHIT assemblies with MetaQUAST.

Usage:

``` bash
bash scripts/bash/modules/06_assembly_quality.sh SAMPLE_ID
```

**Default:**

``` bash
RUN_ASSEMBLY_QC=false
```

Module 06 is disabled by default because assembly-quality assessment can substantially increase total runtime. It can be enabled when required:

``` bash
RUN_ASSEMBLY_QC=true
```

The module is independent of the downstream R workflow and does not alter the primary AMR analyses.

### Module 07 — Contig Annotation

`07_contig_annotation.sh` analyzes assembled contigs with:

- ABRicate-CARD for contig-level ARG detection;
- ABRicate-PlasmidFinder for plasmid replicon detection.

The module applies the configured ABRicate identity and coverage thresholds. Downstream R modules consume these filtered outputs rather than redefining the primary hit-detection thresholds.

Usage:

``` bash
bash scripts/bash/modules/07_contig_annotation.sh SAMPLE_ID
```

## 3. Module 08 and downstream R analysis

`08_launch_R_analysis.sh` validates the R workflow and the upstream results required by the enabled analyses, activates the configured R environment, exports configuration, and launches `run_analysis.R`.

Usage:

``` bash
bash scripts/bash/modules/08_launch_R_analysis.sh
```

`run_analysis.R` orchestrates three downstream analyses:

### Taxonomy analysis

Consumes Bracken profiles and produces taxonomic summaries and relative-abundance outputs at the configured levels.

### Resistome analysis

Consumes RGI read-based results and ABRicate-CARD contig results. It calculates configured ARG richness, prevalence, abundance, drug-class summaries, and related resistome measures.

``` text
ARG richness = number of unique detected ARGs
ARG abundance = normalized RGI mapped-read signal expressed as RPKM
ARG prevalence = proportion of analyzed samples containing a given ARG
```

RPKM (Read per Kilobase Million) is the configured primary abundance metric for read-based ARG analysis, accounting for both sequencing depth and ARG reference length.

The current RGI filtering configuration specifies:

``` text
Minimum mapped reads: 5
Minimum average coverage: 80%
```

The current ABRicate filtering configuration is:

``` text
Minimum identity: 80%
Minimum coverage: 80%
```

These thresholds are configured centrally and applied during downstream resistome processing.

### Plasmid analysis

Consumes ABRicate-PlasmidFinder and ABRicate-CARD results to characterize plasmid replicons, calculate richness and prevalence, and evaluate ARG–plasmid co-occurrence and contig-level co-localization.

``` text
ARG–plasmid co-occurrence = an ARG and a plasmid replicon are detected in the same biological sample
Contig-level co-localization = an ARG and a plasmid replicon are detected on the same assembled contig.
```

Sample-level co-occurrence does **not** demonstrate physical linkage. Contig-level co-localization provides sequence-level evidence but does not by itself establish complete plasmid structure, circularity, mobility, or biological transfer.

## 4. Configuration

The principal configuration file is:

``` text
config/config.sh
```

Important configuration groups include:

- project metadata;
- computational resources (`THREADS`, `R_THREADS`);
- Conda environment names;
- input, reference, database, script, and result directories;
- sequencing-depth input location;
- RGI thresholds;
- ABRicate identity and coverage thresholds;
- taxonomic analysis levels and top-taxon settings;
- resistome and plasmid calculation switches;
- figure-generation options;
- pipeline validation options;
- module execution switches.

The main module switches include:

``` bash
RUN_QC=true
RUN_HOST_DEPLETION=true
RUN_TAXONOMY=true
RUN_READ_RESISTOME=true
RUN_ASSEMBLY=true
RUN_ASSEMBLY_QC=false
RUN_CONTIG_ANNOTATION=true
RUN_R_ANALYSIS=true
```

Individual downstream R analyses can also be enabled or disabled independently.

## 5. Pipeline modes

The configuration supports three execution modes:

``` text
full
    Execute the complete configured workflow.

bioinformatics
    Execute the bioinformatics modules without Module 08.

analysis
    Execute downstream R analysis using available upstream results.
```

The individual module switches remain authoritative within the selected mode.

## 6. Validation and resumability

The installation validator checks the required software and reference resources used by the configured workflow. It should be run before the main pipeline, especially after creating or modifying Conda environments.

The shared Bash utilities also provide file/directory checks, expected-output validation, paired-read validation, Bowtie2 index validation, environment validation, command execution with error handling, timing, and completed-step detection.

When expected outputs already exist, supported steps can be skipped rather than recomputed. This allows downstream problems or configuration changes to be addressed without automatically repeating all upstream processing.

## 7. Configuration provenance

Module 08 writes the active downstream configuration to:

``` text
results/08_downstream_analysis/analysis_configuration.txt
```

For reproducibility, retain the pipeline version, configuration file, software/environment information, reference database versions or identifiers, analytical thresholds, sequencing-depth definition, and relevant output summaries.
