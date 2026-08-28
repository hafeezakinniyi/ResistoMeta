# ResistoMeta Output Reference

This document describes the principal output organization produced by ResistoMeta.

## Output directory

All pipeline results are written beneath:

``` text
results/
```

The processing stages are organized as follows:

``` text
results/
├── 01_quality_control/
├── 02_host_depletion/
├── 03_taxonomy/
├── 04_read_based_resistome/
├── 05_assembly/
├── 06_assembly_quality/
├── 07_contig_annotation/
└── 08_downstream_analysis/
```

## Module 01 — Quality control

``` text
01_quality_control/
├── raw_fastqc/
├── trimmed_fastqc/
├── fastp/
├── trimmed_reads/
└── multiqc/
```

These directories contain raw-read QC, processed-read QC, fastp reports, trimmed reads, and project-level MultiQC reports.

## Module 02 — Host depletion

``` text
02_host_depletion/
├── depleted_reads/
├── fastqc/
└── multiqc/
```

The `depleted_reads/` directory contains the host-depleted paired-end reads used by downstream metagenomic analyses.

## Module 03 — Taxonomic profiling

``` text
03_taxonomy/
├── kraken2/
└── bracken/
```

Kraken2 classification results and Bracken abundance estimates are stored here.

## Module 04 — Read-based resistome

``` text
04_read_based_resistome/
├── rgi/
└── sequencing_depth.csv
```

`sequencing_depth.csv` records the R1 read count, R2 read count, and the pipeline-defined `total_sequences` value used for downstream RPKM normalization.

## Module 05 — Metagenome assembly

``` text
05_assembly/
└── megahit/
    └── SAMPLE_ID/
        └── final.contigs.fa
```

The principal per-sample assembly is `final.contigs.fa`.

## Module 06 — Assembly quality

``` text
06_assembly_quality/
```

This directory contains MetaQUAST results when Module 06 is enabled. **Module 06 is disabled by default** because it can substantially increase runtime.

## Module 07 — Contig annotation

``` text
07_contig_annotation/
├── abricate_card/
└── plasmidfinder/
```

ABRicate-CARD results provide contig-level ARG detections, while ABRicate-PlasmidFinder results provide plasmid replicon detections.

## Module 08 — Downstream R analysis

``` text
08_downstream_analysis/
├── analysis_configuration.txt
├── objects/
├── taxonomy/
├── resistome/
└── plasmids/
```

### Taxonomy outputs

Typical outputs include:

``` text
taxonomy/
├── genus_relative_abundance.tsv
├── genus_sample_summary.tsv
├── genus_summary.tsv
├── phylum_relative_abundance.tsv
├── phylum_sample_summary.tsv
├── phylum_summary.tsv
├── species_relative_abundance.tsv
├── species_sample_summary.tsv
├── species_summary.tsv
├── taxonomy_analysis_summary.txt
├── top_genus.tsv
├── top_phylum.tsv
└── top_species.tsv
```

### Resistome outputs

Typical outputs include:

``` text
resistome/
├── abricate_arg_prevalence.csv
├── abricate_arg_richness.csv
├── abricate_card_high_confidence_hits.csv
├── resistome_sample_summary.csv
├── rgi_arg_drug_class_mapping.csv
├── rgi_arg_prevalence.csv
├── rgi_arg_richness.csv
├── rgi_arg_sample_table.csv
├── rgi_drug_class_abundance.csv
├── rgi_drug_class_arg_counts.csv
├── rgi_drug_class_global_summary.csv
├── rgi_drug_class_prevalence.csv
├── rgi_high_confidence_hits.csv
└── rgi_resistome_load.csv
```

### Plasmid outputs

Typical outputs include:

``` text
plasmids/
├── abricate_arg_contig_table.csv
├── abricate_card_parsed_hits.csv
├── arg_plasmid_contig_colocalization.csv
├── arg_plasmid_contig_colocalization_details.csv
├── arg_plasmid_sample_cooccurrence.csv
├── arg_profile.csv
├── plasmid_analysis_summary.txt
├── plasmidfinder_parsed_hits.csv
├── plasmid_overall_summary.csv
├── plasmid_prevalence.csv
├── plasmid_profile.csv
├── plasmid_richness.csv
├── plasmid_sample_summary.csv
└── plasmid_sample_table.csv
```

The exact set of files depends on the enabled calculation switches in `config/config.sh` and on whether relevant detections are present in the input data.

## Logs

Pipeline execution is logged to:

``` text
logs/pipeline.log
```

The log records the pipeline version, execution mode, start and finish times, resources, module status, output locations, and runtime.

## Configuration record

The downstream R workflow records its active configuration in:

``` text
08_downstream_analysis/analysis_configuration.txt
```

This file should be retained with the results when documenting an analysis.
