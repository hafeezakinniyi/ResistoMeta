# ResistoMeta

**ResistoMeta** is a modular and reproducible workflow for **metagenomic antimicrobial resistance (AMR) characterization** from paired-end sequencing reads. It combines read-level and assembly-level approaches to profile taxonomy, antimicrobial resistance genes (ARGs), and plasmid replicons, followed by downstream analysis in R. ResistoMeta is designed as a general-purpose and extensible workflow for researchers investigating antimicrobial resistance and microbial communities in metagenomic datasets across diverse One Health, environmental, animal, and human-associated settings.

## Workflow

``` text
config/config.sh
      │
      ▼
validate_installation.sh
      │
      ▼
run_pipeline.sh
      │
      ├── 01 Quality Control
      ├── 02 Host Depletion
      ├── 03 Taxonomic Profiling
      ├── 04 Read-based Resistome
      ├── 05 Metagenome Assembly
      ├── 06 Assembly Quality (optional)
      ├── 07 Contig Annotation
      │       ├── ABRicate-CARD
      │       └── ABRicate-PlasmidFinder
      │
      └── 08 Downstream R Analysis
              ├── Taxonomy
              ├── Resistome
              └── Plasmid / ARG localization
```

The detailed workflow, configuration, module dependencies, and output reference are available in [`docs/PIPELINE_DETAILS.md`](docs/PIPELINE_DETAILS.md) and [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## Main components

| Module | Function | Main tools |
|----|----|----|
| 01 | Read QC and preprocessing | FastQC, fastp, MultiQC |
| 02 | Host-read depletion | Bowtie2 |
| 03 | Taxonomic profiling | Kraken2, Bracken |
| 04 | Read-based ARG detection | RGI-CARD BWT workflow |
| 05 | Metagenome assembly | MEGAHIT |
| 06 | Assembly quality assessment | MetaQUAST *(optional)* |
| 07 | Contig AMR and plasmid annotation | ABRicate-CARD, ABRicate-PlasmidFinder |
| 08 | Downstream analysis | R |

## Requirements

ResistoMeta is intended for Linux systems and uses Conda environments for command-line and R dependencies.

Required reference resources are configured in:

``` text
config/config.sh
```

The current configuration includes paths for:

### Host reference

``` bash
HOST_REFERENCE_INDEX="/path/to/bowtie2_host_index"
```

This must point to the Bowtie2 index prefix, not to an individual `.bt2` file.

### Kraken2 database

``` bash
KRAKEN_DATABASE="/path/to/kraken2_database"
```

### ABRicate databases

``` bash
ABRICATE_CARD_DB="card"
ABRICATE_PLASMID_DB="plasmidfinder"
```

See [`docs/RGI_CARD_DATABASE.md`](docs/RGI_CARD_DATABASE.md) for RGI-CARD database setup.

## Installation

Clone the repository and enter the project directory:

``` bash
git clone https://github.com/hafeezakinniyi/ResistoMeta.git
cd ResistoMeta
```

Create the Conda environments defined in `envs/`:

``` bash
conda env create -f envs/main_environment.yml
conda env create -f envs/rgi_environment.yml
conda env create -f envs/r_environment.yml
```

Create the required project directories if they are not already present:

``` bash
mkdir -p \
    raw_data \
    databases \
    logs \
    results
```

Before running the pipeline, configure the paths and analysis options in:

``` text
config/config.sh
```

Then validate the installation:

``` bash
bash scripts/bash/validate_installation.sh
```

## Input data

The pipeline expects paired-end compressed FASTQ files in:

``` text
raw_data/
```

The master launcher discovers samples from files matching an R1 pattern such as:

``` text
SAMPLE_ID_R1_001.fastq.gz
```

The corresponding paired R2 file should be present:

``` text
SAMPLE_ID_R2_001.fastq.gz
```

The sequencing-depth table used by downstream resistome analysis is generated automatically by Module 04 from the actual host-depleted RGI input reads. For paired-end reads, `N` R1 reads and `N` R2 reads represent `N` sequencing fragments; the pipeline therefore records `total_sequences = N`, not `2N`.

## Running the pipeline

The main entry point is:

``` bash
bash scripts/bash/run_pipeline.sh
```

The pipeline validates the installation before execution and records execution information in:

``` text
logs/pipeline.log
```

## Results

Results are organized by processing stage:

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

The detailed output files and their purposes are documented in [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## Design principles

- **Modular execution:** major processing stages are implemented as independent modules.
- **Centralized configuration:** `config/config.sh` is the single configuration source.
- **Environment isolation:** bioinformatics, RGI, and downstream R analysis use separate Conda environments.
- **Read- and assembly-level evidence:** RGI and ABRicate results are retained as complementary AMR evidence.
- **Primary detection versus downstream analysis:** Bash modules perform primary processing; R modules perform downstream summaries and analyses.
- **Validation and resumability:** expected inputs and outputs are checked, and completed processing can be reused where supported by the workflow.

## Documentation

- [`docs/PIPELINE_DETAILS.md`](docs/PIPELINE_DETAILS.md) — detailed workflow, configuration, dependencies, and execution behavior.
- [`docs/OUTPUTS.md`](docs/OUTPUTS.md) — output directories and important result files.
- [`docs/RGI_CARD_DATABASE.md`](docs/RGI_CARD_DATABASE.md) — RGI-CARD database installation and configuration.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

If ResistoMeta contributes to your research, please cite the version used in your analysis.

## License

ResistoMeta is released under the **MIT License**. See [`LICENSE`](LICENSE).

## Author

**Hafeez T. Akinniyi**\
DVM \| MSc, Molecular Biology & Genomics\
GitHub: <https://github.com/hafeezakinniyi>

For questions, suggestions, or collaborations, please open an issue or submit a pull request.
