# daff-viral-assembly

A [Nextflow](https://www.nextflow.io/) pipeline for viral genome assembly and identification from Illumina short-read data. Developed through a collaboration between [QUT eResearch](https://www.qut.edu.au/research/eresearch) and [QCIF](https://www.qcif.edu.au/).

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Input](#input)
- [Pipeline Steps](#pipeline-steps)
- [Parameters](#parameters)
- [Output](#output)
- [Profiles](#profiles)
- [Tool Versions](#tool-versions)

## Overview

`daff-viral-assembly` detects and characterises known and novel viruses from paired-end Illumina RNA-seq data. It performs:

1. Read quality control and trimming
2. rRNA and PhiX decontamination
3. Metagenomic read classification (Kraken2 and Kaiju)
4. De novo viral genome assembly (SPAdes)
5. BLAST-based contig annotation
6. Reference mapping and consensus generation
7. ORF prediction and protein domain annotation
8. Novel virus candidate identification
9. Per-sample HTML report generation

## Requirements

- [Nextflow](https://nextflow.io/) ≥ 21.05.0
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) (recommended for reproducibility)
- The following databases (paths provided as parameters):
  - [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) nucleotide database (e.g. `core_nt`)
  - [TaxonKit](https://bioinf.shenwei.me/taxonkit/) taxonomy database
  - [Kraken2](https://ccb.jhu.edu/software/kraken2/) database (e.g. `core_nt`)
  - [Kaiju](https://github.com/bioinformatics-centre/kaiju) database (e.g. `kaiju_db_nr_euk`)
  - [GeNomad](https://github.com/apcamargo/genomad) database
  - [HMMER](http://hmmer.org/) Pfam database (`Pfam-A.hmm`)
  - rRNA reference sequences (for rRNA filtering)

## Quick Start

1. Install [Nextflow](https://nextflow.io/docs/latest/install.html) (≥ 21.05.0).

2. Install [Singularity](https://sylabs.io/guides/latest/user-guide/quick_start.html) or [Docker](https://docs.docker.com/get-docker/).

3. Create a samplesheet (see [Input](#input) for format).

4. Create a params file (copy and edit `params/params_example.yml`):

   ```yaml
   input: /path/to/samplesheet.csv
   blastn_db: /path/to/ncbi/core_nt
   taxdump: ~/.taxonkit
   blast_threads: 2
   genomad_db: /path/to/genomad/genomad_db
   hmmer_db: /path/to/hmms/Pfam-A.hmm
   kraken2_db: /path/to/kraken_databases/core_nt
   kaiju_dbname: /path/to/kaiju_databases/kaiju_db_nr_euk.fmi
   kaiju_names: /path/to/kaiju_databases/names.dmp
   kaiju_nodes: /path/to/kaiju_databases/nodes.dmp
   rrna_ref: /path/to/rrna_reference.fasta
   outdir: results
   ```

5. Run the pipeline:

   ```bash
   nextflow run main.nf \
     -profile singularity \
     -params-file params/user_params.yml \
     --analyst_name "Your Name" \
     --facility "Your Institution"
   ```

   To resume a failed run, add `-resume`:

   ```bash
   nextflow run main.nf \
     -profile singularity \
     -params-file params/user_params.yml \
     --analyst_name "Your Name" \
     --facility "Your Institution" \
     -resume
   ```

## Input

### Samplesheet

Provide a CSV samplesheet via `--input` (default: `index.csv`). The file must include the following columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier (no spaces) |
| `fastq_1` | Yes | Absolute path to R1 FASTQ file (`.fastq.gz` or `.fq.gz`) |
| `fastq_2` | No | Absolute path to R2 FASTQ file (`.fastq.gz` or `.fq.gz`). Leave blank for single-end data |
| `sample_information` | No | Free-text sample information |
| `sample_type` | No | Sample type descriptor |
| `sample_receipt_date` | No | Date sample was received |
| `storage_location` | No | Physical storage location of the sample |

**Example samplesheet (`index.csv`):**

```csv
sample,fastq_1,fastq_2,sample_information,sample_type,sample_receipt_date,storage_location
Sample_A,/data/Sample_A_R1.fastq.gz,/data/Sample_A_R2.fastq.gz,Tomato leaf,Plant,2024-01-15,Freezer_1
Sample_B,/data/Sample_B_R1.fastq.gz,/data/Sample_B_R2.fastq.gz,Rose stem,Plant,2024-01-16,Freezer_1
```

Multiple FASTQ files from the same sample (e.g. different sequencing lanes) can be listed as separate rows with the same `sample` identifier — they will be concatenated before processing.

## Pipeline Steps

```
Input FASTQ files
      │
      ▼
 CAT_FASTQ ──────── Concatenate files from the same sample
      │
      ▼
 (optional) SEQTK_SAMPLE ── Subsample to 60 M reads (default)
      │
      ├──► FASTQC_RAW ──── Raw read quality control
      │
      ▼
 FASTP ──────────── Adapter trimming and quality filtering
      │
      ├──► FASTQC_TRIM ─── Post-trimming quality control
      │
      ▼
 BBMAP_BBDUK ────── rRNA read removal
      │
      ▼
 BBMAP_BBSPLIT ──── PhiX decontamination
      │
      ├──► KRAKEN2_KRAKEN2 ──► BRACKEN ──► KRAKEN2_TO_KRONA ──► (Krona chart)
      │         │
      │         └──► RETRIEVE_VIRAL_READS_KRAKEN2 ── Extract viral + unclassified reads
      │                         │
      ├──► KAIJU_KAIJU ──────►  │          ── Protein-level classification + Krona chart
      │         │               │
      │         └──► SUMMARISE_READ_CLASSIFICATION
      │                         │
      ▼                         ▼
 SPADES ─────────── De novo assembly (rnaSPAdes)
      │
      ▼
 SEQTK ──────────── Filter contigs < 150 bp
      │
      ▼
 BLASTN ─────────── Megablast search (split input for parallel execution)
      │
      ▼
 EXTRACT_VIRAL_BLAST_HITS ── Identify viral contigs
      │
      ├──► EXTRACT_CONTIGS ─── Separate viral / other contigs
      │         │
      │         ▼
      │    MAPPING_BACK_TO_CONTIGS ──► PILEUP ──► IDENTIFY_ERRORS ──► TRIM_ENDS
      │              │                                                      │
      │              └──────────────────────────── REALIGN ◄────────────────┘
      │                                               │
      │                          SAMTOOLS_CONTIGS ◄───┘
      │                          MOSDEPTH_CONTIGS
      │                          PYFAIDX_CONTIGS
      │
      ├──► BLASTN_ROUND2 ──── Second BLAST on trimmed contigs
      │         │
      │         ▼
      │    EXTRACT_VIRAL_BLAST_HITS_ROUND2 ──► FASTA2TABLE
      │                                              │
      │              ┌───────────────────────────────┘
      │              ▼
      │    EXTRACT_REF_FASTA ──► CLUSTER ──► MAPPING_BACK_TO_REF
      │                                            │
      │                                       SAMTOOLS2 ──► BCFTOOLS ──► BEDTOOLS
      │                                       PYFAIDX / MOSDEPTH
      │                                       REF_COVSTATS ──► FASTA2TABLE2
      │                                       CONTIG_COVSTATS
      │
      ├──► ORFIPY ──────────── ORF prediction from contigs
      │         │
      │         ▼
      │    HMMSCAN ─────────── Protein domain annotation (Pfam)
      │
      ├──► GENOMAD ─────────── Virus / provirus classification
      │
      ▼
 SUMMARISE_RESULTS ── Combine all evidence into per-sample summary table
      │
 NOVELS ──────────── Identify novel virus candidates
      │
 QCREPORT ────────── Run-level QC report
      │
 HTML_REPORT ──────── Per-sample interactive HTML report
```

### Step Descriptions

| Step | Tool | Description |
|------|------|-------------|
| Read concatenation | `cat` | Merges multiple FASTQ files per sample |
| Subsampling | seqtk | Randomly subsample reads to a target depth (default 60 M reads) |
| Raw QC | FastQC | Quality metrics for raw reads |
| Adapter trimming | fastp | Trims adapters, low-quality bases, and short reads |
| Trimmed QC | FastQC | Quality metrics after trimming |
| rRNA removal | BBDuk | Removes ribosomal RNA reads using a reference |
| PhiX removal | BBSplit | Removes PhiX174 spike-in reads |
| Taxonomic classification | Kraken2 + Bracken | k-mer based classification; Bracken re-estimates abundance |
| Taxonomic classification | Kaiju | Protein-level classification against nr/euk database |
| Visualisation | Krona | Interactive pie-chart visualisation of classifications |
| Viral read retrieval | KrakenTools | Extracts unclassified + viral-classified reads |
| De novo assembly | SPAdes (rnaSPAdes mode) | Assembles viral genome contigs |
| Contig filtering | seqtk | Removes contigs shorter than 150 bp |
| BLAST annotation | BLASTN | Searches assembled contigs against the NCBI nucleotide database |
| Hit extraction | custom script | Filters to top viral BLAST hits |
| Reference retrieval | Entrez Direct | Downloads matched reference sequences from NCBI |
| Reference clustering | CD-HIT | Clusters highly similar references |
| Reference mapping | BWA | Aligns cleaned reads back to reference sequences |
| Consensus calling | bcftools | Calls variants and applies them to generate consensus sequences |
| Coverage masking | bedtools | Masks low-coverage regions in the consensus |
| Depth analysis | mosdepth | Calculates per-base sequencing depth |
| ORF prediction | orfipy | Predicts open reading frames from assembled contigs |
| Protein domain annotation | HMMER (hmmscan) | Searches ORFs against Pfam protein domain database |
| Virus/provirus prediction | GeNomad | Classifies contigs as viral or proviral |
| Results summary | custom script | Integrates evidence from all tools into a summary table |
| Novel virus detection | custom script | Identifies potential novel virus candidates |
| QC report | custom script | Generates run-level quality control report |
| HTML report | custom script | Generates interactive per-sample HTML report |

## Parameters

Parameters can be provided via a YAML params file (`-params-file`) or on the command line (`--param_name value`).

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--blastn_db` | Path to BLAST nucleotide database (e.g. NCBI `core_nt`) |
| `--taxdump` | Path to TaxonKit taxonomy database directory |
| `--kraken2_db` | Path to Kraken2 database directory |
| `--kaiju_dbname` | Path to Kaiju `.fmi` index file |
| `--kaiju_names` | Path to Kaiju `names.dmp` file |
| `--kaiju_nodes` | Path to Kaiju `nodes.dmp` file |
| `--genomad_db` | Path to GeNomad database directory |
| `--hmmer_db` | Path to Pfam HMM database file (`Pfam-A.hmm`) |
| `--rrna_ref` | Path to rRNA reference FASTA file (for BBDuk filtering) |
| `--analyst_name` | Name of the analyst (included in reports) |
| `--facility` | Name of the sequencing/analysis facility (included in reports) |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | `index.csv` | Path to samplesheet CSV file |
| `--outdir` | `results` | Output directory |
| `--blast_threads` | `2` | Number of CPU threads for BLASTN |
| `--subsample_enabled` | `true` | Enable/disable read subsampling before processing |
| `--subsample_size` | `60000000` | Target number of reads after subsampling (60 M) |
| `--kraken2_save_classified_reads` | `false` | Save Kraken2-classified reads to output |
| `--kraken2_save_unclassified_reads` | `true` | Save Kraken2-unclassified reads to output |
| `--kraken2_save_readclassifications` | `true` | Save per-read Kraken2 classification output |
| `--save_trimmed_fail` | `false` | Save reads that failed fastp quality trimming |
| `--save_merged` | `false` | Save fastp-merged reads |
| `--publish_dir_mode` | `copy` | Nextflow output publish mode (`copy`, `symlink`, `link`, etc.) |
| `--help` | `false` | Print help message and exit |

## Output

All results are written to `--outdir` (default: `results`). Per-sample results are nested under the sample identifier.

```
results/
├── 01_pipeline_logs/
│   ├── execution_timeline_<timestamp>.html   # Nextflow timeline report
│   ├── execution_report_<timestamp>.html     # Nextflow execution report
│   ├── execution_trace_<timestamp>.txt       # Process trace file
│   ├── pipeline_dag_<timestamp>.html         # Pipeline DAG visualisation
│   └── <timestamp>_nextflow_start_timestamp.txt
│
├── 02_qc_report/
│   ├── run_qc_report_<timestamp>.html        # Run-level QC summary (HTML)
│   └── run_qc_report_<timestamp>.txt         # Run-level QC summary (text)
│
└── <sample_id>/
    ├── 05_read_classification/
    │   ├── <sample_id>_kraken2_report.txt          # Kraken2 report
    │   ├── <sample_id>_kraken2_to_krona.html       # Kraken2 Krona chart
    │   ├── <sample_id>_kaiju_summary.txt           # Kaiju classification summary
    │   ├── <sample_id>_kraken_summary.txt          # Kraken2 classification summary
    │   └── <sample_id>_bracken_report.txt          # Bracken abundance report
    │
    ├── 06_assembly/
    │   └── <sample_id>_spades_scaffolds.fasta      # Filtered assembled contigs (≥ 150 bp)
    │
    ├── 07_annotation/
    │   ├── <sample_id>_megablast_top_viral_hits.txt          # Top viral BLAST hits
    │   ├── <sample_id>_orfs.fasta                            # Predicted ORFs (orfipy)
    │   ├── <sample_id>_hmmscan_per_target_output.txt         # HMMER per-target results
    │   ├── <sample_id>_hmmscan_per_domain_output.txt         # HMMER per-domain results
    │   └── genomad/                                          # GeNomad virus predictions
    │
    ├── 08_mapping_to_contigs/
    │   ├── <sample_id>_contig_aln.sorted.bam                 # BAM alignment to contigs
    │   ├── <sample_id>_contig_aln.sorted.bam.bai             # BAM index
    │   ├── <sample_id>_contig_coverage.txt                   # Per-base contig coverage
    │   ├── <sample_id>_contig_mapq.txt                       # Mapping quality statistics
    │   ├── <sample_id>_contigs.trimmed.fa                    # End-trimmed contig sequences
    │   └── <sample_id>_trim_coords.tsv                       # Trimming coordinates
    │
    ├── 09_mapping_to_ref/
    │   ├── <sample_id>_ref_aln.sorted.bam                    # BAM alignment to references
    │   ├── <sample_id>_ref_aln.sorted.bam.bai                # BAM index
    │   ├── <sample_id>_ref_coverage.txt                      # Per-base reference coverage
    │   ├── <sample_id>_vcf_applied.fasta                     # Consensus with variants applied
    │   ├── <sample_id>_bcftools_masked_consensus.fasta       # Coverage-masked consensus
    │   └── <sample_id>_reference_with_cov_stats_final.txt    # Final reference stats table
    │
    ├── 10_results_summary/
    │   ├── <sample_id>_summary_viral_results.tsv             # Known virus summary table
    │   ├── <sample_id>_novel_virus_candidates.tsv            # Novel virus candidates table
    │   └── <sample_id>_hmm_domain_summary_counts.tsv        # HMM domain summary
    │
    └── 11_report/
        └── <sample_id>_report.html                           # Interactive per-sample report
```

## Profiles

The pipeline includes several execution profiles selectable with `-profile`:

| Profile | Description |
|---------|-------------|
| `docker` | Run processes inside Docker containers |
| `singularity` | Run processes inside Singularity containers (recommended on HPC) |
| `mtdt_test` | Minimal test profile for local testing |
| `peq_test` | PEQ test profile |
| `internal_test` | Internal test profile |

Multiple profiles can be combined with commas, e.g. `-profile singularity,mtdt_test`.

## Tool Versions

| Tool | Version | Purpose |
|------|---------|---------|
| [BBMap](https://sourceforge.net/projects/bbmap/) | 39.37 | rRNA and PhiX decontamination |
| [bcftools](https://samtools.github.io/bcftools/) | 1.19 | Variant calling and consensus generation |
| [bedtools](https://bedtools.readthedocs.io/) | 2.27.1 | Coverage masking |
| [BLASTN](https://blast.ncbi.nlm.nih.gov/) | 2.16.0 | Nucleotide sequence alignment |
| [CD-HIT](https://sites.google.com/view/cd-hit) | 4.8.1 | Reference sequence clustering |
| [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) | 22.4 | NCBI sequence retrieval |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.12.1 | Read quality control |
| [fastp](https://github.com/OpenGene/fastp) | 1.0.1 | Adapter trimming and quality filtering |
| [GeNomad](https://github.com/apcamargo/genomad) | 1.11.2 | Virus / provirus classification |
| [HMMER](http://hmmer.org/) | 3.4.0 | Protein domain annotation |
| [Kaiju](https://github.com/bioinformatics-centre/kaiju) | 1.8.2 | Protein-level taxonomic classification |
| [Kraken2](https://ccb.jhu.edu/software/kraken2/) | 2.1.3 | k-mer based taxonomic classification |
| [KrakenTools](https://github.com/jenniferlu717/KrakenTools) | 1.2.1 | Viral read extraction from Kraken2 output |
| [Krona](https://github.com/marbl/Krona) | 2.8.1 | Interactive taxonomic visualisation |
| [mosdepth](https://github.com/brentp/mosdepth) | 0.3.3 | Sequencing depth analysis |
| [orfipy](https://github.com/urmi-21/orfipy) | 0.0.4 | ORF prediction |
| [pyfaidx](https://github.com/mdshw5/pyfaidx) | 0.8.1.3 | FASTA indexing |
| [Python](https://www.python.org/) | 3.12.9 | Scripting |
| [SAMtools](https://www.htslib.org/) | 1.22.1 | BAM processing |
| [seqtk](https://github.com/lh3/seqtk) | 1.3 | Sequence subsampling and filtering |
| [SPAdes](https://github.com/ablab/spades) | 4.2.0 | De novo genome assembly |
