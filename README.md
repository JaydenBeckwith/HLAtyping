# 🧬 HLA Typing Pipeline

This pipeline performs automated HLA genotyping from Whole Exome Sequencing (WES) data using Dockerised bioinformatics tools.

It supports multiple samples in one run, using a single command-line interface, and is platform-agnostic as long as Docker is installed.

Currently only supports paired-end reads.

---

## 📦 Requirements

- [Docker](https://www.docker.com/products/docker-desktop) installed and running
- Python 3.7+

### Pull Docker Images Automatically

You can use the included `pull_dockers.sh` to pull all required Docker containers:

```bash
./pull_dockers.sh
```

This will ensure all Docker images are available locally before running the pipeline.

---

## Docker Images Used

The following containers will be automatically pulled when the shell script runs:

| Tool          | Docker Image                                      |
|---------------|----------------------------------------------------|
| BWA           | `biocontainers/bwa:v0.7.17_cv1`                    |
| Samtools      | `biocontainers/samtools:v1.3.1_cv4`               |
| FastQC        | `staphb/fastqc:0.11.9`                            |
| Trim Galore   | `quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0` |
| OptiType      | `fred2/optitype`                                  |

---

## Directory Structure

Your input folder should be structured as follows:

```
input_dir/
├── SAMPLE_001/
│   ├── SAMPLE_001_1.fastq.gz
│   ├── SAMPLE_001_2.fastq.gz
├── SAMPLE_002/
│   ├── SAMPLE_002_1.fastq.gz
│   ├── SAMPLE_002_2.fastq.gz
...
```

Each sample folder must contain paired-end FASTQ files.

---

## OptiType Config

A sample `optitype_config.ini` file:

```ini
[mapping]
razers3=/usr/local/bin/razers3
threads=8

[ilp]
solver=cbc
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
```

Place this in your repo or project folder and pass its path using `--optitype_config`.

---

## How to Run

```bash
python3 pipeline.py \
  --input_dir /path/to/samples \
  --genome_fasta /path/to/GRCh38.primary_assembly.genome.fa \
  --optitype_config /path/to/optitype_config.ini
```

### Parameters:
| Argument | Description |
|----------|-------------|
| `--input_dir` | Directory containing sample subfolders |
| `--genome_fasta` | Path to GRCh38 `.fa` file (e.g., `GRCh38.primary_assembly.genome.fa`) |
| `--optitype_config` | Path to OptiType config file |

---

## Output Per Sample

Each sample folder will be populated with:

- `FASTQC/` – Raw and trimmed quality reports
- `Trimmed/` – Trimmed FASTQ reads
- `aligned.sam`, `aligned.bam`, `aligned.sorted.bam` – Alignment files
- `HLA.bam`, `HLA_R1.fastq`, `HLA_R2.fastq` – HLA reads only
- `OptiType/` – Final HLA typing output
  - For the OptiType output you should expect a tsv with only the optimal solution by default with a coverage distribution of the predicted alleles.   

---

## Author
Jayden Beckwith

