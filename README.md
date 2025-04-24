# ImmunoTools
## HLA Typing and Neoantigen Prediction Pipeline


This repository contains two integrated pipelines for analysing Whole Exome Sequencing (WES) data:

1. **HLA Typing Pipeline** – Performs automated HLA genotyping using OptiType.
2. **Neoantigen Pipeline** – Identifies somatic variants and predicts candidate neoantigens using pVACseq.

> **Run the HLA Typing Pipeline first** to obtain patient-specific HLA alleles. These are required inputs for the Neoantigen Pipeline.

---

## 📦 Requirements

- [Docker](https://www.docker.com/products/docker-desktop) installed and running
- Python 3.7+
- Snakemake (optional)

### Pull All Required Docker Images

Use the included shell script to download all necessary containers:

##  Running the Pipeline

### Step 1: Pull Docker images
```bash
./pull_dockers.sh
```

### Step 2: Run Snakemake
```bash
snakemake --use-singularity -j 4
```

> If you're using Docker locally instead of Singularity, you can still use `--use-singularity` and Snakemake will handle it under the hood. Make sure Docker is running and the container names in `pipeline_config.yaml` are correct.

This will execute:
- HLA typing using `hla_pipeline.py`
- Neoantigen prediction using `neoantigen_pipeline.py`, using output HLA alleles

## Main Snakemake Structure
```
ImmunoTools/
├── Snakefile
├── pipeline_config.yaml
├── data/
│   ├── hla/
│   │   ├── SAMPLE_001/
│   │   │   ├── SAMPLE_001_normal_1.fastq.gz
│   │   │   └── SAMPLE_001_normal_2.fastq.gz
│   │   ├── SAMPLE_002/
│   │   │   ├── SAMPLE_002_normal_1.fastq.gz
│   │   │   └── SAMPLE_002_normal_2.fastq.gz
│   └── neo/
│       ├── SAMPLE_001/
│       │   ├── SAMPLE_001_tumor_R1.fastq.gz
│       │   ├── SAMPLE_001_tumor_R2.fastq.gz
│       │   ├── SAMPLE_001_normal_R1.fastq.gz
│       │   └── SAMPLE_001_normal_R2.fastq.gz
│       └── SAMPLE_002/
│           ├── SAMPLE_002_tumor_R1.fastq.gz
│           ├── SAMPLE_002_tumor_R2.fastq.gz
│           ├── SAMPLE_002_normal_R1.fastq.gz
│           └── SAMPLE_002_normal_R2.fastq.gz
├── results/
│   ├── hla/
│   └── neoantigen/
├── pull_dockers.sh
```
This will pull the following versions:

| Tool          | Docker Image & Version                             |
|---------------|-----------------------------------------------------|
| BWA           | `biocontainers/bwa:v0.7.17_cv1`                     |
| Samtools      | `biocontainers/samtools:v1.3.1_cv4`                |
| FastQC        | `staphb/fastqc:0.11.9`                             |
| Trim Galore   | `quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0` |
| OptiType      | `fred2/optitype`                                   |
| Picard        | `broadinstitute/picard`                            |
| GATK          | `broadinstitute/gatk`                              |
| VEP           | `ensemblorg/ensembl-vep`                           |
| pVACtools     | `griffithlab/pvactools`                            |

---

You can also run the pipeline manually using raw commands.

## HLA Typing Pipeline

This pipeline performs automated HLA genotyping using paired-end WES data and Dockerised tools.

### Directory Structure

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

### OptiType Config Example

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

### Run the HLA Typing Pipeline

```bash
python3 hla_pipeline.py \
  --input_dir /path/to/samples \
  --genome_fasta /path/to/GRCh38.primary_assembly.genome.fa \
  --optitype_config /path/to/optitype_config.ini
```

### Output Per Sample

Each sample folder will contain:
- `FASTQC/` – Raw and trimmed quality reports
- `Trimmed/` – Trimmed FASTQ reads
- `aligned.sam`, `aligned.bam`, `aligned.sorted.bam` – Alignment files
- `HLA.bam`, `HLA_R1.fastq`, `HLA_R2.fastq` – HLA reads only
- `OptiType/` – Final HLA typing TSV with predicted Class I alleles

---

## Neoantigen Pipeline

This pipeline processes the WES data to detect somatic mutations and predict neoantigens.

### Required Inputs

- Tumor and normal paired-end FASTQ files
- Reference genome FASTA (GRCh38)
- Class I HLA alleles (from HLA Typing Pipeline)
- Local Ensembl VEP cache directory
- VEP plugins directory (must include `Wildtype.pm` and `Frameshift.pm`)

### Example One-liner

```bash
python3 neoantigen_pipeline.py \
  --tumor_r1 "/path/to/tumor_R1.fastq.gz" \
  --tumor_r2 "/path/to/tumor_R2.fastq.gz" \
  --normal_r1 "/path/to/normal_R1.fastq.gz" \
  --normal_r2 "/path/to/normal_R2.fastq.gz" \
  --genome_fasta "/path/to/GRCh38.primary_assembly.genome.fa" \
  --output_dir "/output/neoantigen" \
  --hla HLA-A*02:13 HLA-A*68:01 HLA-B*40:01 HLA-B*49:01 HLA-C*03:04 HLA-C*06:02 \
  --sampleid 48354 \
  --vep_cache "/path/to/homo_sapiens" \
  --plugin_dir "/path/to/vep_plugins" \
  --threads 4
```

### Output Structure

The Neoantigen Pipeline produces:
- `Trimmed/` – Trimmed FASTQ reads
- `*.bam` – Aligned, sorted, deduplicated BAM files
- `*.vcf.gz` – Raw and filtered somatic variant calls
- `annotated.vcf.gz` – VEP-annotated VCF (with Frameshift/Wildtype annotations)
- `pvacseq_output/` – Neoantigen predictions from pVACseq

### Notes
- Only **Class I HLA alleles** are supported.
- Default peptide lengths are **8, 9, 10, 11**.
- VEP must run successfully **with plugins** before using pVACseq.

---

## Full Example Workflow

```bash
# Stage 1: Run HLA Typing
python3 hla_pipeline.py \
  --input_dir ./input \
  --genome_fasta ./GRCh38.fa/GRCh38.primary_assembly.genome.fa \
  --optitype_config ./optitype_config.ini

# Stage 2: Use HLA alleles in Neoantigen Pipeline
python3 neoantigen_pipeline.py \
  --tumor_r1 ./input/48354-tumour_R1.fastq.gz \
  --tumor_r2 ./input/48354-tumour_R2.fastq.gz \
  --normal_r1 ./input/48354-normal_R1.fastq.gz \
  --normal_r2 ./input/48354-normal_R2.fastq.gz \
  --genome_fasta ./GRCh38.fa/GRCh38.primary_assembly.genome.fa \
  --output_dir ./output/neo_48354 \
  --hla HLA-A*02:13 HLA-A*68:01 HLA-B*40:01 HLA-B*49:01 HLA-C*03:04 HLA-C*06:02 \
  --sampleid 48354 \
  --vep_cache ./homo_sapiens \
  --plugin_dir ./vep_plugins \
  --threads 4
```
### VEP Cache & Plugin Setup
To run the neoantigen pipeline, you will need the Ensembl VEP cache and plugin directory:

- **VEP Cache** (`--vep_cache`):
  - Download from Ensembl using VEP installer:
    ```bash
    vep_install --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --PLUGINS all --DESTDIR ./vep_data
    ```
  - Set `--vep_cache ./vep_data/homo_sapiens`

- **Plugin Directory** (`--plugin_dir`):
  - Usually created during VEP install in the same destination directory:
    ```bash
    --plugin_dir ./vep_data/Plugins
    ```
  - If missing, you can clone manually:
    ```bash
    git clone https://github.com/Ensembl/VEP_plugins.git ./vep_plugins
    ```

---

## Author
Jayden Beckwith

