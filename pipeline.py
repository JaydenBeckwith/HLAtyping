import os
import glob
import argparse
import subprocess
from pathlib import Path

def run_command(cmd, description=""):
    print(f"\n[INFO] {description}")
    try:
        subprocess.run(cmd, check=True)
        print("[SUCCESS] Completed.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to {description}")
        print(e)

def bwa_index_if_missing(genome_fasta):
    bwa_index_files = [f"{genome_fasta}.{ext}" for ext in ["bwt", "pac", "ann", "amb", "sa"]]
    if not all(Path(f).exists() for f in bwa_index_files):
        print("\n[INFO] BWA index files not found. Indexing reference genome...")
        run_command([
            "docker", "run", "--rm",
            "-v", f"{Path(genome_fasta).parent}:/ref",
            "biocontainers/bwa:v0.7.17_cv1",
            "bwa", "index", f"/ref/{Path(genome_fasta).name}"
        ], f"indexing {Path(genome_fasta).name} for BWA")
    else:
        print("[INFO] BWA index already exists. Skipping indexing.")

def main(input_dir, genome_fasta, config_path):
    extract_dir = Path(input_dir).resolve()
    genome_fasta = Path(genome_fasta).resolve()
    config_path = Path(config_path).resolve()

    bwa_index_if_missing(genome_fasta)

    for sample_dir in sorted(extract_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample_name = sample_dir.name
        print(f"\n[INFO] Processing sample: {sample_name}")

        fastqc_dir = sample_dir / "FASTQC"
        trimmed_dir = sample_dir / "Trimmed"
        fastqc_dir.mkdir(exist_ok=True)
        trimmed_dir.mkdir(exist_ok=True)

        fastq_r1_list = list(sample_dir.glob("*_1*.fastq.gz"))
        fastq_r2_list = list(sample_dir.glob("*_2*.fastq.gz"))
        if not fastq_r1_list or not fastq_r2_list:
            print(f"[WARNING] FASTQ files not found for {sample_name}, skipping...")
            continue

        fastq_r1 = fastq_r1_list[0].name
        fastq_r2 = fastq_r2_list[0].name

        # Step 1: FastQC (raw)
        run_command([
            "docker", "run", "--rm", "-v", f"{extract_dir}:/mnt/input", "-v", f"{extract_dir}:/mnt/output",
            "staphb/fastqc:0.11.9", "fastqc", "-o", f"/mnt/output/{sample_name}/FASTQC",
            f"/mnt/input/{sample_name}/{fastq_r1}", f"/mnt/input/{sample_name}/{fastq_r2}"
        ], f"FastQC (raw) for {sample_name}")

        # Step 2: Trim Galore
        run_command([
            "docker", "run", "--rm", "-v", f"{extract_dir}:/mnt/input", "-v", f"{extract_dir}:/mnt/output",
            "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0", "trim_galore", "--paired",
            "-o", f"/mnt/output/{sample_name}/Trimmed",
            f"/mnt/input/{sample_name}/{fastq_r1}", f"/mnt/input/{sample_name}/{fastq_r2}"
        ], f"Trim Galore for {sample_name}")

        # Step 3: FastQC (trimmed)
        trimmed_r1 = os.path.basename(list(trimmed_dir.glob("*_val_1.fq.gz"))[0])
        trimmed_r2 = os.path.basename(list(trimmed_dir.glob("*_val_2.fq.gz"))[0])
        run_command([
            "docker", "run", "--rm", "-v", f"{extract_dir}:/mnt/input", "-v", f"{extract_dir}:/mnt/output",
            "staphb/fastqc:0.11.9", "fastqc", "-o", f"/mnt/output/{sample_name}/FASTQC",
            f"/mnt/output/{sample_name}/Trimmed/{trimmed_r1}", f"/mnt/output/{sample_name}/Trimmed/{trimmed_r2}"
        ], f"FastQC (trimmed) for {sample_name}")

        # Step 4: Align and downstream processing
        container_path = f"/data/{sample_name}"
        trimmed_r1_path = f"{container_path}/Trimmed/{trimmed_r1}"
        trimmed_r2_path = f"{container_path}/Trimmed/{trimmed_r2}"

        # BWA mem
        run_command([
            "docker", "run", "--rm", "-v", f"{extract_dir}:/data", "-v", f"{genome_fasta.parent}:/ref",
            "--entrypoint", "bash", "biocontainers/bwa:v0.7.17_cv1", "-c",
            f"bwa mem -t 8 /ref/{genome_fasta.name} {trimmed_r1_path} {trimmed_r2_path} > {container_path}/aligned.sam"
        ], f"Aligning {sample_name} with BWA")

        # SAM to BAM, sort, index
        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "view", "-bS",
                     f"{container_path}/aligned.sam", "-o", f"{container_path}/aligned.bam"],
                    f"SAM to BAM for {sample_name}")

        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "sort",
                     f"{container_path}/aligned.bam", "-o", f"{container_path}/aligned.sorted.bam"],
                    f"Sorting BAM for {sample_name}")

        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "index",
                     f"{container_path}/aligned.sorted.bam"], f"Indexing BAM for {sample_name}")

        # Subset BAM to HLA
        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "view", "-b",
                     f"{container_path}/aligned.sorted.bam", "chr6:28477797-33448354",
                     "-o", f"{container_path}/HLA.bam"], f"Extract HLA region for {sample_name}")

        # Convert to FASTQ
        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "collate", "-O", "-u",
                     f"{container_path}/HLA.bam"], f"Collate BAM for {sample_name}")

        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "biocontainers/samtools:v1.3.1_cv4", "samtools", "fastq",
                     "-1", f"{container_path}/HLA_R1.fastq", "-2", f"{container_path}/HLA_R2.fastq",
                     "-0", "/dev/null", "-s", "/dev/null", "-n",
                     f"{container_path}/HLA.bam"], f"Convert BAM to FASTQ for {sample_name}")

        # Run OptiType
        run_command(["docker", "run", "--rm", "-v", f"{extract_dir}:/data",
                     "-v", f"{config_path}:/config.ini", "fred2/optitype",
                     "-i", f"{container_path}/HLA_R1.fastq", f"{container_path}/HLA_R2.fastq",
                     "--dna", "--verbose", "-o", f"{container_path}/OptiType", "--config", "/config.ini"],
                    f"Run OptiType for {sample_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HLA Typing Pipeline using Docker")
    parser.add_argument("--input_dir", required=True, help="Path to folder containing sample directories")
    parser.add_argument("--genome_fasta", required=True, help="Path to GRCh38.primary_assembly.genome.fa")
    parser.add_argument("--optitype_config", required=True, help="Path to optitype_config.ini")
    args = parser.parse_args()

    main(args.input_dir, args.genome_fasta, args.optitype_config)
