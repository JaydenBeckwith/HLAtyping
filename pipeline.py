import os
import glob
import subprocess

# Reference paths
EXTRACT_DIR = "C:/Users/meltest/Desktop/Opacin_DNA_HLA"
BASE_DIR = EXTRACT_DIR  # For FastQC/Trim mounting
GENOME_DIR = "C:/Users/meltest/Desktop/GRCh38.fa"
HG38_FASTA = os.path.join(GENOME_DIR, "GRCh38.primary_assembly.genome.fa")
HG38_BASENAME = os.path.basename(HG38_FASTA)
CONFIG_PATH = os.path.abspath("optitype_config.ini")

def run_command(cmd, description=""):
    print(f"\n[INFO] {description}")
    try:
        subprocess.run(cmd, check=True)
        print("[SUCCESS] Completed.")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to {description}")
        print(e)

# Step 0: BWA index if missing
bwa_index_files = [f"{HG38_FASTA}.{ext}" for ext in ["bwt", "pac", "ann", "amb", "sa"]]
if not all(os.path.exists(f) for f in bwa_index_files):
    print("\n[INFO] BWA index files not found. Indexing reference genome...")
    run_command([
        "docker", "run", "--rm",
        "-v", f"{GENOME_DIR}:/ref",
        "biocontainers/bwa:v0.7.17_cv1",
        "bwa", "index", f"/ref/{HG38_BASENAME}"
    ], "indexing GRCh38.primary_assembly.genome.fa for BWA")
else:
    print("[INFO] BWA index already exists. Skipping indexing.")

# Step 1: Loop through samples
for sample_name in os.listdir(EXTRACT_DIR):
    sample_dir = os.path.join(EXTRACT_DIR, sample_name)
    if not os.path.isdir(sample_dir):
        continue
    
    os.makedirs(os.path.join(sample_dir, "FASTQC"), exist_ok=True)
    os.makedirs(os.path.join(sample_dir, "Trimmed"), exist_ok=True)

    fastq_r1 = glob.glob(os.path.join(sample_dir, "*_1*.fastq.gz"))
    fastq_r2 = glob.glob(os.path.join(sample_dir, "*_2*.fastq.gz"))

    if not fastq_r1 or not fastq_r2:
        print(f"[WARNING] FASTQ files not found for {sample_name}, skipping...")
        continue

    fastq_r1 = os.path.basename(fastq_r1[0])
    fastq_r2 = os.path.basename(fastq_r2[0])
    container_sample_path = f"/mnt/input/{sample_name}"

    # Step 1.1: FastQC on raw FASTQ files
    run_command([
        "docker", "run", "--rm", "-it",
        "-v", f"{BASE_DIR}:/mnt/input",
        "-v", f"{BASE_DIR}:/mnt/output",
        "staphb/fastqc:0.11.9",
        "fastqc", "-o", f"/mnt/output/{sample_name}/FASTQC",
        f"/mnt/input/{sample_name}/{fastq_r1}",
        f"/mnt/input/{sample_name}/{fastq_r2}"
    ], f"Running FastQC for {sample_name} (raw files)")

    # Step 1.2: Run Trim Galore
    run_command([
        "docker", "run", "--rm", "-it",
        "-v", f"{BASE_DIR}:/mnt/input",
        "-v", f"{BASE_DIR}:/mnt/output",
        "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0",
        "trim_galore", "--paired", "-o", f"/mnt/output/{sample_name}/Trimmed",
        f"/mnt/input/{sample_name}/{fastq_r1}",
        f"/mnt/input/{sample_name}/{fastq_r2}"
    ], f"Running Trim Galore for {sample_name}")

     # Step 1.3: Dynamically find validated trimmed FASTQ files
    trimmed_dir = os.path.join(sample_dir, "Trimmed")
    trimmed_r1 = os.path.basename(glob.glob(os.path.join(trimmed_dir, "*_val_1.fq.gz"))[0])
    trimmed_r2 = os.path.basename(glob.glob(os.path.join(trimmed_dir, "*_val_2.fq.gz"))[0])

    # Run FastQC on trimmed FASTQ files
    run_command([
        "docker", "run", "--rm", "-it",
        "-v", f"{BASE_DIR}:/mnt/input",
        "-v", f"{BASE_DIR}:/mnt/output",
        "staphb/fastqc:0.11.9",
        "fastqc", "-o", f"/mnt/output/{sample_name}/FASTQC",
        f"/mnt/output/{sample_name}/Trimmed/{trimmed_r1}",
        f"/mnt/output/{sample_name}/Trimmed/{trimmed_r2}"
    ], f"Running FastQC for {sample_name} (trimmed files)")

    # Use validated trimmed FASTQ files for alignment
    container_sample_path = f"/data/{sample_name}"
    trimmed_r1_path = f"{container_sample_path}/Trimmed/{trimmed_r1}"
    trimmed_r2_path = f"{container_sample_path}/Trimmed/{trimmed_r2}"

    # Step 2: Align to GRCh38 using BWA
    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "-v", f"{GENOME_DIR}:/ref",
        "--entrypoint", "bash",
        "biocontainers/bwa:v0.7.17_cv1",
        "-c", f"bwa mem -t 8 /ref/GRCh38.primary_assembly.genome.fa "
              f"{trimmed_r1_path} {trimmed_r2_path} "
              f"> {container_sample_path}/aligned.sam"
    ], f"aligning {sample_name} with BWA and saving to aligned.sam")

    # Step 3: Convert to BAM, sort, and index
    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "view", "-bS",
        f"{container_sample_path}/aligned.sam",
        "-o", f"{container_sample_path}/aligned.bam"
    ], f"converting SAM to BAM for {sample_name}")

    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "sort",
        f"{container_sample_path}/aligned.bam",
        "-o", f"{container_sample_path}/aligned.sorted.bam"
    ], f"sorting BAM for {sample_name}")

    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "index",
        f"{container_sample_path}/aligned.sorted.bam"
    ], f"indexing BAM for {sample_name}")

    # Step 4: Subset BAM to HLA region
    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "view", "-b",
        f"{container_sample_path}/aligned.sorted.bam",
        "chr6:28477797-33448354",
        "-o", f"{container_sample_path}/HLA.bam"
    ], f"extracting HLA region for {sample_name}")

    # Step 5: Convert HLA BAM to paired FASTQ
    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "collate", "-O", "-u",
        f"{container_sample_path}/HLA.bam"
    ], f"collating HLA BAM for {sample_name}")

    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "biocontainers/samtools:v1.3.1_cv4",
        "samtools", "fastq",
        "-1", f"{container_sample_path}/HLA_R1.fastq",
        "-2", f"{container_sample_path}/HLA_R2.fastq",
        "-0", "/dev/null", "-s", "/dev/null", "-n",
        f"{container_sample_path}/HLA.bam"
    ], f"converting HLA BAM to FASTQ for {sample_name}")

    # Step 6: Run OptiType
    run_command([
        "docker", "run", "--rm",
        "-v", f"{EXTRACT_DIR}:/data",
        "-v", f"{CONFIG_PATH}:/config.ini",
        "fred2/optitype",
        "-i", f"{container_sample_path}/HLA_R1.fastq",
              f"{container_sample_path}/HLA_R2.fastq",
        "--dna",
        "--verbose",
        "-o", f"{container_sample_path}/OptiType",
        "--config", "/config.ini"
    ], f"running OptiType for {sample_name}")