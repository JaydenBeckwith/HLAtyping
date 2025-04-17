import os
import subprocess
from pathlib import Path
import argparse
from Bio import SeqIO

def run_command(cmd, description=""):
    print(f"\n[INFO] {description}")
    try:
        subprocess.run(cmd, check=True)
        print("[SUCCESS]")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {description}")
        print(e)
        exit(1)

def bwa_index_if_missing(genome_fasta):
    index_files = [f"{genome_fasta}.{ext}" for ext in ["bwt", "pac", "ann", "amb", "sa"]]
    if not all(Path(f).exists() for f in index_files):
        run_command([
            "docker", "run", "--rm",
            "-v", f"{Path(genome_fasta).parent}:/ref",
            "biocontainers/bwa:v0.7.17_cv1",
            "bwa", "index", f"/ref/{Path(genome_fasta).name}"
        ], f"BWA index for {Path(genome_fasta).name}")
    else:
        print("[INFO] BWA index found.")

def trim_and_qc(fq1, fq2, sample_name, output_dir):
        trimmed_dir = output_dir / "Trimmed"
        trimmed_dir.mkdir(exist_ok=True)

        # Run FastQC on raw reads
        run_command([
            "docker", "run", "--rm", "-v", f"{output_dir}:/data",
            "staphb/fastqc:0.11.9", "fastqc", "-o", "/data",
            f"/data/{fq1.name}", f"/data/{fq2.name}"
        ], f"FastQC before trimming for {sample_name}")

        # Run Trim Galore
        run_command([
            "docker", "run", "--rm", "-v", f"{output_dir}:/data",
            "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0",
            "trim_galore", "--paired", "-o", "/data/Trimmed",
            f"/data/{fq1.name}", f"/data/{fq2.name}"
        ], f"Trim Galore for {sample_name}")

        # Identify trimmed output files
        r1_trimmed = next(trimmed_dir.glob("*_val_1.fq.gz"))
        r2_trimmed = next(trimmed_dir.glob("*_val_2.fq.gz"))

        # Run FastQC on trimmed reads
        run_command([
            "docker", "run", "--rm", "-v", f"{output_dir}:/data",
            "staphb/fastqc:0.11.9", "fastqc", "-o", "/data",
            f"/data/Trimmed/{r1_trimmed.name}", f"/data/Trimmed/{r2_trimmed.name}"
        ], f"FastQC after trimming for {sample_name}")

        return r1_trimmed, r2_trimmed


def align_and_index(fq1, fq2, sample_name, genome_fasta, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    sam_path = output_dir / f"{sample_name}.sam"
    sorted_bam = output_dir / f"{sample_name}.sorted.bam"
    dedup_bam = output_dir / f"{sample_name}.dedup.bam"

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "-v", f"{genome_fasta.parent}:/ref",
        "--entrypoint", "bash", "biocontainers/bwa:v0.7.17_cv1", "-c",
        f"bwa mem -t 4 /ref/{genome_fasta.name} /data/{fq1.name} /data/{fq2.name} > /data/{sam_path.name}"
    ], f"Align {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "view", "-bS",
        f"/data/{sam_path.name}", "-o", f"/data/{sample_name}.bam"
    ], f"SAM to BAM for {sample_name}")

    # Sort BAM
    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "sort",
        f"/data/{sample_name}.bam", "-o", f"/data/{sorted_bam.name}"
    ], f"Sort BAM for {sample_name}")

    # Index sorted BAM
    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "index",
        f"/data/{sorted_bam.name}"
    ], f"Index sorted BAM for {sample_name}")

    # Mark duplicates with Picard
    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "broadinstitute/picard",
        "picard", "MarkDuplicates",
        f"I=/data/{sorted_bam.name}",
        f"O=/data/{dedup_bam.name}",
        f"M=/data/{sample_name}.metrics.txt",
        "CREATE_INDEX=true", "VALIDATION_STRINGENCY=LENIENT"
    ], f"Mark duplicates for {sample_name}")

    return dedup_bam

def run_mutect2(tumor_bam, normal_bam, sample_id, genome_fasta, output_dir):
    raw_vcf = output_dir / f"{sample_id}_raw.vcf.gz"
    filtered_vcf = output_dir / f"{sample_id}_filtered.vcf.gz"

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "-v", f"{genome_fasta.parent}:/ref",
        "broadinstitute/gatk", "gatk", "Mutect2",
        "-R", f"/ref/{genome_fasta.name}",
        "-I", f"/data/{tumor_bam.name}", "-tumor", "TUMOR",
        "-I", f"/data/{normal_bam.name}", "-normal", "NORMAL",
        "-O", f"/data/{raw_vcf.name}"
    ], "Run Mutect2")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "-v", f"{genome_fasta.parent}:/ref",
        "broadinstitute/gatk", "gatk", "FilterMutectCalls",
        "-R", f"/ref/{genome_fasta.name}",
        "-V", f"/data/{raw_vcf.name}",
        "-O", f"/data/{filtered_vcf.name}"
    ], "Filter Mutect2 calls")

    return filtered_vcf

def run_vep(filtered_vcf, output_dir):
    vep_vcf = output_dir / "annotated.vcf"
    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "ensemblorg/ensembl-vep",
        "vep", "-i", f"/data/{filtered_vcf.name}", "-o", f"/data/{vep_vcf.name}",
        "--vcf", "--cache", "--offline", "--assembly", "GRCh38",
        "--symbol", "--canonical", "--distance", "5"
    ], "VEP annotation")
    return vep_vcf

def extract_peptides_from_vep(vep_vcf, output_fasta):
    print("[INFO] Extracting peptides from VEP VCF...")
    peptides = []
    with open(vep_vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            info = fields[7]
            if "missense_variant" in info:
                peptides.append(f"MUT{len(peptides)}")  # Placeholder

    with open(output_fasta, "w") as f_out:
        for i, p in enumerate(peptides):
            f_out.write(f">mut{i}\n{p}\n")

    print(f"[SUCCESS] Wrote peptides to {output_fasta}")
    return output_fasta

def run_mhcflurry(peptide_fasta, hla_list, output_dir):
    alleles_file = output_dir / "alleles.txt"
    peptides_file = output_dir / "peptides.txt"

    # Convert and write alleles
    with open(alleles_file, "w") as f:
        for allele in hla_list:
            formatted = allele.replace("HLA-", "").replace("*", "").replace(":", "")
            f.write(formatted + "")

    # Convert FASTA to plain peptide list
    with open(peptide_fasta) as fasta, open(peptides_file, "w") as txt:
        for record in SeqIO.parse(fasta, "fasta"):
            txt.write(str(record.seq) + "")

    run_command([
        "docker", "run", "--rm", "-v", f"{peptide_fasta.parent}:/data",
        "openvax/mhcflurry",
        "mhcflurry-predict",
        "--alleles-file", "/data/alleles.txt",
        "--peptides", "/data/peptides.txt",
        "--out", "/data/mhcflurry_predictions.csv"
    ], "MHCflurry binding prediction")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tumor_r1", required=True)
    parser.add_argument("--tumor_r2", required=True)
    parser.add_argument("--normal_r1", required=True)
    parser.add_argument("--normal_r2", required=True)
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--hla", nargs="+", required=True)
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    genome_fasta = Path(args.genome_fasta).resolve()

    tumor_r1 = Path(args.tumor_r1)
    tumor_r2 = Path(args.tumor_r2)
    normal_r1 = Path(args.normal_r1)
    normal_r2 = Path(args.normal_r2)

    bwa_index_if_missing(genome_fasta)

    tumor_r1_trimmed, tumor_r2_trimmed = trim_and_qc(tumor_r1, tumor_r2, "tumor", output_dir)
    tumor_bam = align_and_index(tumor_r1_trimmed, tumor_r2_trimmed, "tumor", genome_fasta, output_dir)
    normal_r1_trimmed, normal_r2_trimmed = trim_and_qc(normal_r1, normal_r2, "normal", output_dir)
    normal_bam = align_and_index(normal_r1_trimmed, normal_r2_trimmed, "normal", genome_fasta, output_dir)

    filtered_vcf = run_mutect2(tumor_bam, normal_bam, "sample", genome_fasta, output_dir)
    vep_vcf = run_vep(filtered_vcf, output_dir)

    peptides_fasta = extract_peptides_from_vep(vep_vcf, output_dir / "peptides.fasta")
    run_mhcflurry(peptides_fasta, args.hla, output_dir)

if __name__ == "__main__":
    main()
