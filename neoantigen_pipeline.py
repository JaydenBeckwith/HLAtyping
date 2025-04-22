import os
import subprocess
from pathlib import Path
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip

def run_command(cmd, description=""):
    """Run a shell command and handle errors with logging."""
    print(f"\n[INFO] {description}")
    try:
        subprocess.run(cmd, check=True)
        print("[SUCCESS]")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] {description}")
        print(e)
        exit(1)

def bwa_index_if_missing(genome_fasta):
    """Run BWA index on the reference genome if index files are missing."""
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
    """Trim adapter sequences from FASTQ files and run FastQC before and after trimming."""
    trimmed_dir = output_dir / "Trimmed"
    trimmed_dir.mkdir(exist_ok=True)

    r1 = list(trimmed_dir.glob(f"*{sample_name}*_val_1.fq.gz"))
    r2 = list(trimmed_dir.glob(f"*{sample_name}*_val_2.fq.gz"))

    if r1 and r2:
        print(f"[INFO] Trimmed files for {sample_name} already exist. Skipping QC and trimming.")
        r1_trimmed, r2_trimmed = r1[0], r2[0]
    else:
        run_command([
            "docker", "run", "--rm",
            "-v", f"{fq1.parent.resolve()}:/input",
            "-v", f"{output_dir.resolve()}:/data",
            "staphb/fastqc:0.11.9", "fastqc", "-o", "/data",
            f"/input/{fq1.name}", f"/input/{fq2.name}"
        ], f"FastQC before trimming for {sample_name}")

        run_command([
            "docker", "run", "--rm",
            "-v", f"{fq1.parent.resolve()}:/input",
            "-v", f"{output_dir.resolve()}:/data",
            "quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0",
            "trim_galore", "--paired", "-o", "/data/Trimmed",
            f"/input/{fq1.name}", f"/input/{fq2.name}"
        ], f"Trim Galore for {sample_name}")

        r1_trimmed = next(trimmed_dir.glob("*_val_1.fq.gz"))
        r2_trimmed = next(trimmed_dir.glob("*_val_2.fq.gz"))

        run_command([
            "docker", "run", "--rm", "-v", f"{output_dir}:/data",
            "staphb/fastqc:0.11.9", "fastqc", "-o", "/data",
            f"/data/Trimmed/{r1_trimmed.name}", f"/data/Trimmed/{r2_trimmed.name}"
        ], f"FastQC after trimming for {sample_name}")

    return r1_trimmed, r2_trimmed

def align_and_index(fq1, fq2, sample_name, genome_fasta, output_dir):
    """Align reads with BWA, sort, add read groups, mark duplicates, and index BAM."""
    fq1 = output_dir / "Trimmed" / fq1.name
    fq2 = output_dir / "Trimmed" / fq2.name
    rg_bam = output_dir / f"{sample_name}.rg.bam"
    dedup_bam = output_dir / f"{sample_name}.dedup.bam"
    dedup_bai = output_dir / f"{sample_name}.dedup.bam.bai"

    if dedup_bam.exists() and dedup_bai.exists():
        print(f"[INFO] Alignment already completed for {sample_name}. Skipping.")
        return dedup_bam

    sam_path = output_dir / f"{sample_name}.sam"
    sorted_bam = output_dir / f"{sample_name}.sorted.bam"

    run_command([
        "docker", "run", "--rm",
        "-v", f"{output_dir}:/data",
        "-v", f"{genome_fasta.parent}:/ref",
        "--entrypoint", "bash", "biocontainers/bwa:v0.7.17_cv1", "-c",
        f"bwa mem -t 4 /ref/{genome_fasta.name} /data/Trimmed/{fq1.name} /data/Trimmed/{fq2.name} > /data/{sam_path.name}"
    ], f"Align {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "view", "-bS",
        f"/data/{sam_path.name}", "-o", f"/data/{sample_name}.bam"
    ], f"SAM to BAM for {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "sort",
        f"/data/{sample_name}.bam", "-o", f"/data/{sorted_bam.name}"
    ], f"Sort BAM for {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "broadinstitute/picard",
        "java", "-jar", "/usr/picard/picard.jar", "AddOrReplaceReadGroups",
        f"I=/data/{sorted_bam.name}", f"O=/data/{rg_bam.name}",
        "RGID=1", f"RGLB=lib1", f"RGPL=illumina", f"RGPU=unit1", f"RGSM={sample_name}",
        "VALIDATION_STRINGENCY=LENIENT"
    ], f"Add read groups for {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "broadinstitute/picard",
        "java", "-jar", "/usr/picard/picard.jar", "MarkDuplicates",
        f"I=/data/{rg_bam.name}", f"O=/data/{dedup_bam.name}",
        f"M=/data/{sample_name}.metrics.txt",
        "CREATE_INDEX=true", "VALIDATION_STRINGENCY=LENIENT"
    ], f"Mark duplicates for {sample_name}")

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data",
        "biocontainers/samtools:v1.3.1_cv4", "samtools", "index",
        f"/data/{dedup_bam.name}"
    ], f"Index dedup BAM for {sample_name}")

    return dedup_bam

def run_mutect2(tumor_bam, normal_bam, sample_id, genome_fasta, output_dir):
    """Run Mutect2 variant calling and filtering using GATK."""
    raw_vcf = output_dir / f"{sample_id}_raw.vcf.gz"
    filtered_vcf = output_dir / f"{sample_id}_filtered.vcf.gz"

    run_command([
        "docker", "run", "--rm", "-v", f"{output_dir}:/data", "-v", f"{genome_fasta.parent}:/ref",
        "broadinstitute/gatk", "gatk", "Mutect2",
        "-R", f"/ref/{genome_fasta.name}",
        "-I", f"/data/{tumor_bam.name}", "-tumor", "tumor",
        "-I", f"/data/{normal_bam.name}", "-normal", "normal",
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

def filter_vcf_standard_chroms(vcf_path: Path, output_path: Path):
    """Filter VCF (plain or .gz) to keep only standard chromosomes (1–22, X, Y, MT)."""
    standard_chroms = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}

    # Automatically handle .gz files
    open_func = gzip.open if vcf_path.suffix == ".gz" else open
    mode = "rt" if vcf_path.suffix == ".gz" else "r"

    with open_func(vcf_path, mode) as infile, output_path.open("w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                chrom = line.split()[0].replace("chr", "")
                if chrom in standard_chroms:
                    outfile.write(line)

    print(f"[INFO] Filtered VCF written to: {output_path}")

def run_vep(filtered_vcf_path: Path, output_dir: Path, vep_cache: Path, plugin_dir: Path, fasta_path: Path):
    """
    Run VEP with local cache, plugins (Wildtype & Frameshift), and reference FASTA.

    Args:
        filtered_vcf_path (Path): Path to the filtered VCF input file.
        output_dir (Path): Output directory where annotated VCF will be saved.
        vep_cache (Path): Path to the local VEP cache directory.
        plugin_dir (Path): Path to the directory containing VEP plugin files.
        fasta_path (Path): Path to the reference genome FASTA file (e.g., GRCh38.primary_assembly.genome.fa).
    """
    vep_vcf = output_dir / "annotated.vcf.gz"

    run_command([
        "docker", "run", "--rm",
        "-v", f"{filtered_vcf_path.parent.as_posix()}:/vcfdir",
        "-v", f"{output_dir.as_posix()}:/data3",
        "-v", f"{vep_cache.as_posix()}:/opt/vep/.vep",
        "-v", f"{plugin_dir.as_posix()}:/plugins",
        "-v", f"{fasta_path.parent.as_posix()}:/ref",
        "ensemblorg/ensembl-vep",
        "vep", "-i", f"/vcfdir/{filtered_vcf_path.name}",
               "-o", "/data3/annotated.vcf.gz",
        "--format", "vcf",
        "--vcf", "--verbose", "--assembly", "GRCh38",
        "--symbol", "--canonical", "--distance", "5",
        "--plugin", "Wildtype,/ref/" + fasta_path.name,
        "--plugin", "Frameshift,/ref/" + fasta_path.name,
        "--fasta", f"/ref/{fasta_path.name}",
        "--offline", "--cache",
        "--dir_plugins", "/plugins",
        "--dir_cache", "/opt/vep/.vep",
        "--tsl", "--biotype", "--hgvs",
        "--terms", "SO",
        "--force_overwrite",
        "--compress_output", "bgzip",
        "--fork", "4"
    ], "VEP annotation using local cache, plugins, and reference FASTA")

    print(f"[DEBUG] Files in output directory after VEP: {list(output_dir.iterdir())}")
    return vep_vcf

def run_pvacseq(sample_id: str, hla_list: list, output_dir: Path, threads: int):
    """
    Runs pVACseq using the annotated VCF in the output directory.
    """
    vcf_file = output_dir / "annotated.vcf.gz"
    if not vcf_file.exists():
        raise FileNotFoundError(f"Expected annotated VCF not found at: {vcf_file}")

    hla_str = ",".join(hla_list)

    run_command([
        "docker", "run", "--rm",
        "-v", f"{output_dir}:/data3",
        "griffithlab/pvactools",
        "pvacseq", "run",
        f"/data3/{vcf_file.name}", sample_id, hla_str,
        "MHCflurry", "/data3/pvacseq_output",
        "-t", str(threads)
    ], "Run pVACseq for neoantigen prediction")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tumor_r1")
    parser.add_argument("--tumor_r2")
    parser.add_argument("--normal_r1")
    parser.add_argument("--normal_r2")
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--hla", nargs="+", required=True)
    parser.add_argument("--sampleid", required=True)
    parser.add_argument("--vep_cache", required=True)
    parser.add_argument("--plugin_dir", required=True)
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for pVACseq")
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True) 

    genome_fasta = Path(args.genome_fasta).resolve()
    vep_cache = Path(args.vep_cache).resolve()
    plugin_dir = Path(args.plugin_dir).resolve()

    tumor_r1 = Path(args.tumor_r1)
    tumor_r2 = Path(args.tumor_r2)
    normal_r1 = Path(args.normal_r1)
    normal_r2 = Path(args.normal_r2)

    bwa_index_if_missing(genome_fasta)

    with ThreadPoolExecutor(max_workers=4) as executor:
        future_trim_tumor = executor.submit(trim_and_qc, tumor_r1, tumor_r2, "tumor", output_dir)
        future_trim_normal = executor.submit(trim_and_qc, normal_r1, normal_r2, "normal", output_dir)
        tumor_r1_trimmed, tumor_r2_trimmed = future_trim_tumor.result()
        normal_r1_trimmed, normal_r2_trimmed = future_trim_normal.result()

    with ThreadPoolExecutor(max_workers=4) as executor:
        future_align_tumor = executor.submit(align_and_index, tumor_r1_trimmed, tumor_r2_trimmed, "tumor", genome_fasta, output_dir)
        future_align_normal = executor.submit(align_and_index, normal_r1_trimmed, normal_r2_trimmed, "normal", genome_fasta, output_dir)
        tumor_bam = future_align_tumor.result()
        normal_bam = future_align_normal.result()

    filtered_vcf = run_mutect2(tumor_bam, normal_bam, args.sampleid, genome_fasta, output_dir)
    vcf_for_vep = output_dir / "filtered_for_vep.vcf"
    filter_vcf_standard_chroms(filtered_vcf, vcf_for_vep)
    run_vep(vcf_for_vep, output_dir, vep_cache, plugin_dir, genome_fasta)
    run_pvacseq(args.sampleid, args.hla, output_dir, args.threads)


if __name__ == "__main__":
    main()