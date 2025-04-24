import yaml

# Load config
with open("pipeline_config.yaml") as f:
    config = yaml.safe_load(f)

samples = config["samples"]
paths = config["paths"]
docker_images = config["docker_images"]

rule all:
  input:
    expand("results/neoantigen/{sample}/", sample=samples)

rule run_hla:
  input:
    fasta = paths["genome_fasta"],
    config_file = paths["optitype_config"],
    input_dir = lambda wc: f"data/hla/{wc.sample}"
  output:
    hla_file = "results/hla/{sample}_hla_alleles.txt"
  container:
    docker_images["hla_typing"]
  shell:
    """
    python3 hla_pipeline.py \
      --input_dir {input.input_dir} \
      --genome_fasta {input.fasta} \
      --optitype_config {input.config_file} > {output.hla_file}
    """

rule run_neo:
  input:
    tumor_r1 = lambda wc: f"data/neo/{wc.sample}_tumor_R1.fastq.gz",
    tumor_r2 = lambda wc: f"data/neo/{wc.sample}_tumor_R2.fastq.gz",
    normal_r1 = lambda wc: f"data/neo/{wc.sample}_normal_R1.fastq.gz",
    normal_r2 = lambda wc: f"data/neo/{wc.sample}_normal_R2.fastq.gz",
    fasta = paths["genome_fasta"],
    vep_cache = paths["vep_cache"],
    plugin_dir = paths["plugin_dir"],
    hla_file = "results/hla/{sample}_hla_alleles.txt"
  output:
    directory("results/neoantigen/{sample}/")
  container:
    docker_images["neoantigen_prediction"]
  threads: 4
  shell:
    """
    python3 neoantigen_pipeline.py \
      --tumor_r1 {input.tumor_r1} \
      --tumor_r2 {input.tumor_r2} \
      --normal_r1 {input.normal_r1} \
      --normal_r2 {input.normal_r2} \
      --genome_fasta {input.fasta} \
      --output_dir results/neoantigen/{wildcards.sample}/ \
      --hla $(cat {input.hla_file}) \
      --sampleid {wildcards.sample} \
      --vep_cache {input.vep_cache} \
      --plugin_dir {input.plugin_dir} \
      --threads {threads}
    """

