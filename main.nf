// main.nf - Nextflow pipeline for ImmunoTools

params.samples = ["48354"]
params.config = "pipeline_config.yaml"

// Load YAML config
config = file(params.config).text
config_parsed = config.toYaml().loadYaml()

def genome_fasta = config_parsed.paths.genome_fasta
def optitype_config = config_parsed.paths.optitype_config

def vep_cache = config_parsed.paths.vep_cache
def plugin_dir = config_parsed.paths.plugin_dir

def hla_container = config_parsed.docker_images.hla_typing

def neo_container = config_parsed.docker_images.neoantigen_prediction

Channel.fromList(params.samples)
    .set { sample_ids }

process HLA_Typing {
    container hla_container
    input:
    val sample_id

    output:
    path "results/hla/${sample_id}_hla_alleles.txt" 

    script:
    """
    python3 hla_pipeline.py \
      --input_dir data/hla/${sample_id} \
      --genome_fasta ${genome_fasta} \
      --optitype_config ${optitype_config} > results/hla/${sample_id}_hla_alleles.txt
    """
}

process Neoantigen_Prediction {
    container neo_container
    input:
    val sample_id
    path "results/hla/${sample_id}_hla_alleles.txt"

    output:
    path "results/neoantigen/${sample_id}/"

    script:
    """
    python3 neoantigen_pipeline.py \
      --tumor_r1 data/neo/${sample_id}/${sample_id}_tumor_R1.fastq.gz \
      --tumor_r2 data/neo/${sample_id}/${sample_id}_tumor_R2.fastq.gz \
      --normal_r1 data/neo/${sample_id}/${sample_id}_normal_R1.fastq.gz \
      --normal_r2 data/neo/${sample_id}/${sample_id}_normal_R2.fastq.gz \
      --genome_fasta ${genome_fasta} \
      --output_dir results/neoantigen/${sample_id}/ \
      --hla \$(cat results/hla/${sample_id}_hla_alleles.txt) \
      --sampleid ${sample_id} \
      --vep_cache ${vep_cache} \
      --plugin_dir ${plugin_dir} \
      --threads 4
    """
}

workflow {
    sample_ids | HLA_Typing | Neoantigen_Prediction
}
