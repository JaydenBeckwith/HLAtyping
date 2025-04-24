#!/bin/bash

# pull_dockers.sh
# Usage: bash pull_dockers.sh

set -e

echo "Pulling Docker containers..."

# HLA Typing
docker pull biocontainers/bwa:v0.7.17_cv1
docker pull biocontainers/samtools:v1.3.1_cv4
docker pull staphb/fastqc:0.11.9
docker pull quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0
docker pull fred2/optitype

# Neoantigen Pipeline
docker pull broadinstitute/picard
docker pull broadinstitute/gatk
docker pull ensemblorg/ensembl-vep
docker pull griffithlab/pvactools

echo "Done pulling all required containers."
