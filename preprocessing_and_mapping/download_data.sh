#!/bin/bash

: ' This script downloads raw fastq files, metadata and references and stores them for the QC and Alignment pipeline'

# create folder to store data
mkdir -p data/fastqs/raw
mkdir -p data/references/{fasta,annotation}

# data repository: BioStudies ArrayExpress
FASTQS_URL=''
METADATA_URL=''
# chosen genome and annotation version: GENCODE M25
GENOME_URL='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz'
ANNOTATION_URL='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz'

# download 
#fastqs and metadata
wget ${FASTQS_URL} -P data/fastqs/raw/
wget ${METADATA_URL} -P data/fastqs/
# references
wget ${GENOME_URL} -P data/references/fasta/ 
wget ${ANNOTATION_URL} -P data/references/annotation/
