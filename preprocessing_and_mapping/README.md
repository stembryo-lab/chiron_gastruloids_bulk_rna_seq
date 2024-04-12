# PREPROCESSING AND MAPPING 
The structure of the pipeline is:

- `requirements.txt`: software and versions used in the pipleline.
- `download_data.sh`: download data used in the pipeline.
- `fastqc_and_trimming.sh`: QC and trimming of raw fastqs.
- `mapping.sh`: mapping reads to genome.
- `countdata_construct.ipynb`: mapped reads concatenation to construct bulk countdata matrix.

## 0. Download data and references
Firstly, run:
'''
./download_data.sh
'''
The script downloads the raw sequencing fastq reads and metadata used in the project, as well as the genome v25 (.fasta) and annotation (.gtf) used to align the reads.
Folder structure created needs to be preserved for the downstream analysis.

## 1. FastQC and Trimming
To perform the quality control with FastQC/MultiQC and adapter trimming with TrimGalore! on the reads, run `fastqc_and_trimming.sh`.
For it, **four** arguments need to be provided to run the pipeline (space separated).

1. `ANALYSIS_NAME`: provide a name for the analysis that will be used in the report name.
2. `true/false`: run FastQC on the raw fastq files (set to `true`) or not (set to `false`).
3. `true/false`: run TrimGalore! on the raw fastq files (set to `true`) or not (set to `false`).
4. `true/false`: run FastQC on the trimmed fastq files (set to `true`) or not (set to `false`).
   
For example (in this case, run all the processes):
'''
./fastqc_and_trimming.sh chir_gastruloids true true true
'''

## 2. Mapping
To align the trimmed reads with to the genome with STAR, run `mapping.sh`
In this case, **two** arguments need to be provided.

1. `ANALYSIS_NAME`: provide a name for the analysis that will be used in the report name.
2. `true/false`: perform the alignment of the reads (set to `true`) or not (set to `false`).

For example (in this case, do align reads):
'''
./fastqc_and_trimming.sh chir_gastruloids true 
'''
If willing to align the raw reads (default are trimmed reads), modify the `$PATH2READS` path to *data/fastqs/raw*.

## 3. Countdata construct
Execute the `countdata_construct.ipynb` script to concatenate mapped read counts and construct the bulk countdata object (which is used for the downstream analysis pipeline).

