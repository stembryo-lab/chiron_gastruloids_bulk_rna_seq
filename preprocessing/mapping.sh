#!/bin/bash

#SBATCH -J star # A job-name
#SBATCH -p haswell # Partition 
#SBATCH -n 12 # Number of cores 
#SBATCH --mem 50000 #(GB per node)
#SBATCH -o ./jobs.out/slurm.%j.out
#SBATCH -e ./jobs.out/slurm.%j.err

: 'Mapping trasncriptomic reads to genome with STAR'

interactive
module load gzip/1.10-GCCcore-11.2.0 
module load STAR/2.7.10a-GCC-11.2.0

# global options
DATASET_NAME=$1
DO_MAPPING=$2

# paths
PATH2READS='data/fastqs/trimmed'
PATH2FASTA='data/references/fasta/GRCm38.p6.genome.fa'
PATH2ANNOTATION='data/references/annotation/gencode.vM25.primary_assembly.annotation.gtf'
PATH2GENOME_IDX='data/references/fasta/index'

MAPPING_OUT='star_out'
MULTIQC_OUT='star_out/reports'

# STAR parameters
THREADS=12
# index generation
GENOME_SA_INDEXBASES=12
GENOME_SA_SPARSE=3
#mapping
LOAD_MODE='NoSharedMemory'
READ_FORMAT='zcat'
GENE_ID='gene_name'
FILTER_MULTIMAP=10
OUTSAM_FORMAT='BAM Unsorted'
QUANTMODE='GeneCounts' 

cat << \
.
            MAPPING
=========================================

# PIPELINE GLOBAL OPTIONS
- Dataset Name: $DATASET_NAME
- Perform alignment: $DO_MAPPING

# DATA FOLDER
Raw fastq folder: $PATH2READS

# REFERENCES
Fasta: $(echo $PATH2FASTA | sed 's:.*/::')
Anntotation: $(echo $PATH2ANNOTATION | sed 's:.*/::')

# OUTPUT FOLDERS
Mapping output: $MAPPING_OUT
MultiQC report: $MULTIQC_OUT

# MAPPING PARAMETERS
# Memory
Threads: $THREADS
# memory for index generation
SA pre-indexingn length: $GENOME_SA_INDEXBASES
Indexing distance : $GENOME_SA_SPARSE
#mapping
Genome loading mode: $LOAD_MODE
Read format= $READ_FORMAT
Gene id in annotation: $GENE_ID
Multimappers filter: $FILTER_MULTIMAP
SAM/BAM Output format: $OUTSAM_FORMAT 
STAR quantMode: $QUANTMODE  


=========================================
.

if [ $# -ne 2 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    echo "2) "True/False" for performing Star Alignment"
.
fi


# 1. GENOME INDEX GENERATION
if [[ ! -f ${PATH2GENOME_IDX}/genomeParameters.txt ]] 
then
    mkdir -p ${PATH2GENOME_IDX}
    echo '
    Starting to index the genome...
    '

    STAR \
    --runMode genomeGenerate \
    --runThreadN $THREADS \
    --genomeDir $PATH2GENOME_IDX \
    --genomeFastaFiles $PATH2FASTA \
    --sjdbGTFfile $PATH2ANNOTATION \
    --genomeSAindexNbases $GENOME_SA_INDEXBASES \
    --genomeSAsparseD $GENOME_SA_SPARSE
else
    echo 'Genome index already computed.'
fi 


# 2.ALIGNMENT
for file in $PATH2READS/*.fq.gz
do
    sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1 | cut -d '_' -f 1-2)
    echo $sample
    if [[ ! -f $PATH2READS/${sample}_R1_001_val_1.fq.gz ]]
    r1=$PATH2READS/${sample}_R1_001_trimmed.fq.gz
    then
        if [[ ! -f $STAR_OUT/${condition_id}/${replicate_id}/Log.final.out ]]; then
            STAR \
            --genomeDir $PATH2GENOME_IDX \
            --genomeLoad $LOAD_MODE \
            --runThreadN $THREADS \
            --readFilesCommand $READ_FORMAT \
            --readFilesIn $r1  \
            --sjdbGTFfile $PATH2ANNOTATION \
            --sjdbGTFtagExonParentGene $GENE_ID \
            --outFilterMultimapNmax $FILTER_MULTIMAP \
            --outFileNamePrefix $MAPPING_OUT/$sample/ \
            --outSAMtype $OUTSAM_FORMAT \
            --quantMode $QUANTMODE

        fi

        # paired-end reads    

    elif [[ -f $PATH2READS/${sample}_R1_001_val_1.fq.gz ]]
    then
        r1=$PATH2READS/${sample}_R1_001_val_1.fq.gz
        r2=$PATH2READS/${sample}_R2_001_val_2.fq.gz

        if [[ ! -f $STAR_OUT/${condition_id}/${replicate_id}/Log.final.out ]]; then
            STAR \
            --genomeDir $PATH2GENOME_IDX \
            --genomeLoad $LOAD_MODE \
            --runThreadN $THREADS \
            --readFilesCommand $READ_FORMAT \
            --readFilesIn $r1 $r2 \
            --sjdbGTFfile $PATH2ANNOTATION \
            --sjdbGTFtagExonParentGene $GENE_ID \
            --outFilterMultimapNmax $FILTER_MULTIMAP \
            --outFileNamePrefix $MAPPING_OUT/$sample/ \
            --outSAMtype $OUTSAM_FORMAT \
            --quantMode $QUANTMODE

        fi
    fi
done


# ####### 2. QUALITY CONTROL ON STARSOLO RESULTS: MULTIQC
# cd $SCRIPTS_DIR/$STAR_OUT   

# # create folder to store counts.csv and alignment reports
# mkdir -p Gene_counts MultiQC

# #QC on the alignment results (summary.csv file)
# for condition in *h/; do
#     condition_id=${condition::-1}  
#     for replicate in ${condition_id}/*; do
#         replicate_id=$(echo $replicate | sed 's:.*/::')

#         # copying alignment reports files to folder 
#         cp ${condition_id}/${replicate_id}/Log.final.out MultiQC/${condition_id}_${replicate_id}_Log.final.out
#         # copying counts.csv to Gene_counts folder
#         cp ${condition_id}/${replicate_id}/ReadsPerGene.out.tab Gene_counts/${condition_id}_${replicate_id}.tab

#     done
# done


# #MultiQC 
# module load MultiQC/1.9-foss-2019b-Python-3.7.4 

# cd MultiQC/

# #run MultiQC on alignment reports
# multiqc . -n $dataset_name -f -s

