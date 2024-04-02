#!/bin/bash

#SBATCH -J star # A job-name
#SBATCH -p haswell # Partition 
#SBATCH -n 12 # Number of cores 
#SBATCH --mem 50000 #(GB per node)
#SBATCH -o ./Jobs_info/slurm.%j.out
#SBATCH -e ./Jobs_info/slurm.%j.err

interactive
module load gzip/1.10-GCCcore-11.2.0 
module load STAR/2.7.10a-GCC-11.2.0

### input paths (to modify by user)

SCRIPTS_DIR='/homes/users/ppascual/scratch/Gulbenkian_André/Alignment' #relative to your path
RAW_FILES_DIR="./Data/Fastq" # path to where reads are stored

cd $SCRIPTS_DIR

### check input parameters
if [ $# -ne 1 ]
then
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    exit
fi

dataset_name=$1

# 0. GENOME INDEX GENERATION
# If not previously generated, aligner needs to use the reference genome and annotation file in order to generate an INDEX for the later alignment
GENOME_FILE='./Data/References/Mus_musculus/fasta/GRCm39.genome.fa.gz'
ANNOTATION_FILE='./Data/References/Mus_musculus/annotation/gencode.vM33.primary_assembly.basic.annotation.gtf.gz'
GENOME_INDEX_DIR='./Data/References/Mus_musculus/genome_index/'
mkdir $GENOME_INDEX_DIR


: ' ### READ ALIGNMENT: STAR + MultiQC ###
    1. Genome Index generation
    2. Results folder creation
    3. Read Mapping
    4. QC Reports
'

# 1. GENOME INDEX GENERATION
if [[ ! -f ${GENOME_INDEX_DIR}/genomeParameters.txt ]] 
then
    STAR \
    --runMode genomeGenerate \
    --runThreadN 12 \
    --genomeDir $GENOME_INDEX_DIR \
    --genomeFastaFiles $GENOME_FILE \
    --sjdbGTFfile $ANNOTATION_FILE \
    --genomeSAindexNbases 12 \
    --genomeSAsparseD 3 
fi 

#2. RESULTS FOLDER CREATION
# Create folders to store results per condition
cd $SCRIPTS_DIR/$RAW_FILES_DIR

STAR_OUT='STAR_Results'
for condition in XAV*h/; do
    condition_id=${condition::-1}  
    for file in ${condition_id}/*R1_001.fastq.gz; do
    
        replicate_id=$(echo $file | sed 's:.*/::')
        replicate_id=${replicate_id%_R1_001.fastq.gz}

        cd $SCRIPTS_DIR
        # create directory for each condition
        mkdir -p $STAR_OUT/${condition_id}/${replicate_id}/

        cd $RAW_FILES_DIR
    done
done

# 3. READ MAPPING
#Running STAR aligner
for condition in XAV*/; do
    condition_id=${condition::-1}  
    for file in ${condition_id}/AF*_R1_001.fastq.gz; do
        replicate_id=$(echo $file | sed 's:.*/::')
        replicate_id=${replicate_id%_R1_001.fastq.gz}

        # Input read files
        # single-end reads
        if [[ ! -f ${condition_id}/${replicate_id}_R2_001.fastq.gz ]]
        then
            r1=$(ls ${condition_id}/${replicate_id}_R1_001_trimmed.fq.gz)

            cd $SCRIPTS_DIR

            # runnning STARsolo aligner, parameters to resemble the most as possible as CellRanger
            if [[ ! -f $STAR_OUT/${condition_id}/${replicate_id}/Log.final.out ]]; then
                STAR \
                --genomeDir $GENOME_INDEX_DIR \
                --genomeLoad NoSharedMemory \
                --runThreadN 12 \
                --readFilesCommand zcat \
                --readFilesIn $RAW_FILES_DIR/${r1}  \
                --sjdbGTFfile $ANNOTATION_FILE \
                --sjdbGTFtagExonParentGene gene_name \
                --outFilterMultimapNmax 10 \
                --outFileNamePrefix $STAR_OUT/${condition_id}/${replicate_id}/ \
                --outSAMtype BAM Unsorted \
                --quantMode GeneCounts

                cd $RAW_FILES_DIR
            fi

        # paired-end reads    
        elif [[ -f ${condition_id}/${replicate_id}_R2_001.fastq.gz ]]
        then
            r1=$(ls ${condition_id}/${replicate_id}_R1_001.fastq.gz)
            r2=$(ls ${condition_id}/${replicate_id}_R2_001.fastq.gz)

            echo $r1 $r2
            cd $SCRIPTS_DIR

            # runnning STARsolo aligner, parameters to resemble the most as possible as CellRanger
            if [[ ! -f $STAR_OUT/${condition_id}/${replicate_id}/Log.final.out ]]; then
                STAR \
                --genomeDir $GENOME_INDEX_DIR \
                --genomeLoad NoSharedMemory \
                --runThreadN 12 \
                --readFilesCommand zcat \
                --readFilesIn $RAW_FILES_DIR/${r1} $RAW_FILES_DIR/${r2} \
                --sjdbGTFfile $ANNOTATION_FILE \
                --sjdbGTFtagExonParentGene gene_name \
                --outFilterMultimapNmax 10 \
                --outFileNamePrefix $STAR_OUT/${condition_id}/${replicate_id}/ \
                --outSAMtype BAM Unsorted \
                --quantMode GeneCounts

                cd $RAW_FILES_DIR
            fi

        elif [[ -f ${condition_id}/${replicate_id}_R2_001_val_2.fq.gz ]]
        then
            r1=$(ls ${condition_id}/${replicate_id}_R1_001.fastq.gz)
            r2=$(ls ${condition_id}/${replicate_id}_R2_001.fastq.gz)

            cd $SCRIPTS_DIR

            # runnning STARsolo aligner, parameters to resemble the most as possible as CellRanger
            if [[ ! -f $STAR_OUT/${condition_id}/${replicate_id}/Log.final.out ]]; then
                STAR \
                --genomeDir $GENOME_INDEX_DIR \
                --genomeLoad NoSharedMemory \
                --runThreadN 12 \
                --readFilesCommand zcat \
                --readFilesIn $RAW_FILES_DIR/${r1} $RAW_FILES_DIR/${r2} \
                --sjdbGTFfile $ANNOTATION_FILE \
                --sjdbGTFtagExonParentGene gene_name \
                --outFilterMultimapNmax 10 \
                --outFileNamePrefix $STAR_OUT/${condition_id}/${replicate_id}/ \
                --outSAMtype BAM Unsorted \
                --quantMode GeneCounts

                cd $RAW_FILES_DIR
            fi
            echo $r1
        fi
    done
done


####### 2. QUALITY CONTROL ON STARSOLO RESULTS: MULTIQC
cd $SCRIPTS_DIR/$STAR_OUT   

# create folder to store counts.csv and alignment reports
mkdir -p Gene_counts MultiQC

#QC on the alignment results (summary.csv file)
for condition in *h/; do
    condition_id=${condition::-1}  
    for replicate in ${condition_id}/*; do
        replicate_id=$(echo $replicate | sed 's:.*/::')

        # copying alignment reports files to folder 
        cp ${condition_id}/${replicate_id}/Log.final.out MultiQC/${condition_id}_${replicate_id}_Log.final.out
        # copying counts.csv to Gene_counts folder
        cp ${condition_id}/${replicate_id}/ReadsPerGene.out.tab Gene_counts/${condition_id}_${replicate_id}.tab

    done
done


#MultiQC 
module load MultiQC/1.9-foss-2019b-Python-3.7.4 

cd MultiQC/

#run MultiQC on alignment reports
multiqc . -n $dataset_name -f -s

# cd $SCRIPTS_DIR/Results/Alignment
# module unload MultiQC/1.9-foss-2019b-Python-3.7.4 


# # Sort BAM files and get raw counts file
# cd $SCRIPTS_DIR/Results/BAM_files
# module load  SAMtools/1.12-GCC-9.3.0
# samtools index -@ 8 *.bam


# for condition in */; do
#     condition_id=${condition::-1}
#     for cell_num in $condition/d*; do
#         cellnum_id=$(echo $cell_num | sed 's:.*/::')
#         for replicate in $cell_num/${condition_id}*; do
#             replicate_id=$(echo $replicate | sed 's:.*/::')

#             # sort BAM files with SAMtools
#             module load  SAMtools/1.12-GCC-9.3.0
#             samtools sort -@ 8 ${condition_id}/${cellnum_id}/${replicate_id}/${replicate_id}_${version}_Aligned.out.bam -o ${condition_id}/${cellnum_id}/${replicate_id}/sorted_${replicate_id}_${version}_Aligned.out.bam 
#             samtools index -@ 8 ${condition_id}/${cellnum_id}/${replicate_id}/sorted_${replicate_id}_${version}_Aligned.out.bam 
#             module unload SAMtools/1.12-GCC-9.3.0
#         done
#     done
# done

########################################################

: ' STARsolo configuration regarding each Chemistry version used

10x v1
Whitelist, 737K-april-2014_rc.txt
CB length, 14
UMI start, 15
UMI length, 10 

10X v2
Whitelist, 737K-august-2016.txt
CB length, 16
UMI start, 17
UMI length, 10

10x v3
Whitelist, 3M-Feb_2018_V3.txt
CB length, 16
UMI start, 17
UMI length, 12 '