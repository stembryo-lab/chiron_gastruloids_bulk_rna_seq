#!/bin/bash

: 'QC and Adapter Trimming on fastqs'

# input parameters
DATASET_NAME=$1
DO_FASTQC=$2 #perform FastQC
D0_TRIM=$3 #perform Trimming

# paths 
PATH2RAW='data/fastqs/raw' # path to raw files
PATH2TRIMMED='data/fastqs/trimmed' # path to trimmed files
PATH2FASTQC='data/fastqs/fastqc' # path to FastQC reports

# parameters
TRIM_Q=20
TRIM_TYPE='--illumina'

cat << \
.
            QC & TRIMMING
=========================================

# PIPELINE GLOBAL OPTIONS
- Dataset Name: $DATASET_NAME
- Perform FastQC: $DO_FASTQC
- Perform Trimming: $D0_TRIM

# DATA FOLDER
Raw fastq folder: $PATH2RAW

# OUTPUT FOLDERS
Trimmed reads output folder: $PATH2TRIMMED
FastQC and MultiQC reports folder: $PATH2FASTQC

# TRIM GALORE PARAMETERS
Quality threshold: $TRIM_Q
Adapter type: $TRIM_TYPE

=========================================
.


if [ $# -ne 3 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    echo "2) "True/False" for performing FastQC"   
    echo "3) "True/False" for performing Trimming"
.
fi

: '
## QC ##

Run FastQC for every file in every condition. 
Then, MultiQC merges all FastQC reports.
'
# FastQC
mkdir -p ${PATH2FASTQC}

if [[ $DO_FASTQC==True ]]
then
    echo '
PERFORMING FASTQC ANALYSIS
Starting FastQC analysis...
    '
    for file in ${PATH2RAW}/*.fastq.gz
    do
        sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

        # perform FastQC 
        if [[ ! -f $PATH2FASTQC/${sample}_fastqc.html ]]
        then
            for file in ${PATH2RAW}/*.fastq.gz
            do
                fastqc -o $PATH2FASTQC $file
                echo 'Analysis on '${sample}' done'
            done
        else
            echo 'FastQC analysis already found for '${file}
        fi
    done
    echo '
==> FASTQC ANALYSIS FINISHED <==
--------------------------------

Number of files analyzed: '$(ls ${PATH2RAW} -1 | wc -l)'
Reports stored in: '${PATH2FASTQC}'

'
else
    echo 'No FastQC performed on raw fastqs'
fi

# MultiQC
echo 'MULTIQC
Merging FastQC reports...
'
multiqc ${PATH2FASTQC}/* -n ${DATASET_NAME}_raw -f -o $PATH2FASTQC 
echo '
==> MULTIQC ANALYSIS FINISHED <==
--------------------------------
Report name: '${DATASET_NAME}'.html
Report stored in: '${PATH2FASTQC} 

: '
## TRIMMING ##

Trim Illumina Adapter from read using TrimGalore!
'

# cd $SCRIPTS_DIR/$RAW_FILES_DIR

# for condition in *; do
#     for file in ${condition}/*R1*.fastq.gz; do
#         sample_id=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

#         #TRIMMING: Trim Galore (default parameters)

#         module load Trim_Galore/0.6.7-GCCcore-10.3.0 
 
#         if [[ $read_type == single ]]
#         then
#             cdna_read=${sample_id}.fastq.gz
#         elif [[ $read_type == paired ]]
#         then
#             cdna_read=${sample_id%R1_001}R2_001.fastq.gz
#         fi

#         # Run adapter trimmer on cDNA read
#         if [[ ! -f ${cdna_read%.fastq.gz}_trimmed.fq.gz  ]]; then
#             trim_galore ${condition}/${cdna_read} ${TRIM_TYPE} ${TRIM_Q} -o ${condition}
#         fi

#         module unload Trim_Galore/0.6.7-GCCcore-10.3.0
#         module load FastQC/0.11.9-Java-11
#     done
# done


# # FastQC on trimmed files
#  if [[ ! -f $FASTQC_OUT/${sample_id}_fastqc.html ]]; then
# for trimmed_file in *trimmed.fq.gz
# do
#     fastqc -o $FASTQC_OUT $trimmed_file
# done
# #load MultiQC
# module load MultiQC/1.9-foss-2019b-Python-3.7.4

# cd $FASTQC_OUT

# # Merge FastQC reports with MultiQC
# multiqc . -n $dataset_name -f -s
     

