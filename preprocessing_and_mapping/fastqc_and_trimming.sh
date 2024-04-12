#!/bin/bash

#SBATCH -J qc_trim # A job-name
#SBATCH -n 1 # number of cores
#SBATCH -p haswell # Partition 
#SBATCH --mem 16000 # Memory request (8Gb) 
#SBATCH -o ./jobs.out/slurm.%j.out
#SBATCH -e ./jobs.out/slurm.%j.err

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4    

# input parameters
DATASET_NAME=$1
DO_FASTQC_RAW=$2 #perform FastQC on raw fastqs
DO_TRIMMING=$3 # perform trimming
DO_FASTQC_TRIMMED=$4 #perform FastQC on trimmed fastqs

# paths 
PATH2RAW='data/fastqs/raw' # path to raw files
PATH2TRIMMED='data/fastqs/trimmed' # path to trimmed files
PATH2FASTQC='data/fastqs/fastqc' # path to FastQC reports

# parameters
TRIM_Q=20 # read quality threshold
TRIM_TYPE='--illumina'

cat << \
.
            QC & TRIMMING
=========================================

# PIPELINE GLOBAL OPTIONS
- Dataset Name: $DATASET_NAME
- Perform FastQC on raw files: $DO_FASTQC_RAW
- Perform Trimming: $D0_TRIMMING
- Perform FastQC on trimmed files: $DO_FASTQC_TRIMMED

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

if [ $# -ne 4 ]
then
    cat << \
.
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    echo "2) "true/false" for performing FastQC on raw files"   
    echo "3) "true/false" for performing Trimming"
    echo "4) "true/false" for performing FastQC on trimmed files"   
.
fi

: '
## QC ##

Run FastQC for every file in every condition. 
Then, MultiQC merges all FastQC reports.
'
# FastQC
mkdir -p ${PATH2FASTQC}

if [ "$DO_FASTQC_RAW" = true ]
then
    echo '
    PERFORMING FASTQC ANALYSIS
    Starting FastQC analysis...
    '
    for file in ${PATH2RAW}/*.fastq.gz
    do
        sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

        for file in ${PATH2RAW}/*.fastq.gz
        do
            fastqc -o $PATH2FASTQC $file
            echo 'Analysis on '${sample}' done'
        done
    done
    echo '
    ==> FASTQC ANALYSIS FINISHED <==
    --------------------------------

    Number of files analyzed: '$(ls ${PATH2RAW} -1 | wc -l)'
    Reports stored in: '${PATH2FASTQC}'

    '
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

elif [ "$DO_FASTQC_RAW" = false ]
then
    echo 'No FastQC performed on raw fastqs'
fi

module unload FastQC/0.11.9-Java-11
module unload MultiQC/1.9-foss-2019b-Python-3.7.4   

: '
## TRIMMING ##

Trim Illumina Adapter from read using TrimGalore!
'
if [ "$DO_TRIMMING" = true ]
then
    module load Trim_Galore/0.6.7-GCCcore-10.3.0 
    echo '
    PERFORMING ADAPTER TRIMMING
    Running Trim Galore!...
    '
    for file in ${PATH2RAW}/*R1*.fastq.gz
    do
        sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1 | cut -d '_' -f 1-2)

        # single-end samples (no R2)
        if [[ ! -f $PATH2RAW/${sample}_R2_001.fastq.gz ]]
        then
            trim_galore --quality $TRIM_Q $TRIM_TYPE $PATH2RAW/${sample}_R1_001.fastq.gz -o $PATH2TRIMMED
            echo 'Trimming for '${sample} ' done'
    
        elif [[ -f $PATH2RAW/${sample}_R2_001.fastq.gz ]]
        then
            trim_galore --quality $TRIM_Q $TRIM_TYPE --paired $PATH2RAW/${sample}_R1_001.fastq.gz $PATH2RAW/${sample}_R2_001.fastq.gz -o $PATH2TRIMMED
            echo 'Trimming for '${sample} ' done'
        fi
    done
    module unload Trim_Galore/0.6.7-GCCcore-10.3.0
    echo '
    ==> ADAPTER TRIMMING FINISHED <==
    --------------------------------
    Number of files analyzed: '$(ls ${PATH2RAW} -1 | wc -l)'
    Paired-end samples: '$(ls ${PATH2RAW}/*R2* -1 | wc -l)'
    Reports stored in: '${PATH2FASTQC}'
    '
elif [ "$DO_TRIMMING" = false ]
then
    echo 'No trimming performed on fastqs'
fi


if [ "$DO_FASTQC_TRIMMED" = true ]
then
    echo '
    PERFORMING FASTQC ANALYSIS
    Starting FastQC analysis...
    '
    for file in ${PATH2TRIMMED}/*.fq.gz
    do
        sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

        fastqc -o $PATH2FASTQC $file
        echo 'Analysis on '${sample}' done'

    done
    echo '
    ==> FASTQC ANALYSIS FINISHED <==
    --------------------------------

    Number of files analyzed: '$(ls ${PATH2TRIMMED} -1 | wc -l)'
    Reports stored in: '${PATH2FASTQC}'

    '
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
elif [ "$DO_FASTQC_TRIMMED" = false ]
then
    echo 'No FastQC performed on trimmed fastqs'
fi

