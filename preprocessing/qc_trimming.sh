#!/bin/bash

: 'Performing QC and Trimming on the FASTQs before mapping'

module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4    

DATASET_NAME=$1

### check input parameters
if [ $# -ne 2 ]
then
    echo "Please, give:"
    echo "1) Naming for the dataset to analyze"   
    exit
fi

# paths to data
PATH2RAW='data/fastqs/raw'
PATH2TRIMMED='data/fastqs/trimmed'
PATH2FASTQC='data/fastqs/fastqc'

# parameters
TRIM_Q=20
TRIM_TYPE='--illumina'

: '
## FastQC + MultiQC ##

Run FastQC for every file in every condition. 
Then, MultiQC merges all FastQC reports.
'
# FastQC
for file in ${PATH2RAW}/*; do
    sample=$(echo $file | sed 's:.*/::' | cut -d '.' -f 1)

    if [[ ! -d $FASTQC_OUT ]]; then
        mkdir $FASTQC_OUT
    fi

    # perform FastQC
    if [[ ! -f $FASTQC_OUT/${sample_id}_fastqc.html ]]; then
        for file in ${condition}/*.fastq.gz; do
            fastqc -o $FASTQC_OUT $file
        done
    fi

    cd $SCRIPTS_DIR/$RAW_FILES_DIR
done

#perform MultiQC on FastQC files
cd $FASTQC_OUT

multiqc . -n $dataset_name -f

module unload FastQC/0.11.9-Java-11
module unload MultiQC/1.9-foss-2019b-Python-3.7.4   

: '
2. ADAPTER TRIMMING + QC ON TRIMMED FILES: TrimGalore! + FastQC + MultiQC

File containing cDNA (R1 in case of single-end; R2 in case of paired-end) is trimmed using Trim Galore with default parameters 
Once trimmed, FastQC and MultiQC is performed.
'

cd $SCRIPTS_DIR/$RAW_FILES_DIR

# for condition in *; do
#     for file in ${condition}/*R1*.fastq.gz; do
#         sample_id=$(echo $file | sed 's:.*/::')
#         sample_id=$(echo $sample_id | cut -d '.' -f 1)

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
#             trim_galore ${cdna_read} --illumina -q 20 
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
     

