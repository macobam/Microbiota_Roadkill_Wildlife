#!/bin/bash

#PROGRAM: QualityControlFastqFiles.sh 
#OBJETIVE: Assess the quality of FastQ files obtained from Illumina sequencing using the FastQC tool
#AUTHOR: Manuel Alejandro Coba-Males
#DATE: 2023/11/01

# User-provided paths for analysis
DataFolder=$1  # Path to the data folder which contains fastq files
ResultsQC_folder=$2  # Path to the folder to save FastQC results

# Export the path to the FastQC executable file
export fastqc=$HOME/Programs/FastQC/fastqc

# Check if the data directory exists
if [[ ! -d "$DataFolder" ]]; then
    echo -e "\nERROR: Directory $DataFolder does not exist."
    exit 1
fi

# Check if there are .fastq.gz files in the provided directory
if [[ $(find ${DataFolder} -maxdepth 1 -type f -iname "*.fastq.gz") ]]
then
	# Create the quali control results directory
	mkdir -p $ResultsQC_folder
	
	# Loop to run FastQC on each fastq file in the provided directory
	for fastqfile in ${DataFolder}*.fastq.gz
	do
		echo -e "\n***** Beginning fasta quality control of $(basename $fastqfile) *****"
		$fastqc ${fastqfile} --outdir=${ResultsQC_folder}
		echo "***** $(basename $fastqfile) processing complete *****"
	done
	
	# Print the number of files successfully executed with FastQC
	echo -e "\n##############################################"
	echo "Successfully executed FastQC on $(find $DataFolder -iname "*.fastq.gz" | wc -l) fastq files"
	echo "The results are available in: ${ResultsQC_folder}"
	echo -e "##############################################\n"
else
	echo -e "\nERROR: No fastq.gz files found in the specified directory."
	exit 1
fi
