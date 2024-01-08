#!/bin/bash

#PROGRAM: WorkEnvironment.sh
#OBJETIVE: Prepare the work environment by creating the necessary folders and download files before running the metabarcoding analysis 
#AUTHOR: Manuel Alejandro Coba-Males
#DATE: 2023/11/06

# Paths to the project directory
Microbiota_Roadkill_Wildlife=$HOME/Microbiota_Roadkill_Wildlife/

# Paths to the reference directories
Results=${Microbiota_Roadkill_Wildlife}Results/
Microbiota=${Results}Microbiota/
Plots=${Microbiota}Plots/
RData=${Microbiota}RData/
Mothur=${Results}Mothur/
SILVA_NR_v132=${Mothur}SILVA_NR_v132/
CustomizeV3V4_SILVA=${Mothur}CustomizeV3V4_SILVA/
Analysis=${Mothur}Analysis/
MothurResults=${Mothur}MothurResults/

# Export the path to the edirect folder
export PATH=$PATH:$HOME/Programs/edirect/

# Create folders
mkdir -p \
	${Results} \
	${Microbiota} \
	${Plots} \
	${RData} \
	${Mothur} \
	${SILVA_NR_v132} \
	${CustomizeV3V4_SILVA} \
	${Analysis} \
	${MothurResults} \

# Download the oligo file containing primer sequences used for the V3-V4 region of the 16S rRNA gene for customize SILVA reference alignment
wget -P ${Mothur} https://gist.githubusercontent.com/macobam/f24b4cb098ac3176c80b0929633fc62f/raw/784d7ccd3c828b72e6436d2a4b767b1ccea00682/v3v4.oligos

# Download 16S rRNA gene sequence for E. coli (Accession No: NR_024570.1) from NCBI repository using Entrez
esearch -db nucleotide -query "NR_024570.1" | efetch -format fasta > ${Mothur}ecoli16S_refseq.fasta

#Download the reference alignment file of SILVA v132
wget -P ${SILVA_NR_v132} https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v132.tgz
tar xzvf ${SILVA_NR_v132}silva.nr_v132.tgz -C ${SILVA_NR_v132}
