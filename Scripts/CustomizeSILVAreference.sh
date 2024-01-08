#!/bin/bash

#PROGRAM: CustomizeSILVAreference.sh
#OBJETIVE: Customize the SILVA_NR_v132 refence alignment to retain only V3-V4 region information of 16S rRNA sequences
#AUTHOR: Manuel Alejandro Coba-Males
#DATE: 2023/11/07

# Path to the directories
CustomizeV3V4_SILVA=$HOME/Microbiota_Roadkill_Wildlife/Results/Mothur/CustomizeV3V4_SILVA/
MothurBatchFiles=$HOME/Microbiota_Roadkill_Wildlife/Scripts/MothurBatchFiles/

# Export the path to the mothur program folder
export PATH=$PATH:$HOME/Programs/mothur-1.45/

# Run in the batch mode the file for customizing the SILVA reference
mothur ${MothurBatchFiles}customizeV3V4region.batch 

# Rename the final customize SILVA reference alignment
mv ${CustomizeV3V4_SILVA}silva.nr_v132.pcr.align ${CustomizeV3V4_SILVA}customize_silva.nr_v132.align

# Rename and move the mothur output logfile
mv ./mothur.[0-9]*.logfile ${CustomizeV3V4_SILVA}CustomizeSILVA_v3v4.logfile
