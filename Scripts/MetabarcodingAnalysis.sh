#!/bin/bash

#PROGRAM: MetabarcodingAnalysis.sh
#OBJETIVE: Analyze metabarcoding data from gut samples in Amphisbaena bassleri and Crotophaga ani
#AUTHOR: Manuel Alejandro Coba-Males
#DATE: 2023/11/07

# Path to the directories
Analysis=$HOME/Microbiota_Roadkill_Wildlife/Results/Mothur/Analysis/
MothurResults=$HOME/Microbiota_Roadkill_Wildlife/Results/Mothur/MothurResults/
MothurBatchFiles=$HOME/Microbiota_Roadkill_Wildlife/Scripts/MothurBatchFiles/

# Export the path to the mothur program folder
export PATH=$PATH:$HOME/Programs/mothur-1.45/

# Run in the batch mode the file for customizing the SILVA database
mothur ${MothurBatchFiles}metabarcoding.batch

# Rename the output mothur files
mv ${Analysis}roadkill.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.01.cons.taxonomy ${MothurResults}TaxTable.taxonomy
mv ${Analysis}roadkill.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared ${MothurResults}OtuTable.shared

# Rename and move the output mothur logfile
mv ./mothur.[0-9]*.logfile ${Analysis}MetabarcodingAnalysis.logfile
