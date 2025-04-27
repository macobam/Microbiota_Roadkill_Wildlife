[![Workflow](https://img.shields.io/badge/Tool-Mothur-yellow?logo=data)](https://mothur.org/wiki/miseq_sop/)
![Analyzed by R](https://img.shields.io/badge/Analyzed%20by-R-blue?logo=r)
![Analyzed by Bash](https://img.shields.io/badge/Shell-Bash-lightgrey?logo=gnu-bash)
![License: CC BY 4.0](https://img.shields.io/badge/License-CC--BY%204.0-lightblue.svg)

# Wildlife Microbiota Roadkill Analysis

### This repository describes the necessary steps to analyze metabarcoding data from gut samples of roadkill wildlife in the Tropical Andes, specifically focusing on the species *Amphisbaena bassleri* and *Crotophaga ani*. The data was obtained through Illumina sequencing targeting the V<sub>3</sub> - V<sub>4</sub> region of the 16S rRNA gene.

---

## Considerations

* ***Note 1:*** The necessary FastQ files for this project should be placed within the `Data/` directory. You can download the FastQ files for this study from the Sequence Read Archive (SRA) under accession number [**PRJNA1061813**](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1061813/).
* ***Note 2:*** All scripts are located within `Scripts/`.
* ***Note 3:*** Verify and adjust the paths as needed for the environment variables according to your own directories for both batch files inside `Scripts/MothurBatchFiles/`.
* ***Note 4:*** To ensure the correct execution of all scripts, it is necessary to have the following programs installed:
  * **mothur v.1.45.0** (https://github.com/mothur/mothur/releases/tag/v.1.45.0)
  * **FastQC** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  * **Entrez Direct** (https://www.ncbi.nlm.nih.gov/books/NBK179288/)

**Create a project folder**
```
mkdir ~/Microbiota_Roadkill_Wildlife/
```

## 1. Sequencing Data Quality Control

This step performs a set of analyses on one or more raw fastq format sequence files. The script developed for this analysis requires two mandatory input arguments.

*Argument 1:* Path to the folder containing the fastq files: `path/to/the/folder/with/fastqfiles/`.

*Argument 2:* Path to the folder where the results will be saved: `path/to/the/results/folder/`.
```
bash QualityControlFastqFiles.sh ~/Microbiota_Roadkill_Wildlife/Data/ ~/Microbiota_Roadkill_Wildlife/QC_Results/
```

## 2. Set up the Working Environment for the Created Project

This step creates the necessary directories and downloads reference files for subsequent analyses.
```
bash WorkEnvironment.sh
```

## 3. Customize SILVA Reference Alignment

This step customizes the reference alignment of the SILVA non-redundant database v132 to retain only sequences corresponding to the V3-V4 region of the 16S rRNA gene. The mothur pipeline used in this stage can be found in `Scripts/MothurBatchFiles/customizeV3V4region.batch`.
```
bash CustomizeSILVAreference.sh
```

***Note:*** This step needs to be performed only once before the analysis.

## 4. Metabarcoding Data Analysis

This step analyzes the sequencing data to perform taxonomic assignment in each sample. The mothur pipeline used can be found in `Scripts/MothurBatchFiles/metabarcoding.batch`.
```
bash MetabarcodingAnalysis.sh
```

This generates multiple files in the `Results/Mothur/Analysis/` folder. However, only two important output files that will be used for subsequent analysis in R are:

* `Results/Mothur/MothurResults/OtuTable.shared`
* `Results/Mothur/MothurResults/TaxTable.taxonomy`

## 5. Metabarcoding Results Analysis

This final step requires Rstudio to analyze the results obtained from mothur. In `Scripts/MicrobiotaAnalysis/`, the script `MicrobiotaAnalysis.R` examines, for each species and sample: community composition, alpha diversity, beta diversity, and the core microbiota.

The main results are in `Results/Microbiota/Plots/`. Additionally, `Results/Microbiota/RData/` stores the data sets used, which can be imported into Rstudio using the `load("file_name.RData")` function.

## Citation

Please cite our paper if you find it useful for your work.

Coba-Males, M.A., Díaz M., Molina C.A., Medrano-Vizcaíno P., Brito-Zapata D., Martin-Solano S., Ocaña-Mayorga, S., Carrillo-Bilbao, G. A., Narváez, W., Arrivillaga-Henríquez, J., González-Suárez, M., Enríquez, S. , & Poveda, A. (2024) Gut bacterial communities in roadkill animals: A pioneering study of two species in the Amazon region in Ecuador. *PLoS ONE 19*(12): e0313263. [DOI: 10.1371/journal.pone.0313263](https://doi.org/10.1371/journal.pone.0313263)