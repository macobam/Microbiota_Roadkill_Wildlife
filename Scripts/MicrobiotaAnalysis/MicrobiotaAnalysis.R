#Title: Microbiota Analysis in Roadkill Wildlife
#Author: Manuel Alejandro Coba-Males
#Date: 2023/12/06

################################################################################
# Some of the microbial diversity analyses are based on the R programming code 
# used by (Díaz et al., 2023) in https://doi.org/10.3389/fmicb.2023.1154815
# This served as a methodological foundation for adjusting and customizing 
# to enable the analysis of our own dataset.
###############################################################################

# Loading libraries ----

library(ComplexHeatmap)
library(cowplot)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggtext)
library(grid)
library(gridExtra)
library(MASS)
library(metagMisc)
library(microbiome)
library(microbiomeutilities)
library(phyloseq)
library(ranacapa)
library(RColorBrewer)
library(readxl)
library(sfsmisc)
library(stats)
library(stringr)
library(tidyverse)
library(VennDiagram)
library(viridis)

################################################################################ 
# IMPORTANT NOTE: Please consider the full paths to the input and output files.
################################################################################

# Set working directory ----

setwd("C:/Users/macobam/Desktop/Microbiota_Roadkill_Wildlife")

# Data importing from mothur results ----

# Outputs files names

sharedfile <- "./Results/Mothur/MothurResults/OtuTable.shared"

taxfile <- "./Results/Mothur/MothurResults/TaxTable.taxonomy"

# Creating a phyloseq object

mothur_raw_results <- import_mothur(mothur_shared_file = sharedfile,
                                    mothur_constaxonomy_file = taxfile)

# Adding taxonomy information to phyloseq object

colnames(tax_table(mothur_raw_results)) <- c("Kingdom", 
                                             "Phylum", 
                                             "Class", 
                                             "Order", 
                                             "Family", 
                                             "Genus")

# Importing metadata

metadata <- read_xlsx("./Scripts/MicrobiotaAnalysis/MetadataRoadkillWildlife.xlsx", 
                      sheet = "Metadata") %>%
  as.data.frame()

rownames(metadata) <- metadata$SampleID  

# Adding metadata information to phyloseq object

sample_data(mothur_raw_results) <- sample_data(metadata)

# Pre-processing of mothur raw results ----

# Agglomerate OTU to genus

mothur_results <- tax_glom(mothur_raw_results, "Genus")

mothur_results

# Remove other kingdoms

table(tax_table(mothur_results)[,1])

mothur_results <- subset_taxa(mothur_results, 
                              Kingdom != "unknown" & Kingdom != "Archaea")

mothur_results

# Remove unclassified bacteria

table(tax_table(mothur_results)[,2])

mothur_results <- subset_taxa(mothur_results, 
                              Phylum != "Bacteria_unclassified")

mothur_results

# Remove singletons

sum(taxa_sums(mothur_results) == 1)

mothur_results <- prune_taxa(taxa_sums(mothur_results) > 1, 
                             mothur_results)

mothur_results

# Sequencing depth ----

# Number of reads per sample

SeqDepth <- colSums(otu_table(mothur_results))

SeqDepth 

# Adding sequencing depth data to metadata in the phyloseq object

sample_data(mothur_results)$SeqDepth <- SeqDepth

# Bar plot of sequencing depth per sample

plot_seqdepth <- 
  ggpubr::ggbarplot(meta(mothur_results), 
                    x = "SampleID", 
                    y = "SeqDepth", 
                    fill = "Specie",
                    color = NA) +
  theme_bw() + 
  ggtitle("Sequencing depth per sample") +
  xlab("Sample ID") +
  ylab("Number of reads") +
  rotate_x_text(90) +  
  scale_y_continuous(expand = c(0,0),
                     limits = c(0, 190000)) + 
  scale_fill_aaas() +
  facet_wrap(~ Specie,
             scales = "free_x",
             nrow = 1) +
  theme(plot.title=element_text(face = "bold",
                                size = 15,
                                hjust = 0.5,
                                vjust = 0.5,
                                margin = margin(0, 0, 0.5, 0, "cm")),
        strip.text = element_text(face = "bold.italic",
                                  size = 10),
        strip.background = element_rect(fill = "#F0FFF0"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold",
                                    size = 10),
        legend.title.align = 0.5,
        legend.text = element_text(face = "italic",
                                   size = 10),
        axis.text.x = element_text(face = "bold", 
                                   size = 10),
        axis.text.y = element_text(face = "bold", 
                                   size = 10))

plot_seqdepth

# Stats sequencing depth by species

meta(mothur_results) %>% 
  group_by(Specie) %>% 
  summarise(mean = mean(SeqDepth), 
            sd = sd(SeqDepth), 
            median = median(SeqDepth), 
            Total = sum(SeqDepth),
            .groups = "drop")

# Separating samples per species ----

Species <- phyloseq_sep_variable(mothur_results, "Specie")

Amphisbaena_results <- Species$`Amphisbaena bassleri`

Crotophaga_results <- Species$`Crotophaga ani`

# Function to count the total number of OTUs within a specific higher taxonomic level from a deeper taxonomic level ----

count_taxlevel <- function (phyloseq_obj, HigherTaxLevel, DeeperTaxLevel) {
  tax_table <- tax_table(phyloseq_obj)
  
  tax_df <- as.data.frame(tax_table) %>%
    rownames_to_column(var = "OTU_ID") %>%
    filter(!is.na(.data[[HigherTaxLevel]]))
  
  count_per_taxlevel <- tax_df %>%
    group_by(.data[[HigherTaxLevel]]) %>%
    summarize(total = n_distinct(.data[[DeeperTaxLevel]]), .groups = 'drop') %>% 
    arrange(desc(total))
  
  print(count_per_taxlevel)
}

# Functions for Relative Abundance Analysis ----

# Create a table of OTUs count from phyloseq object

otu_counts <- function(phyloseq_object){
  otu_table(phyloseq_object) %>%
    psmelt() %>%
    as_tibble() %>%
    arrange(Sample, OTU) %>%
    rename(SampleID = Sample, Count = Abundance) %>%
    dplyr::select(SampleID, OTU, Count)
}

# Create a table of taxonomy from phyloseq object

taxonomy <- function(phyloseq_object) {
  tax_table(phyloseq_object) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(OTU = rownames(.)) %>%
    as_tibble() %>%
    arrange(OTU) %>%
    dplyr::select(OTU, Kingdom, Phylum, Class, Order, Family, Genus)
}

# Merge OTU count table with taxonomy table

otu_rel_abund <- function(data, otu, tax, samples, tax_rank) {
  inner_join(data, otu, by = "SampleID") %>%
    inner_join(., tax, by = "OTU") %>%
    dplyr::group_by(SampleID) %>%
    mutate(Rel_Abund = Count / sum(Count)) %>%
    ungroup() %>%
    dplyr::select(-Count) %>%
    pivot_longer(cols = c("Kingdom", 
                          "Phylum", 
                          "Class", 
                          "Order", 
                          "Family", 
                          "Genus", 
                          "OTU"),
                 names_to = "Level",
                 values_to = "Taxon") %>%
    mutate(SampleID = factor(SampleID, levels = samples)) %>%
    filter(Level == tax_rank) %>%
    dplyr::group_by(SampleID, Hours_since_death, Taxon) %>%
    rename(Death_hours = Hours_since_death) %>%
    summarise(Rel_Abund = sum(Rel_Abund), .groups = "drop") %>%
    dplyr::group_by(SampleID, Death_hours, Taxon) %>%
    summarise(Mean_Rel_Abund = 100*mean(Rel_Abund), .groups = "drop") %>%
    arrange(SampleID, desc(Mean_Rel_Abund))
}

# Create a taxa pool by filtering taxa with relative abundance above a specified cut-off value

taxon_pool <- function(otu_abund, cutoff) {
  dplyr::group_by(otu_abund, Taxon) %>%
    summarise(Pool = max(Mean_Rel_Abund) < cutoff, 
              Mean = mean(Mean_Rel_Abund),
              .groups = "drop")
}

relative_abundance <- function(rel_abun, tax_pool){
  inner_join(rel_abun, tax_pool, by = "Taxon") %>%
    mutate(Taxon = if_else(Pool, "Others", Taxon)) %>%
    dplyr::group_by(SampleID, Death_hours, Taxon) %>%
    summarise(Mean_Rel_Abund = sum(Mean_Rel_Abund), 
              Mean = sum(Mean),
              .groups = "drop") %>%
    mutate(Taxon = factor(Taxon),
           Taxon = fct_reorder(Taxon, Mean, .desc = TRUE),
           Taxon = fct_shift(Taxon, n = 1))
}

# Select specific variables of interest from the metadata table to create a new data table

data <- metadata %>% 
  dplyr::select(SampleID, Specie, Hours_since_death, Landscape) %>% 
  as_tibble()

# Relative abundance analysis at phylum level ----

# Amphisbaena bassleri ----

# Count unique taxa within each taxonomic level

cat("Unique taxa for A. bassleri:", "\n",
    "Kingdom: ", length(unique(tax_table(Amphisbaena_results)[, "Kingdom"])), "\n",
    "Phylum: ", length(unique(tax_table(Amphisbaena_results)[, "Phylum"])), "\n",
    "Class: ", length(unique(tax_table(Amphisbaena_results)[, "Class"])), "\n",
    "Order: ", length(unique(tax_table(Amphisbaena_results)[, "Order"])), "\n",
    "Family: ", length(unique(tax_table(Amphisbaena_results)[, "Family"])), "\n",
    "Genus: ", length(unique(tax_table(Amphisbaena_results)[, "Genus"])), "\n")

# Count the number of families within each phylum

count_taxlevel(Amphisbaena_results, 
               HigherTaxLevel = "Phylum",
               DeeperTaxLevel = "Family")

# Relative abundance at the phylum level

amph_otu_count <- otu_counts(Amphisbaena_results)

amph_taxonomy <- taxonomy(Amphisbaena_results)

amph_phylum_otu_rel_abund <- otu_rel_abund(data = data,
                                           otu = amph_otu_count,
                                           tax = amph_taxonomy,
                                           samples = c("SW001", 
                                                       "SW002", 
                                                       "SW003", 
                                                       "SW004"),
                                           tax_rank = "Phylum")

amph_phylum_tax_pool <- taxon_pool(otu_abund = amph_phylum_otu_rel_abund,
                                   cutoff = 3)

amph_phylum_rel_abund <- relative_abundance(rel_abun = amph_phylum_otu_rel_abund,
                                            tax_pool = amph_phylum_tax_pool)

# Bar plot of relative abundance at the phylum level

plot_amph_rel_abund_phylum <- 
  ggplot(amph_phylum_rel_abund, 
         aes(x = SampleID, 
             y = Mean_Rel_Abund, 
             fill = Taxon)) + 
  geom_col(position = "fill") +
  labs(fill = "Phylum") +
  xlab("Sample ID") +
  ylab("Relative Abundance") + 
  theme_classic2() +
  rotate_x_text(90) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(breaks = c("Firmicutes",
                               "Verrucomicrobia",
                               "Bacteroidetes",
                               "Proteobacteria",
                               "Cyanobacteria",
                               "Others"),
                    values = c("#7570B3", 
                               "#66A61E", 
                               "#1B9E77", 
                               "#E7298A", 
                               "#D95F02", 
                               "#B5B5B5"))

plot_amph_rel_abund_phylum

# Crotophaga ani ----

# Count unique taxa within each taxonomic level

cat("Unique taxa for C. ani:", "\n",
    "Kingdom: ", length(unique(tax_table(Crotophaga_results)[, "Kingdom"])), "\n",
    "Phylum: ", length(unique(tax_table(Crotophaga_results)[, "Phylum"])), "\n",
    "Class: ", length(unique(tax_table(Crotophaga_results)[, "Class"])), "\n",
    "Order: ", length(unique(tax_table(Crotophaga_results)[, "Order"])), "\n",
    "Family: ", length(unique(tax_table(Crotophaga_results)[, "Family"])), "\n",
    "Genus: ", length(unique(tax_table(Crotophaga_results)[, "Genus"])), "\n")

# Count the number of families within each phylum

count_taxlevel(Crotophaga_results, 
               HigherTaxLevel = "Phylum",
               DeeperTaxLevel = "Family") %>% print(., n = 21)

# Relative abundance at the phylum level

crot_otu_count <- otu_counts(Crotophaga_results)

crot_taxonomy <- taxonomy(Crotophaga_results)

crot_phylum_otu_rel_abund <- otu_rel_abund(data = data,
                                           otu = crot_otu_count,
                                           tax = crot_taxonomy,
                                           samples = c("SW005", 
                                                       "SW006", 
                                                       "SW007", 
                                                       "SW008",
                                                       "SW009"),
                                           tax_rank = "Phylum")

crot_phylum_tax_pool <- taxon_pool(otu_abund = crot_phylum_otu_rel_abund,
                                   cutoff = 3)

crot_phylum_rel_abund <- relative_abundance(rel_abun = crot_phylum_otu_rel_abund,
                                     tax_pool = crot_phylum_tax_pool)

# Bar plot of relative abundance at the phylum level

plot_crot_rel_abund_phylum <- 
  ggplot(crot_phylum_rel_abund, 
         aes(x = SampleID, 
             y = Mean_Rel_Abund, 
             fill = Taxon)) + 
  geom_col(position = "fill") +
  labs(fill = "Phylum") +
  xlab("Sample ID") +
  ylab("Relative Abundance") + 
  theme_classic2() +
  rotate_x_text(90) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(breaks = c("Firmicutes",
                               "Actinobacteria",
                               "Proteobacteria",
                               "Epsilonbacteraeota",
                               "Chlamydiae",
                               "Bacteroidetes",
                               "Others"),
                    values = c("#7570B3", 
                               "#E6AB02", 
                               "#E7298A", 
                               "#FDC086", 
                               "#B3CDE3", 
                               "#1B9E77", 
                               "#B5B5B5"))

plot_crot_rel_abund_phylum

# Alpha diversity ----

# Amphisbaena bassleri ----

# Alpha Diversity Indexes

amph_alpha_div <- estimate_richness(Amphisbaena_results,  
                                    measures = c("Observed",
                                                 "Shannon",
                                                 "Simpson",
                                                 "Chao1"))

rownames(amph_alpha_div) <- sample_data(Amphisbaena_results)$SampleID

amph_alpha_div

# Rarefaction curves

sample_data(Amphisbaena_results)$Hours_since_death <- 
  as(sample_data(Amphisbaena_results)$Hours_since_death, "character")

rarefaction_amph <- ranacapa::ggrare(Amphisbaena_results, 
                                     step = 600,
                                     se = F, 
                                     plot = F,
                                     color = "Hours_since_death")

rarefaction_amph$data$Sample <- rarefaction_amph$data$SampleID 

rarefac_amph_lab <- rarefaction_amph$data %>% 
  group_by(Sample) %>% 
  summarise(Size = max(Size), 
            `.S` = max(`.S`), 
            Hours_since_death = max(Hours_since_death))

plot_amph_rarefaction <- 
  rarefaction_amph + 
  geom_line(linewidth = 0.6) +
  ggrepel::geom_text_repel(data = rarefac_amph_lab, 
                           aes(x = Size, y =`.S`, label = Sample),
                           size = 3,
                           nudge_x = 0.1,
                           nudge_y = 0.01,
                           box.padding = 0.2,
                           segment.color = NA) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  labs(colour = "Time since death") +
  theme_light() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align = 0.5,
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +
  scale_colour_manual(breaks = c("0", 
                                 "2", 
                                 "6"), 
                      labels = c("0 hours", 
                                 "2 hours", 
                                 "6 hours"),
                      values = c("#00008B", 
                                 "#FF3030", 
                                 "#7FFF00"))

plot_amph_rarefaction

# Crotophaga ani ----

# Alpha Diversity Indexes

crot_alpha_div <- estimate_richness(Crotophaga_results,  
                                    measures = c("Observed",
                                                 "Shannon",
                                                 "Simpson",
                                                 "Chao1"))

rownames(crot_alpha_div) <- sample_data(Crotophaga_results)$SampleID

crot_alpha_div

# Rarefaction curves

sample_data(Crotophaga_results)$Hours_since_death <- 
  as(sample_data(Crotophaga_results)$Hours_since_death, "character")

rarefaction_crot <- ranacapa::ggrare(Crotophaga_results, 
                                     step = 600,
                                     se = F, 
                                     plot = F,
                                     color = "Hours_since_death")

rarefaction_crot$data$Sample <- rarefaction_crot$data$SampleID 

rarefac_crot_lab <- rarefaction_crot$data %>% 
  group_by(Sample) %>% 
  summarise(Size = max(Size), 
            `.S`=max(`.S`), 
            Hours_since_death = max(Hours_since_death))

plot_crot_rarefaction <- 
  rarefaction_crot + 
  geom_line(linewidth = 0.6) +
  ggrepel::geom_text_repel(data = rarefac_crot_lab, 
                           aes(x = Size, y = `.S`, label = Sample),
                           size = 3,  
                           nudge_x = 0.1,  
                           nudge_y = 0.01,  
                           box.padding = 0.2,
                           segment.color = NA) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  labs(colour = "Time since death") +
  theme_light() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align = 0.5,
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +
  scale_colour_manual(breaks = c("1", 
                                 "2", 
                                 "6", 
                                 "48"), 
                      labels = c("1 hours", 
                                 "2 hours", 
                                 "6 hours", 
                                 "48 hours"),
                      values = c("#FF66CC", 
                                 "#FF3030", 
                                 "#7FFF00", 
                                 "#FF9900"))

plot_crot_rarefaction

# Beta diversity ----

# Amphisbaena bassleri ----

# Principal Coordinate Analysis (PCoA)

# Normalize sequencing depth across samples

set.seed(231206)

amph_rarefied <- rarefy_even_depth(Amphisbaena_results, 
                                   sample.size = min(sample_sums(Amphisbaena_results)),
                                   rngseed = TRUE)

# Calculate PCoA with Bray-Curtis distance

amph_pcoa <- ordinate(amph_rarefied, method = "PCoA", distance = "bray")

# Plotting PCoA

plot_amph_pcoa <- 
  plot_ordination(physeq = amph_rarefied, 
                  ordination = amph_pcoa, 
                  color = "Hours_since_death",
                  shape = "SampleID") +
  geom_point(size = 3) +
  labs(color = "Hours since death",
       shape = "Sample ID") +
  theme_light() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align = 0.5,
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +        
  scale_shape_manual(values = c("SW001" = 15,
                                "SW002" = 16,
                                "SW003" = 17,
                                "SW004" = 18)) + 
  scale_colour_manual(breaks = c("0", 
                                 "2", 
                                 "6"), 
                      labels = c("0 hours", 
                                 "2 hours", 
                                 "6 hours"),
                      values = c("#00008B", 
                                 "#FF3030", 
                                 "#7FFF00")) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2))

plot_amph_pcoa

# Crotophaga ani ----

# Normalize sequencing depth across samples

set.seed(231206)

crot_rarefied <- rarefy_even_depth(Crotophaga_results, 
                                   sample.size = min(sample_sums(Crotophaga_results)),
                                   rngseed = TRUE)

# Calculate PCoA with Bray-Curtis distance

crot_pcoa <- ordinate(crot_rarefied, method = "PCoA", distance = "bray")

# Plotting PCoA

plot_crot_pcoa <- 
  plot_ordination(physeq = crot_rarefied, 
                  ordination = crot_pcoa,
                  color = "Hours_since_death",
                  shape = "SampleID") +
  geom_point(size = 3) +
  labs(color = "Hours since death",
       shape = "Sample ID") +
  theme_light() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align = 0.5,
        axis.text.x = element_text(face = "bold", size = 10),
        axis.text.y = element_text(face = "bold", size = 10)) +
  scale_shape_manual(values = c("SW005" = 15,
                                "SW006" = 16,
                                "SW007" = 17,
                                "SW008" = 18,
                                "SW009" = 8)) +
  scale_colour_manual(breaks = c("1", 
                                 "2", 
                                 "6", 
                                 "48"), 
                      labels = c("1 hours", 
                                 "2 hours", 
                                 "6 hours", 
                                 "48 hours"),
                      values = c("#FF66CC", 
                                 "#FF3030", 
                                 "#7FFF00", 
                                 "#FF9900")) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2))

plot_crot_pcoa

# Function to prepare data for constructing a community matrix ----

preprocess_data <- function(physeq_object, cols_samples_names) {
  
  # Extract raw OTU table
  OTU_raw <- as.data.frame(otu_table(physeq_object))
  
  # Extract raw taxonomy table
  tax_raw <- as.data.frame(tax_table(physeq_object))
  
  # Join both raw taxonomy with OTU table
  comm_raw <- cbind(tax_raw, OTU_raw)
  
  # Pivot and normalize data
  data_gt <- pivot_longer(comm_raw, 
                          cols = all_of(cols_samples_names),
                          values_to = "num_hits", 
                          names_to = "Sample") %>% 
    group_by(Sample) %>% 
    mutate(X._hits = num_hits/sum(num_hits))
  
  return(data_gt)
}

# Function to build a community matrix from OTU and taxonomy data ----

communi_matrix <- function(dat, level_tax, field, tlrn = 1e-5, rt = 1) {
  
  community_dt1 <-  dat %>%
    rename(Otu_level = all_of(level_tax)) %>%
    group_by(Sample, Otu_level)
  
  switch(field, 
         Per = {community_dt1 <- community_dt1 %>% 
           summarise(Aggregated = sum(X._hits, na.rm = T), .groups = "drop")},
         Count = {community_dt1 <- community_dt1 %>%
           summarise(Aggregated = sum(num_hits, na.rm = T), .groups = "drop")})
  
  community_dt <- community_dt1 %>%
    group_by(Sample, Otu_level) %>%
    summarise(Aggregated = sum(Aggregated, na.rm = T), .groups = "drop") %>% 
    filter(Aggregated > tlrn) %>%
    spread(all_of(Otu_level), Aggregated) %>%
    mutate_at(vars(-group_cols()), list(~ifelse(is.na(.), 0, .))) %>%
    as.data.frame()
  
  rownames(community_dt) <- as.character(community_dt$Sample)
  
  community_mt <- community_dt
  
  community_mt <- community_mt[,-1]
  
  ifelse(rt == 1, 
         return(community_dt),
         return(as.matrix(community_mt[, -c(length(community_dt), 
                                            length(community_dt)-1)])))
}

# Function to convert a heatmap object into a plot ----

heatmap_to_ggplot <- function(heatmap_object) {
  ggplot_new <- grid.grabExpr(draw(heatmap_object)) %>%
    ggpubr::as_ggplot()
}

# Heatmaps ----

# Amphisbaena bassleri ----

# Prepare data for constructing a community matrix

data_gt_amph <- preprocess_data(physeq_object = Amphisbaena_results,
                                cols_samples_names = c("SW001", 
                                                       "SW002", 
                                                       "SW003", 
                                                       "SW004"))

# Build a community matrix from OTU and taxonomy data

com_amph_hm <- communi_matrix(dat = data_gt_amph,
                              level_tax = "Family",
                              field = "Count",
                              tlrn = 1e-5,
                              rt = 2)

# Transform the abundances of families to a logarithmic scale

log_sums_amph <- as.data.frame(colSums(log(com_amph_hm + 1))) %>% 
  rename('log(x+1)' = 'colSums(log(com_amph_hm + 1))') %>%
  arrange(desc(`log(x+1)`))

# Obtaining the abundances of the 25 most dominant families

log_sums_amph[25, ]

# Reduce the families to the 25 most dominant

reduced_hm_amph <- log(com_amph_hm + 1)[, colSums(log(com_amph_hm + 1)) > 10.3] %>% 
  t()

# Rename some unclassified and uncultured taxa 

rownames(reduced_hm_amph) <- str_replace_all(rownames(reduced_hm_amph),
                                             c("(.*)_unclassified" = "Unclassified \\1",
                                               "uncultured" = "Uncultured Bacteria",
                                               "_" = " "))

# Heatmap annotations

# Extract data for postmortem hours

postmortem_hours_amph <- sample_data(Amphisbaena_results)$Hours_since_death

# Set up colors and legends for postmortem hours dataset

color_postmortem_amph <- c("0" = "#F6A97A", 
                           "2" = "#D44292", 
                           "6" = "#952EA0") 

legend_postmortem_amph <- c("0 hours", 
                            "2 hours", 
                            "6 hours")

# Extract data for landscape

landscape_amph <- sample_data(Amphisbaena_results)$Landscape

# Set up colors and legends for landscape dataset

color_landscape_amph <- c("Altered area" = "#FF0000",
                          "Unaltered area" = "#7FFF00")

# Create annotations for heatmap

hm_annot_amph <- 
  HeatmapAnnotation("Time since death" = postmortem_hours_amph,
                    "Landscape" = landscape_amph,
                    col = list("Time since death" = color_postmortem_amph,
                               "Landscape" = color_landscape_amph),
                    annotation_legend_param = list("Time since death" = list(title = "Time since death",
                                                                             at = names(color_postmortem_amph),
                                                                             labels = legend_postmortem_amph),
                                                   "Landscape" = list(title = "Landscape",
                                                                             at = names(color_landscape_amph))),
                    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_legend = TRUE,
                    annotation_name_side = "right")

# Create the heatmap object

amph_heatmap <- 
  Heatmap(reduced_hm_amph,
          name = 'log OTUs Abundance',
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(legend_direction = "horizontal", 
                                      legend_width = unit(4, "cm")),
          top_annotation = hm_annot_amph,
          row_names_gp = grid::gpar(fontsize = 10),
          column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
          column_dend_height = unit(1, "cm"), 
          row_dend_width = unit(1.5, "cm"),
          col = mako(n = 20),
          cluster_columns = TRUE,
          use_raster = FALSE, 
          row_names_side = "right", 
          row_dend_side = "left",
          width = unit(35, "mm"))

# Draw the heatmap

plot_amph_heatmap <- draw(amph_heatmap, 
                          heatmap_legend_side = "bottom")

# Convert the heatmap object into a plot

plot_amph_heatmap <- heatmap_to_ggplot(plot_amph_heatmap)

plot_amph_heatmap

# Crotophaga ani ----

# Prepare data for constructing a community matrix

data_gt_crot <- preprocess_data(physeq_object = Crotophaga_results,
                                cols_samples_names = c("SW005", 
                                                       "SW006", 
                                                       "SW007", 
                                                       "SW008", 
                                                       "SW009"))
# Build a community matrix from OTU and taxonomy data

com_crot_hm <- communi_matrix(dat = data_gt_crot,
                              level_tax = "Family",
                              field = "Count",
                              tlrn = 1e-5,
                              rt = 2)

# Transform the abundances of families to a logarithmic scale

log_sums_crot <- as.data.frame(colSums(log(com_crot_hm + 1))) %>% 
  rename('log(x+1)' = 'colSums(log(com_crot_hm + 1))') %>%
  arrange(desc(`log(x+1)`))

# Obtaining the abundances of the 25 most dominant families

log_sums_crot[25, ]

# Reduce the families to the 25 most dominant

reduced_hm_crot <- log(com_crot_hm + 1)[, colSums(log(com_crot_hm + 1)) > 17.9] %>% 
  t()

# Rename some unclassified and uncultured taxa 

rownames(reduced_hm_crot) <- str_replace_all(rownames(reduced_hm_crot),
                                             c("(.*)_unclassified" = "Unclassified \\1",
                                               "uncultured" = "Uncultured Bacteria",
                                               "_" = " "))

# Heatmap annotations

# Extract data for postmortem hours

postmortem_hours_crot <- sample_data(Crotophaga_results)$Hours_since_death

# Set up colors and legends for postmortem hours dataset

color_postmortem_crot <- c("1" = "#F66D7A", 
                           "2" = "#D44292", 
                           "6" = "#952EA0",
                           "48" = "#4B2991") 

legend_postmortem_crot <- c("1 hours", 
                            "2 hours", 
                            "6 hours",
                            "48 hours")

# Extract data for landscape

landscape_crot <- sample_data(Crotophaga_results)$Landscape

# Set up colors and legends for landscape dataset

color_landscape_crot <- c("Altered area" = "#FF0000",
                          "Unaltered area" = "#7FFF00")

# Create annotations for heatmap

hm_annot_crot <- 
  HeatmapAnnotation("Time since death" = postmortem_hours_crot,
                    "Landscape" = landscape_crot,
                    col = list("Time since death" = color_postmortem_crot,
                               "Landscape" = color_landscape_crot),
                    annotation_legend_param = list("Time since death" = list(title = "Time since death",
                                                                             at = names(color_postmortem_crot),
                                                                             labels = legend_postmortem_crot),
                                                   "Landscape" = list(title = "Landscape",
                                                                      at = names(color_landscape_crot))),
                    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_legend = TRUE,
                    annotation_name_side = "right")

# Create the heatmap object

crot_heatmap <- 
  Heatmap(reduced_hm_crot,
          name = 'log OTUs Abundance',
          show_heatmap_legend = TRUE,
          heatmap_legend_param = list(legend_direction = "horizontal", 
                                      legend_width = unit(4, "cm")),
          top_annotation = hm_annot_crot,
          row_names_gp = grid::gpar(fontsize = 10),
          column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
          column_dend_height = unit(1, "cm"), 
          row_dend_width = unit(1.5, "cm"),
          col = mako(n = 20),
          cluster_columns = TRUE,
          use_raster = FALSE, 
          row_names_side = "right", 
          row_dend_side = "left",
          width = unit(35, "mm"))

# Draw the heatmap

plot_crot_heatmap <- draw(crot_heatmap, 
                          heatmap_legend_side = "bottom")

# Convert the heatmap object into a plot

plot_crot_heatmap <- heatmap_to_ggplot(plot_crot_heatmap)

plot_crot_heatmap

# Function to filter OTUs based on prevalence and abundance ----

venn_core_microbiota <- function (physeq_object, group_name, abund_value, prev_value) {
  
  comp <- microbiome::transform(physeq_object, "compositional")
  
  comp.n <- format_to_besthit(comp)
  
  taxa_names(comp.n) <- str_replace_all(taxa_names(comp.n),
                                        c("uncultured" = "Uncultured Bacteria",
                                          "_" = " "))
  
  var_names <- unique(as.character(meta(comp.n)[[group_name]]))
  
  list_core <- list()
  
  for (i in var_names) {
    filter_expr <- paste0("subset_samples(comp.n, ", group_name, " == '", i, "')")
    sub <- eval(parse(text = filter_expr))
    core_microbiome <- core_members(sub, 
                                    detection = abund_value, 
                                    prevalence = prev_value)
    list_core[[i]] <- core_microbiome
  }
  
  return(list_core)
}

# Core Gut Microbiota ----

# Amphisbaena bassleri ----

# Identifying the taxa that compose the core gut microbiota

amph_list_core <- venn_core_microbiota(physeq_object = Amphisbaena_results,
                                       group_name = "SampleID",
                                       abund_value = 0.01/100,
                                       prev_value = 90/100)

# Venn diagram of the core gut microbiota

plot_amph_core_microbiota <-
  venn.diagram(x = amph_list_core,
               filename = NULL,
               disable.logging = TRUE,
               fill = c("#00008B", 
                        "#FF0000", 
                        "#FFFF00", 
                        "#00FF00"),
               alpha = 0.5,
               cat.fontface = 2,
               cex = 1,
               cat.cex = 1,
               print.mode = c("raw"),
               fontface = 2)

plot_amph_core_microbiota <- ggpubr::as_ggplot(plot_amph_core_microbiota)

plot_amph_core_microbiota

# Visualization of taxa from the core gut microbiota in a data frame

amph_core_microbiota <- Reduce(intersect, amph_list_core)

amph_core_microbiota <- sub(":.*", "", amph_core_microbiota)

amph_core_microbiota <- amph_taxonomy[amph_taxonomy$OTU %in% amph_core_microbiota, ]

amph_core_microbiota

# Crotophaga ani ----

# Identifying the taxa that compose the core gut microbiota

crot_list_core <- venn_core_microbiota(physeq_object = Crotophaga_results,
                                       group_name = "SampleID",
                                       abund_value = 0.01/100,
                                       prev_value = 90/100)

# Venn diagram for the core gut microbiota

plot_crot_core_microbiota <- 
  venn.diagram(x = crot_list_core,
               filename = NULL,
               disable.logging = TRUE,
               fill = c("#FFA500", 
                        "#FFEC8B", 
                        "#9AC0CD", 
                        "#C0FF3E", 
                        "#FF6EB4"),
               alpha = 0.5,
               cat.fontface = 2,
               cat.just = list(c(NA, 1), 
                               c(-0.5,-4), 
                               c(1,-0.5), 
                               c(NA, NA), 
                               c(1, -3)),
               cex = 1,
               cat.cex = 1,
               print.mode = c("raw"),
               fontface = 2)

plot_crot_core_microbiota <- ggpubr::as_ggplot(plot_crot_core_microbiota)

plot_crot_core_microbiota

# Visualization of taxa from the core gut microbiota in a data frame

crot_core_microbiota <- Reduce(intersect, crot_list_core)

crot_core_microbiota <- sub(":.*", "", crot_core_microbiota)

crot_core_microbiota <- crot_taxonomy[crot_taxonomy$OTU %in% crot_core_microbiota, ]

crot_core_microbiota

# Merge plots ----

# Microbiota composition of A. bassleri

FigMicrobiotaComposition_Amph <- grid.arrange(plot_amph_rel_abund_phylum, plot_amph_heatmap,
                                              widths = c(1, 0.01, 1.5),
                                              layout_matrix = rbind(c(1, NA, 2)))

Fig1 <- ggpubr::as_ggplot(FigMicrobiotaComposition_Amph) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.4),
                  y = c(1, 1)) +  
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

Fig1

# Microbiota composition of C. ani

FigMicrobiotaComposition_Crot <- grid.arrange(plot_crot_rel_abund_phylum, plot_crot_heatmap,
                                              widths = c(1, 0.01, 1.5),
                                              layout_matrix = rbind(c(1, NA, 2)))

Fig2 <- ggpubr::as_ggplot(FigMicrobiotaComposition_Crot) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.4),
                  y = c(1, 1)) +  
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))


Fig2

# PCoA

FigPCoA <- grid.arrange(plot_amph_pcoa, plot_crot_pcoa,
                           widths = c(1, 0.05, 1),
                           layout_matrix = rbind(c(1, NA, 2)))

Fig3 <- ggpubr::as_ggplot(FigPCoA) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.5),
                  y = c(1, 1)) +  
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))


Fig3

# Core Gut Microbiota

FigCoreGut <- grid.arrange(plot_amph_core_microbiota, plot_crot_core_microbiota,
                           widths = c(1, 0.05, 1),
                           layout_matrix = rbind(c(1, NA, 2)))

Fig4 <- ggpubr::as_ggplot(FigCoreGut) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.5),
                  y = c(1, 1)) +  
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))


Fig4

# Save plots ----

# Microbiota composition of Amphisbaena bassleri (Fig2)

ggsave("./Results/Microbiota/Plots/Fig1_MicrobiotaCompositionAmphisbaena.png", 
       plot = Fig1, 
       width = 5.85, 
       height = 3.75, 
       dpi = 1000, 
       scale = 1.6)

# Microbiota composition of Crotophaga ani (Fig3)

ggsave("./Results/Microbiota/Plots/Fig2_MicrobiotaCompositionCrotophaga.png", 
       plot = Fig2, 
       width = 5.85, 
       height = 3.75, 
       dpi = 1000, 
       scale = 1.6)

# Beta diversity using PCoA (Fig4)

ggsave("./Results/Microbiota/Plots/Fig3_PCoA.png", 
       plot = Fig3, 
       width = 6, 
       height = 2.5, 
       dpi = 1000, 
       scale = 1.5)

# Core gut microbiota using venn diagrams for A. bassleri and C. ani

ggsave("./Results/Microbiota/Plots/Fig4_CoreGutMicrobiota.png", 
       plot = Fig4, 
       width = 3, 
       height = 1.5, 
       dpi = 1000, 
       scale = 3)