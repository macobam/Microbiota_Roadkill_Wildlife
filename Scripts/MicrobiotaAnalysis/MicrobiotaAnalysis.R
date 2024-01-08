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

# Note: Please consider the full paths to the input and output files.

# Data importing from mothur results ----

# Outputs files names

sharedfile <- "MothurResults/OtuTable.shared"

taxfile <- "MothurResults/TaxTable.taxonomy"

# Creating a phyloseq object

mothur_raw_results = import_mothur(mothur_shared_file = sharedfile,
                                   mothur_constaxonomy_file = taxfile)

# Adding taxonomy information to phyloseq object

colnames(tax_table(mothur_raw_results)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# Importing metadata

metadata <- read_xlsx("MetadataRoadkillWildlife.xlsx", sheet = "Metadata") %>% 
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

mothur_results <- subset_taxa(mothur_results, Kingdom != "unknown" & Kingdom != "Archaea")

mothur_results

# Remove unclassified bacteria

table(tax_table(mothur_results)[,2])

mothur_results <- subset_taxa(mothur_results, Phylum != "Bacteria_unclassified")

mothur_results

# Remove singletons

sum(taxa_sums(mothur_results) == 1)

mothur_results <- prune_taxa(taxa_sums(mothur_results) > 1, mothur_results)

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
             nrow=1) +
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
        legend.title.align=0.5,
        legend.text = element_text(face = "italic",
                                   size = 10),
        axis.text.x = element_text(face = "bold", 
                                   size = 10),
        axis.text.y = element_text(face = "bold", 
                                   size = 10))

plot_seqdepth

# Stats sequencing depth by specie

meta(mothur_results) %>% 
  group_by(Specie) %>% 
  summarise(mean = mean(SeqDepth), 
            sd = sd(SeqDepth), 
            median = median(SeqDepth), 
            Total = sum(SeqDepth),
            .groups = "drop")

# Separating samples per species ----

Species = phyloseq_sep_variable(mothur_results, "Specie")

Amphisbaena_results = Species$`Amphisbaena bassleri`

Crotophaga_results = Species$`Crotophaga ani`

# Function to count the number of families by Phylum ----

count_families <- function (phyloseq_obj) {
  tax_table <- tax_table(phyloseq_obj)
  
  tax_df <- as.data.frame(tax_table) %>%
    rownames_to_column(var = "OTU_ID") %>%
    filter(!is.na(Phylum))
  
  count_per_phylum <- tax_df %>%
    group_by(Phylum) %>%
    summarize(num_families = n_distinct(Family)) %>% 
    arrange(desc(num_families))
  
  print(count_per_phylum)
}  

# Functions for Relative Abundance Analysis ----

# Create Otu count table

otu_counts <- function(phyloseq_object){
  otu_table(phyloseq_object) %>%
    psmelt() %>%
    as_tibble() %>%
    arrange(Sample, OTU) %>%
    rename(SampleID=Sample,Count=Abundance) %>%
    dplyr::select(SampleID, OTU, Count)
}

# Create taxonomy table

taxonomy <- function(phyloseq_object) {
  tax_table(phyloseq_object) %>% 
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(OTU = rownames(.)) %>%
    as_tibble() %>%
    arrange(OTU) %>%
    dplyr::select(OTU, Kingdom, Phylum, Class, Order, Family, Genus)
}

# Join Otu count table and Taxonomy table

otu_rel_abund <- function(data, otu, tax, samples, tax_rank) {
  inner_join(data, otu, by = "SampleID") %>%
    inner_join(., tax, by = "OTU") %>%
    dplyr::group_by(SampleID) %>%
    mutate(Rel_Abund = Count / sum(Count)) %>%
    ungroup() %>%
    dplyr::select(-Count) %>%
    pivot_longer(cols=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "OTU"),
                 names_to="Level",
                 values_to = "Taxon") %>%
    mutate(SampleID = factor(SampleID, levels = samples)) %>%
    filter(Level==tax_rank) %>%
    dplyr::group_by(SampleID, Hours_since_death, Taxon) %>%
    rename(Death_hours = Hours_since_death) %>%
    summarise(Rel_Abund = sum(Rel_Abund), .groups = "drop") %>%
    dplyr::group_by(SampleID, Death_hours, Taxon) %>%
    summarise(Mean_Rel_Abund = 100*mean(Rel_Abund), .groups = "drop") %>%
    arrange(SampleID, desc(Mean_Rel_Abund))
}

# Create a taxon pool with a cutoff value for relative abundance

taxon_pool <-function(otu_abund, cutoff) {
  dplyr::group_by(otu_abund, Taxon) %>%
    summarise(Pool=max(Mean_Rel_Abund) < cutoff, 
              Mean=mean(Mean_Rel_Abund),
              .groups="drop")
}

relative_abundance <- function(rel_abun, tax_pool){
  inner_join(rel_abun, tax_pool, by="Taxon") %>%
    mutate(Taxon = if_else(Pool, "Others", Taxon)) %>%
    dplyr::group_by(SampleID, Death_hours, Taxon) %>%
    summarise(Mean_Rel_Abund = sum(Mean_Rel_Abund), 
              Mean = sum(Mean),
              .groups = "drop") %>%
    mutate(Taxon = factor(Taxon),
           Taxon = fct_reorder(Taxon, Mean, .desc=TRUE),
           Taxon = fct_shift(Taxon, n = 1))
}

# Selection of specific interest variables from metadata table to create a new data table

data <- metadata %>% 
  dplyr::select(SampleID, Specie, Hours_since_death, Landscape) %>% 
  as_tibble()

# Relative abundance analysis at Phylum level ----

# Amphisbaena bassleri ----

# Count unique taxons at Phylum level

length(unique(tax_table(Amphisbaena_results)[, "Phylum"]))

# Number of families by Phylum

count_families(Amphisbaena_results)

# Relative abundance at Phylum level

amph_otu_count <- otu_counts(Amphisbaena_results)

amph_taxonomy <- taxonomy(Amphisbaena_results)

amph_otu_rel_abund <- otu_rel_abund(data = data,
                                    otu = amph_otu_count,
                                    tax = amph_taxonomy,
                                    samples = c("SW001", 
                                                "SW002", 
                                                "SW003", 
                                                "SW004"),
                                    tax_rank = "Phylum")

amph_tax_pool <- taxon_pool(otu_abund = amph_otu_rel_abund,
                            cutoff = 3)

amph_rel_abund <- relative_abundance(rel_abun = amph_otu_rel_abund,
                                     tax_pool = amph_tax_pool)

# Bar plot of relative abundance

plot_amph_rel_abund_phylum <- 
  ggplot(amph_rel_abund, 
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
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold"),
        legend.title.align=0.5,
        legend.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold")) +
  scale_fill_manual(breaks = c("Bacteroidetes",
                               "Cyanobacteria",
                               "Firmicutes", 
                               "Proteobacteria",
                               "Verrucomicrobia",
                               "Others"),
                    values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#B5B5B5"))

plot_amph_rel_abund_phylum

# Crotophaga ani ----

# Count unique taxons at Phylum level

length(unique(tax_table(Crotophaga_results)[, "Phylum"]))

# Number of families by Phylum

count_families(Crotophaga_results) %>% print(., n=21)

# Relative abundance at Phylum level

crot_otu_count <- otu_counts(Crotophaga_results)

crot_taxonomy <- taxonomy(Crotophaga_results)

crot_otu_rel_abund <- otu_rel_abund(data = data,
                                    otu = crot_otu_count,
                                    tax = crot_taxonomy,
                                    samples = c("SW005", 
                                                "SW006", 
                                                "SW007", 
                                                "SW008",
                                                "SW009"),
                                    tax_rank = "Phylum")

crot_tax_pool <- taxon_pool(otu_abund = crot_otu_rel_abund,
                            cutoff = 3)

crot_rel_abund <- relative_abundance(rel_abun = crot_otu_rel_abund,
                                     tax_pool = crot_tax_pool)

# Bar plot of relative abundance

plot_crot_rel_abund_phylum <- 
  ggplot(crot_rel_abund, 
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
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold"),
        legend.title.align=0.5,
        legend.spacing.y = unit(0.3, "cm"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold")) +
  scale_fill_manual(breaks = c("Actinobacteria",
                               "Bacteroidetes",
                               "Chlamydiae",
                               "Epsilonbacteraeota",
                               "Firmicutes",
                               "Proteobacteria",
                               "Others"),
                    values = c("#E6AB02", "#1B9E77", "#B3CDE3", "#FDC086", "#7570B3", "#E7298A", "#B5B5B5"))

plot_crot_rel_abund_phylum

# Alpha diversity ----

# Amphisbaena bassleri ----

# Measures indexes

amph_alpha_div<-estimate_richness(Amphisbaena_results,  
                                  measures = c("Chao1", 
                                               "Simpson", 
                                               "Shannon"))

rownames(amph_alpha_div)<-sample_data(Amphisbaena_results)$SampleID

amph_alpha_div

# Kruskal-Wallis test in Simpson and Shannon index

# Simpson index
kruskal.test(amph_alpha_div$Simpson, g=factor(rownames(amph_alpha_div)))

# Shannon index
kruskal.test(amph_alpha_div$Shannon, g=factor(rownames(amph_alpha_div)))

# Linear regression with Chao1 index as a function of the time since death

sample_data(Amphisbaena_results)$Hours_since_death <- 
  as(sample_data(Amphisbaena_results)$Hours_since_death, "double")

tmp_amph <- plot_richness(Amphisbaena_results)

amph_lm_data <- tmp_amph$data %>% filter(variable %in% "Chao1")

amph_rlm <- rlm(value ~ Hours_since_death, data = amph_lm_data, maxit = 25)

f.robftest(amph_rlm, var="Hours_since_death")

plot_amph_chao1 <- 
  ggplot(amph_lm_data, 
         aes(x = Hours_since_death, 
             y = value)) +
  stat_smooth(formula = y ~ x, 
              method = "rlm", 
              col = "#8B0A50", 
              linewidth = 0.5, 
              method.args = list(maxit = 25)) +
  geom_point(aes(colour = as.factor(Hours_since_death)), size = 3) + 
  xlab("Time since death (hours)") +
  ylab("Alpha diversity \n(Chao1)\n") +
  labs(colour = "Time since death") +  
  annotate("text", label = "p-value = 0.64", x = 4, y = 210, size = 3) +
  theme_classic() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +
  scale_colour_manual(breaks = c("0", "2", "6"), 
                      labels = c("0 hours", "2 hours", "6 hours"),
                      values = c("#00008B", "#FF3030", "#7FFF00"))

plot_amph_chao1

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
  summarise(Size=max(Size), `.S`=max(`.S`), Hours_since_death=max(Hours_since_death))

plot_amph_rarefaction <- 
  rarefaction_amph + 
  geom_line(linewidth = 0.6) +
  ggrepel::geom_text_repel(data=rarefac_amph_lab, aes(x=Size, y=`.S`, label=Sample),
                           size = 3,
                           nudge_x = 0.1,
                           nudge_y = 0.01,
                           box.padding = 0.2,
                           segment.color = NA) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  labs(colour = "Time since death") +
  theme_light() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +
  scale_colour_manual(breaks = c("0", "2", "6"), 
                      labels = c("0 hours", "2 hours", "6 hours"),
                      values = c("#00008B", "#FF3030", "#7FFF00"))

plot_amph_rarefaction

# Crotophaga ani ----

# Measures indexes

crot_alpha_div<-estimate_richness(Crotophaga_results,  
                                  measures = c("Chao1", 
                                               "Simpson", 
                                               "Shannon"))

rownames(crot_alpha_div)<-sample_data(Crotophaga_results)$SampleID

crot_alpha_div

# Kruskal-Wallis test in Simpson and Shannon index

# Simpson index
kruskal.test(crot_alpha_div$Simpson, g=factor(rownames(crot_alpha_div)))

# Shannon index
kruskal.test(crot_alpha_div$Shannon, g=factor(rownames(crot_alpha_div)))

# Linear regression with Chao1 index as a function of the time since death

sample_data(Crotophaga_results)$Hours_since_death <- 
  as(sample_data(Crotophaga_results)$Hours_since_death, "double")

tmp_crot <- plot_richness(Crotophaga_results)

crot_lm_data <- tmp_crot$data %>% filter(variable %in% "Chao1")

crot_rlm <- rlm(value ~ Hours_since_death, data = crot_lm_data, maxit = 25)

f.robftest(crot_rlm, var="Hours_since_death")

plot_crot_chao1 <- 
  ggplot(crot_lm_data, 
         aes(x = Hours_since_death, 
             y = value)) +
  stat_smooth(formula = y ~ x, method = "rlm", col = "#8B0A50", linewidth = 0.5) +
  geom_point(aes(colour = as.factor(Hours_since_death)), size = 3) + 
  xlab("Time since death (hours)") +
  ylab("Alpha diversity \n(Chao1)\n") +
  labs(colour = "Time since death") +  
  annotate("text", label = "p-value = 0.06", x = 35, y = 210, size = 3) +
  theme_classic() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +
  scale_colour_manual(breaks = c("1", "2", "6", "48"), 
                      labels = c("1 hours", "2 hours", "6 hours", "48 hours"),
                      values = c("#FF66CC", "#FF3030", "#7FFF00", "#FF9900"))

plot_crot_chao1

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
  summarise(Size=max(Size), `.S`=max(`.S`), Hours_since_death=max(Hours_since_death))

plot_crot_rarefaction <- 
  rarefaction_crot + 
  geom_line(linewidth = 0.6) +
  ggrepel::geom_text_repel(data=rarefac_crot_lab, aes(x=Size, y=`.S`, label=Sample),
                           size = 3,  
                           nudge_x = 0.1,  
                           nudge_y = 0.01,  
                           box.padding = 0.2,
                           segment.color = NA) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  labs(colour = "Time since death") +
  theme_light() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +
  scale_colour_manual(breaks = c("1", "2", "6", "48"), 
                      labels = c("1 hours", "2 hours", "6 hours", "48 hours"),
                      values = c("#FF66CC", "#FF3030", "#7FFF00", "#FF9900"))

plot_crot_rarefaction

# Beta diversity ----

# Amphisbaena bassleri ----

# Principal Coordinate Analysis (PCoA)

# Normalize the sequencing depth

set.seed(231206)

amph_rarefied <- rarefy_even_depth(Amphisbaena_results, 
                                   sample.size = min(sample_sums(Amphisbaena_results)),
                                   rngseed = TRUE)

# Calculate PCoA with Bray distance

amph_pcoa <- ordinate(amph_rarefied, method = "PCoA", distance = "bray")

# Ploting PCoA

plot_amph_pcoa <- 
  plot_ordination(physeq = amph_rarefied, 
                  ordination = amph_pcoa, 
                  color = "Hours_since_death",
                  shape = "SampleID") +
  geom_point(size = 3) +
  labs(color = "Hours since death",
       shape = "Sample ID") +
  theme_light() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +        
  scale_shape_manual(values = c("SW001" = 15,
                                "SW002" = 16,
                                "SW003" = 17,
                                "SW004" = 18)) + 
  scale_colour_manual(breaks = c("0", "2", "6"), 
                      labels = c("0 hours", "2 hours", "6 hours"),
                      values = c("#00008B", "#FF3030", "#7FFF00")) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2))

plot_amph_pcoa

#Crotophaga ani ----

# Normalize the sequencing depth

set.seed(231206)

crot_rarefied <- rarefy_even_depth(Crotophaga_results, 
                                   sample.size = min(sample_sums(Crotophaga_results)),
                                   rngseed = TRUE)

# Calculate PCoA with Bray distance

crot_pcoa <- ordinate(crot_rarefied, method = "PCoA", distance = "bray")

# Ploting PCoA

plot_crot_pcoa <- 
  plot_ordination(physeq = crot_rarefied, 
                  ordination = crot_pcoa,
                  color = "Hours_since_death",
                  shape = "SampleID") +
  geom_point(size = 3) +
  labs(color = "Hours since death",
       shape = "Sample ID") +
  theme_light() +
  theme(axis.title.y = element_text(face= "bold"),
        axis.title.x = element_text(face= "bold"),
        legend.title = element_text(face = "bold", size = 10),
        legend.title.align=0.5,
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)) +
  scale_shape_manual(values = c("SW005" = 15,
                                "SW006" = 16,
                                "SW007" = 17,
                                "SW008" = 18,
                                "SW009" = 8)) +
  scale_colour_manual(breaks = c("1", "2", "6", "48"), 
                      labels = c("1 hours", "2 hours", "6 hours", "48 hours"),
                      values = c("#FF66CC", "#FF3030", "#7FFF00", "#FF9900")) +
  guides(shape = guide_legend(order = 1),
         color = guide_legend(order = 2))

plot_crot_pcoa

# Function to prepare data previously to build community matrix ----

preprocess_data <- function(physeq_object, cols_samples_names) {
  
  # Extract raw otu table
  OTU_raw <- as.data.frame(otu_table(physeq_object))
  
  # Extract raw taxonomy table
  tax_raw <- as.data.frame(tax_table(physeq_object))
  
  # Join both raw taxonomy and otu table
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

# Function to build a community matrix ----

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

# Function to convert heatmap object into plot ----

heatmap_to_ggplot <- function(heatmap_object) {
  ggplot_new <- grid.grabExpr(draw(heatmap_object)) %>%
    ggpubr::as_ggplot()
}

# Heatmaps ----

# Amphisbaena bassleri ----

data_gt_amph <- preprocess_data(physeq_object = Amphisbaena_results,
                                cols_samples_names = c("SW001", 
                                                       "SW002", 
                                                       "SW003", 
                                                       "SW004"))

com_amph_hm <- communi_matrix(dat = data_gt_amph,
                              level_tax = "Family",
                              field = "Count",
                              tlrn = 1e-5,
                              rt = 2)

log_sums <- colSums(log(com_amph_hm + 1))

log_sums_amph <- as.data.frame(colSums(log(com_amph_hm + 1))) %>% 
  rename('log(x+1)' = 'colSums(log(com_amph_hm + 1))') %>%
  arrange(desc(`log(x+1)`))

log_sums_amph[25, ]

reduced_hm_amph <- log(com_amph_hm + 1)[, colSums(log(com_amph_hm + 1)) > 10.3] %>% 
  t()

rownames(reduced_hm_amph) <- str_replace_all(rownames(reduced_hm_amph),
                                             c("(.*)_unclassified" = "Unclassified \\1",
                                               "uncultured" = "Uncultured Bacteria",
                                               "_" = " "))

# Heatmap annotation

# Extract postmortem hours data

postmortem_hours_amp <- sample_data(Amphisbaena_results)$Hours_since_death

# Define colors and legends

color_postmortem_amph <- c("0" = "#F6A97A", 
                           "2" = "#D44292", 
                           "6" = "#952EA0") 

legend_postmortem_amph <- c("0 hours", 
                            "2 hours", 
                            "6 hours")

hm_annot_amph <- 
  HeatmapAnnotation("Time since death" = postmortem_hours_amp,
                    col = list("Time since death" = color_postmortem_amph),
                    annotation_legend_param = list("Time since death" = list(title = "Time since death",
                                                                             at = names(color_postmortem_amph),
                                                                             labels = legend_postmortem_amph)),
                    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_legend = FALSE,
                    annotation_name_side = "right")

# Create heatmap object

amph_heatmap <- 
  Heatmap(reduced_hm_amph,
          column_title = "A",
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
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
          cluster_columns = FALSE,
          use_raster = FALSE, 
          row_names_side = "right", 
          row_dend_side = "left",
          width = unit(35, "mm"))

# Draw heatmap

plot_amph_heatmap <- draw(amph_heatmap, 
                          heatmap_legend_side = "bottom")

# Convert heatmap to plot

plot_amph_heatmap <- heatmap_to_ggplot(plot_amph_heatmap)

plot_amph_heatmap

# Crotophaga ani ----

data_gt_crot <- preprocess_data(physeq_object = Crotophaga_results,
                                cols_samples_names = c("SW005", 
                                                       "SW006", 
                                                       "SW007", 
                                                       "SW008", 
                                                       "SW009"))

com_crot_hm <- communi_matrix(dat = data_gt_crot,
                              level_tax = "Family",
                              field = "Count",
                              tlrn = 1e-5,
                              rt = 2)

log_sums <- colSums(log(com_crot_hm + 1))

log_sums_crot <- as.data.frame(colSums(log(com_crot_hm + 1))) %>% 
  rename('log(x+1)' = 'colSums(log(com_crot_hm + 1))') %>%
  arrange(desc(`log(x+1)`))

log_sums_crot[25, ]

reduced_hm_crot <- log(com_crot_hm + 1)[, colSums(log(com_crot_hm + 1)) > 17.9] %>% 
  t()

rownames(reduced_hm_crot) <- str_replace_all(rownames(reduced_hm_crot),
                                             c("(.*)_unclassified" = "Unclassified \\1",
                                               "uncultured" = "Uncultured Bacteria",
                                               "_" = " "))

# Heatmap annotation

# Extract postmortem hours data

postmortem_hours_crot <- sample_data(Crotophaga_results)$Hours_since_death

# Define colors and legends

color_postmortem_crot <- c("1" = "#F66D7A", 
                           "2" = "#D44292", 
                           "6" = "#952EA0",
                           "48" = "#4B2991") 

legend_postmortem_crot <- c("1 hours", 
                            "2 hours", 
                            "6 hours",
                            "48 hours")

hm_annot_crot <- 
  HeatmapAnnotation("Time since death" = postmortem_hours_crot,
                    col = list("Time since death" = color_postmortem_crot),
                    annotation_legend_param = list("Time since death" = list(title = "Time since death",
                                                                             at = names(color_postmortem_crot),
                                                                             labels = legend_postmortem_crot)),
                    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
                    show_legend = FALSE,
                    annotation_name_side = "left")

# Create heatmap object

crot_heatmap <- 
  Heatmap(reduced_hm_crot,
          column_title = "B",
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
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
          cluster_columns = FALSE,
          use_raster = FALSE, 
          row_names_side = "left", 
          row_dend_side = "right",
          width = unit(35, "mm"))

# Draw heatmap

plot_crot_heatmap <- draw(crot_heatmap, 
                          heatmap_legend_side = "bottom")

# Convert heatmap to plot

plot_crot_heatmap <- heatmap_to_ggplot(plot_crot_heatmap)

plot_crot_heatmap

# Create a 'Time since death' legend ----

colors_time_hm <- c("#F6A97A",
                    "#F66D7A",
                    "#D44292", 
                    "#952EA0",
                    "#4B2991")

legend_time_hm <- c("0 hours",
                    "1 hours", 
                    "2 hours", 
                    "6 hours",
                    "48 hours")

lgd_time_death = Legend(labels = legend_time_hm, 
                        title = "Time since death", 
                        legend_gp = gpar(fill = colors_time_hm),
                        title_gp = gpar(fontsize = 10, 
                                        fontface = "bold"),
                        labels_gp = gpar(fontsize = 10),
                        legend_width = unit(0.5, "mm"))

lgd_time_death <- heatmap_to_ggplot(lgd_time_death)

lgd_time_death

# Function filter OTUs based on prevalence and abundance ----

venn_core_microbiome <- function (physeq_object, group_name, abund_value, prev_value) {
  
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
    core_microbiome <- core_members(sub, detection = abund_value, prevalence = prev_value)
    list_core[[i]] <- core_microbiome
  }
  
  return(list_core)
}

# Core Microbiome ----

# Amphisbaena bassleri ----

amph_list_core <- venn_core_microbiome(physeq_object = Amphisbaena_results,
                                       group_name = "SampleID",
                                       abund_value = 0.001,
                                       prev_value = 0.9)

plot_amph_core_microbiome <-
  venn.diagram(x = amph_list_core,
               filename = NULL,
               disable.logging = TRUE,
               fill = c("#00008B", "#FF0000", "#FFFF00", "#00FF00"),
               alpha = 0.5,
               cat.fontface = 2,
               cex = 1,
               cat.cex = 1,
               print.mode = c("raw"),
               fontface = 2)

plot_amph_core_microbiome <- ggpubr::as_ggplot(plot_amph_core_microbiome)

plot_amph_core_microbiome

amph_core_microbiome <- Reduce(intersect, amph_list_core)

amph_core_microbiome

# Crotophaga ani ----

crot_list_core <- venn_core_microbiome(physeq_object = Crotophaga_results,
                                       group_name = "SampleID",
                                       abund_value = 0.001,
                                       prev_value = 0.9)

plot_crot_core_microbiome <- 
  venn.diagram(x = crot_list_core,
               filename = NULL,
               disable.logging = TRUE,
               fill = c("#00F5FF", "#FF3E96", "#FFFF00", "#66CD00", "#68228B"),
               alpha = 0.5,
               cat.fontface = 2,
               cat.just=list(c(NA, 1), c(-0.5,-4) , c(1,-0.5) , c(NA, NA) , c(1, -3)),
               cex = 1,
               cat.cex = 1,
               print.mode = c("raw"),
               fontface = 2)

plot_crot_core_microbiome <- ggpubr::as_ggplot(plot_crot_core_microbiome)

plot_crot_core_microbiome

crot_core_microbiome <- Reduce(intersect, crot_list_core)

crot_core_microbiome

# Merge plots ----

# Relative abundance

FigureRelAbund <- grid.arrange(plot_amph_rel_abund_phylum, plot_crot_rel_abund_phylum,
                               widths = c(1, 0.05, 1),
                               layout_matrix = rbind(c(1, NA, 2)))

Fig2 <- ggpubr::as_ggplot(FigureRelAbund) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.51),
                  y = c(1, 1))

Fig2

# Heatmaps

FigureHeatmaps <- grid.arrange(plot_amph_heatmap, lgd_time_death, plot_crot_heatmap,
                               widths = c(1, 0.1, 1),
                               layout_matrix = rbind(c(1, 2, 3),
                                                     c(1, NA, 3),
                                                     c(1, NA, 3)))

Fig3 <- ggpubr::as_ggplot(FigureHeatmaps) 

Fig3

# Alpha diversity - Chao1 linear regression

FigureChao1 <- grid.arrange(plot_amph_chao1, NULL, plot_crot_chao1, 
                            nrow=3,
                            heights = c(1, 0.01, 1))

Fig4 <- ggpubr::as_ggplot(FigureChao1) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0),
                  y = c(1, 0.5))

Fig4

# Rarefaction curves

FigureRarefCurves <- grid.arrange(plot_amph_rarefaction, NULL, plot_crot_rarefaction, 
                                  nrow=3,
                                  heights = c(1, 0.01, 1))

Fig5 <- ggpubr::as_ggplot(FigureRarefCurves) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0),
                  y = c(1, 0.5))

Fig5

# PCoA

FigurePCoA <- grid.arrange(plot_amph_pcoa, plot_crot_pcoa,
                           widths = c(1, 0.05, 1),
                           layout_matrix = rbind(c(1, NA, 2)))

Fig6 <- ggpubr::as_ggplot(FigurePCoA) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.52),
                  y = c(1, 1))

Fig6

# Core microbiome

FigureCoreMicrobiome <- grid.arrange(plot_amph_core_microbiome, plot_crot_core_microbiome,
                                     widths = c(0.5, 0.05, 0.5),
                                     layout_matrix = rbind(c(1, NA, 2)))

Fig7 <- ggpubr::as_ggplot(FigureCoreMicrobiome) +
  draw_plot_label(label = LETTERS[1:2],
                  size = 12,
                  x = c(0, 0.53),
                  y = c(1, 1))

Fig7

# Save plots ----

ggsave("Plots/Fig_SeqDepth.jpeg", plot = plot_seqdepth, width = 6, height = 4, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig2_RelAbundPhyla.jpeg", plot = Fig2, width = 6, height = 3.5, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig3_HeatmapFamilies.jpeg", plot = Fig3, width = 6.75, height = 4, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig4_Chao1Regression.jpeg", plot = Fig4, width = 3.5, height = 2.75, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig5_RarefCurves.jpeg", plot = Fig5, width = 4.5, height = 3.5, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig6_PCoA.jpeg", plot = Fig6, width = 6, height = 2.5, dpi = 1000, scale = 1.5)
ggsave("Plots/Fig7_CoreMicrobiome.jpeg", plot = Fig7, width = 3, height = 1.5, dpi = 1000, scale = 3)
