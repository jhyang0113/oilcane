## Module and Data load
# Load needed tools
library(devtools)
library(vegan)
library(ggplot2)
library(reshape2)
library(corrplot)
library(wesanderson)
library(phyloseq)
library(pairwiseAdonis)
library(dplyr)
library(rstatix)
library(tidyverse)
library(knitr)
library(ggrepel)
library(ggpubr)
library(corrplot)
library(ggplot2)
library(GGally)
library(ggpubr)
library(venn)
library(pvclust)
library(microbiome)
library(berryFunctions)

# Setup working directory
require("knitr")
opts_knit$set(root.dir = "~/Desktop/R_analysis/Oilcane_data/")

# Read otu, tax, metadata files
otu_ssu <- as.matrix(read.csv(file="otu_table_ssu_no_L5.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
tax_ssu <- as.matrix(read.csv(file="tax_table_ssu.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
metadata_ssu <- read.csv(file="meta_table_ssu_no_L5.csv", stringsAsFactors=FALSE, header=TRUE)
ps_ssu <- phyloseq(otu_table(otu_ssu, taxa_are_rows=TRUE), tax_table(tax_ssu))
row.names(metadata_ssu) <- metadata_ssu$Samples
metadata_ssu <- sample_data(metadata_ssu)
sample_data(ps_ssu) <- metadata_ssu

otu_its <- as.matrix(read.csv(file="otu_table_its_edited_no_L5.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
tax_its <- as.matrix(read.csv(file="tax_table_its.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
metadata_its <- read.csv(file="meta_table_its_edited_no_L5.csv", stringsAsFactors=FALSE, header=TRUE)
ps_its <- phyloseq(otu_table(otu_its, taxa_are_rows=TRUE), tax_table(tax_its))
row.names(metadata_its) <- metadata_its$Samples
metadata_its <- sample_data(metadata_its)
sample_data(ps_its) <- metadata_its

# Check the input number of ASVs before pretreatment
ps_ssu
ps_its

## Data pretreatment
# Filtering the raw data step 1: read counts per sample
sample_sum_df_ssu <- data.frame(sum = sample_sums(ps_ssu))
ggplot(sample_sum_df_ssu, aes(x = sum)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

sample_sum_df_its <- data.frame(sum = sample_sums(ps_its))
ggplot(sample_sum_df_its, aes(x = sum)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Filtering the raw data step 2: compute feature prevalence in a dataframe
prevdf_ssu = apply(X = otu_table(ps_ssu),
               MARGIN = ifelse(taxa_are_rows(ps_ssu), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf_its = apply(X = otu_table(ps_its),
               MARGIN = ifelse(taxa_are_rows(ps_its), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Filtering the raw data step 3: add taxonomy
prevdf_ssu = data.frame(Prevalence = prevdf_ssu,
                    TotalAbundance = taxa_sums(ps_ssu), tax_table(ps_ssu))
plyr::ddply(prevdf_ssu, "Phylum", function(df1_ssu){cbind(mean(df1_ssu$Prevalence),sum(df1_ssu$Prevalence))}) -> dfprev_ssu
kable(dfprev_ssu)

prevdf_its = data.frame(Prevalence = prevdf_its,
                    TotalAbundance = taxa_sums(ps_its), tax_table(ps_its))
plyr::ddply(prevdf_its, "Phylum", function(df1_its){cbind(mean(df1_its$Prevalence),sum(df1_its$Prevalence))}) -> dfprev_its
kable(dfprev_its)

# Filtering the raw data step 4: remove singleton
filterPhyla_ssu = c("Crenarchaeota","Thermoplasmatota", "Nanoarchaeota","Halobacterota", "Euryarchaeota", "Micrarchaeota", "Nanohaloarchaeota", "Aenigmarchaeota", NA)
(ps1_ssu = subset_taxa(ps_ssu, !Phylum %in% filterPhyla_ssu))
filterPhyla2_ssu <- c("Eukaryota", "Chloroplast")
ps1_ssu <- subset_taxa(ps1_ssu, !Kingdom %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Phylum %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Class %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Order %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Family %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Genus %in% filterPhyla2_ssu)
ps1_ssu

filterPhyla_its = c(NA)
(ps1_its = subset_taxa(ps_its, !Phylum %in% filterPhyla_its))
filterPhyla2_its <- c(NA)
ps1_its <- subset_taxa(ps1_its, !Kingdom %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Phylum %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Class %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Order %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Family %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Genus %in% filterPhyla2_its)
ps1_its

# Filtering the raw data step 8: read counts per filtered sample
sample_sum_df_2_ssu <- data.frame(sum_2_ssu = sample_sums(ps1_ssu))
ggplot(sample_sum_df_2_ssu, aes(x = sum_2_ssu)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

sample_sum_df_2_its <- data.frame(sum_2_its = sample_sums(ps1_its))
ggplot(sample_sum_df_2_its, aes(x = sum_2_its)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# Filtering the raw data step 9: make new data frame with new phyloseq objective
metadata_ps1_ssu <- sample_data(metadata_ssu)
sample_data(ps1_ssu) <- metadata_ssu
sampledf_ps1_ssu <- data.frame(sample_data(ps1_ssu))

metadata_ps1_its <- sample_data(metadata_its)
sample_data(ps1_its) <- metadata_its
sampledf_ps1_its <- data.frame(sample_data(ps1_its))

# Subsampling (each location and compartment)
# By location
ps1_ssu_g <- subset_samples(ps1_ssu, loc.time == "greenhouse")
sampledf_ps1_ssu_g <- data.frame(ps1_ssu_g@sam_data)
ps1_ssu_g
ps1_ssu_g_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g) > 0, ps1_ssu_g)
sampledf_ps1_ssu_g_filtered <- data.frame(ps1_ssu_g_filtered@sam_data)
ps1_ssu_g_filtered

ps1_ssu_f <- subset_samples(ps1_ssu, loc.time == "field")
sampledf_ps1_ssu_f <- data.frame(ps1_ssu_f@sam_data)
ps1_ssu_f
ps1_ssu_f_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f) > 0, ps1_ssu_f)
sampledf_ps1_ssu_f_filtered <- data.frame(ps1_ssu_f_filtered@sam_data)
ps1_ssu_f_filtered

ps1_its_g <- subset_samples(ps1_its, loc.time == "greenhouse")
sampledf_ps1_its_g <- data.frame(ps1_its_g@sam_data)
ps1_its_g
ps1_its_g_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g) > 0, ps1_its_g)
sampledf_ps1_its_g_filtered <- data.frame(ps1_its_g_filtered@sam_data)
ps1_its_g_filtered

ps1_its_f <- subset_samples(ps1_its, loc.time == "field")
sampledf_ps1_its_f <- data.frame(ps1_its_f@sam_data)
ps1_its_f
ps1_its_f_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f) > 0, ps1_its_f)
sampledf_ps1_its_f_filtered <- data.frame(ps1_its_f_filtered@sam_data)
ps1_its_f_filtered

# by compartment
ps1_ssu_g_leaf <- subset_samples(ps1_ssu_g, type == "leaf")
sampledf_ps1_ssu_g_leaf <- data.frame(ps1_ssu_g_leaf@sam_data)
ps1_ssu_g_leaf
ps1_ssu_g_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_leaf) > 0, ps1_ssu_g_leaf)
sampledf_ps1_ssu_g_leaf_filtered <- data.frame(ps1_ssu_g_leaf_filtered@sam_data)
ps1_ssu_g_leaf_filtered

ps1_ssu_g_stem <- subset_samples(ps1_ssu_g, type == "stem")
sampledf_ps1_ssu_g_stem <- data.frame(ps1_ssu_g_stem@sam_data)
ps1_ssu_g_stem
ps1_ssu_g_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_stem) > 0, ps1_ssu_g_stem)
sampledf_ps1_ssu_g_stem_filtered <- data.frame(ps1_ssu_g_stem_filtered@sam_data)
ps1_ssu_g_stem_filtered

ps1_ssu_g_root <- subset_samples(ps1_ssu_g, type == "root")
sampledf_ps1_ssu_g_root <- data.frame(ps1_ssu_g_root@sam_data)
ps1_ssu_g_root
ps1_ssu_g_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_root) > 0, ps1_ssu_g_root)
sampledf_ps1_ssu_g_root_filtered <- data.frame(ps1_ssu_g_root_filtered@sam_data)
ps1_ssu_g_root_filtered

ps1_ssu_g_root_associated_soils <- subset_samples(ps1_ssu_g, type == "root associated soils")
sampledf_ps1_ssu_g_root_associated_soils <- data.frame(ps1_ssu_g_root_associated_soils@sam_data)
ps1_ssu_g_root_associated_soils
ps1_ssu_g_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_root_associated_soils) > 0, ps1_ssu_g_root_associated_soils)
sampledf_ps1_ssu_g_root_associated_soils_filtered <- data.frame(ps1_ssu_g_root_associated_soils_filtered@sam_data)
ps1_ssu_g_root_associated_soils_filtered

ps1_ssu_g_soil <- subset_samples(ps1_ssu_g, type == "soil")
sampledf_ps1_ssu_g_soil <- data.frame(ps1_ssu_g_soil@sam_data)
ps1_ssu_g_soil
ps1_ssu_g_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_soil) > 0, ps1_ssu_g_soil)
sampledf_ps1_ssu_g_soil_filtered <- data.frame(ps1_ssu_g_soil_filtered@sam_data)
ps1_ssu_g_soil_filtered

ps1_ssu_f_leaf <- subset_samples(ps1_ssu_f, type == "leaf")
sampledf_ps1_ssu_f_leaf <- data.frame(ps1_ssu_f_leaf@sam_data)
ps1_ssu_f_leaf
ps1_ssu_f_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_leaf) > 0, ps1_ssu_f_leaf)
sampledf_ps1_ssu_f_leaf_filtered <- data.frame(ps1_ssu_f_leaf_filtered@sam_data)
ps1_ssu_f_leaf_filtered

ps1_ssu_f_stem <- subset_samples(ps1_ssu_f, type == "stem")
sampledf_ps1_ssu_f_stem <- data.frame(ps1_ssu_f_stem@sam_data)
ps1_ssu_f_stem
ps1_ssu_f_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_stem) > 0, ps1_ssu_f_stem)
sampledf_ps1_ssu_f_stem_filtered <- data.frame(ps1_ssu_f_stem_filtered@sam_data)
ps1_ssu_f_stem_filtered

ps1_ssu_f_root <- subset_samples(ps1_ssu_f, type == "root")
sampledf_ps1_ssu_f_root <- data.frame(ps1_ssu_f_root@sam_data)
ps1_ssu_f_root
ps1_ssu_f_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_root) > 0, ps1_ssu_f_root)
sampledf_ps1_ssu_f_root_filtered <- data.frame(ps1_ssu_f_root_filtered@sam_data)
ps1_ssu_f_root_filtered

ps1_ssu_f_root_associated_soils <- subset_samples(ps1_ssu_f, type == "root associated soils")
sampledf_ps1_ssu_f_root_associated_soils <- data.frame(ps1_ssu_f_root_associated_soils@sam_data)
ps1_ssu_f_root_associated_soils
ps1_ssu_f_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_root_associated_soils) > 0, ps1_ssu_f_root_associated_soils)
sampledf_ps1_ssu_f_root_associated_soils_filtered <- data.frame(ps1_ssu_f_root_associated_soils_filtered@sam_data)
ps1_ssu_f_root_associated_soils_filtered

ps1_ssu_f_soil <- subset_samples(ps1_ssu_f, type == "soil")
sampledf_ps1_ssu_f_soil <- data.frame(ps1_ssu_f_soil@sam_data)
ps1_ssu_f_soil
ps1_ssu_f_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_soil) > 0, ps1_ssu_f_soil)
sampledf_ps1_ssu_f_soil_filtered <- data.frame(ps1_ssu_f_soil_filtered@sam_data)
ps1_ssu_f_soil_filtered

ps1_its_g_leaf <- subset_samples(ps1_its_g, type == "leaf")
sampledf_ps1_its_g_leaf <- data.frame(ps1_its_g_leaf@sam_data)
ps1_its_g_leaf
ps1_its_g_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_leaf) > 0, ps1_its_g_leaf)
sampledf_ps1_its_g_leaf_filtered <- data.frame(ps1_its_g_leaf_filtered@sam_data)
ps1_its_g_leaf_filtered

ps1_its_g_stem <- subset_samples(ps1_its_g, type == "stem")
sampledf_ps1_its_g_stem <- data.frame(ps1_its_g_stem@sam_data)
ps1_its_g_stem
ps1_its_g_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_stem) > 0, ps1_its_g_stem)
sampledf_ps1_its_g_stem_filtered <- data.frame(ps1_its_g_stem_filtered@sam_data)
ps1_its_g_stem_filtered

ps1_its_g_root <- subset_samples(ps1_its_g, type == "root")
sampledf_ps1_its_g_root <- data.frame(ps1_its_g_root@sam_data)
ps1_its_g_root
ps1_its_g_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_root) > 0, ps1_its_g_root)
sampledf_ps1_its_g_root_filtered <- data.frame(ps1_its_g_root_filtered@sam_data)
ps1_its_g_root_filtered

ps1_its_g_root_associated_soils <- subset_samples(ps1_its_g, type == "root associated soils")
sampledf_ps1_its_g_root_associated_soils <- data.frame(ps1_its_g_root_associated_soils@sam_data)
ps1_its_g_root_associated_soils
ps1_its_g_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_root_associated_soils) > 0, ps1_its_g_root_associated_soils)
sampledf_ps1_its_g_root_associated_soils_filtered <- data.frame(ps1_its_g_root_associated_soils_filtered@sam_data)
ps1_its_g_root_associated_soils_filtered

ps1_its_g_soil <- subset_samples(ps1_its_g, type == "soil")
sampledf_ps1_its_g_soil <- data.frame(ps1_its_g_soil@sam_data)
ps1_its_g_soil
ps1_its_g_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_soil) > 0, ps1_its_g_soil)
sampledf_ps1_its_g_soil_filtered <- data.frame(ps1_its_g_soil_filtered@sam_data)
ps1_its_g_soil_filtered

ps1_its_f_leaf <- subset_samples(ps1_its_f, type == "leaf")
sampledf_ps1_its_f_leaf <- data.frame(ps1_its_f_leaf@sam_data)
ps1_its_f_leaf
ps1_its_f_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_leaf) > 0, ps1_its_f_leaf)
sampledf_ps1_its_f_leaf_filtered <- data.frame(ps1_its_f_leaf_filtered@sam_data)
ps1_its_f_leaf_filtered

ps1_its_f_stem <- subset_samples(ps1_its_f, type == "stem")
sampledf_ps1_its_f_stem <- data.frame(ps1_its_f_stem@sam_data)
ps1_its_f_stem
ps1_its_f_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_stem) > 0, ps1_its_f_stem)
sampledf_ps1_its_f_stem_filtered <- data.frame(ps1_its_f_stem_filtered@sam_data)
ps1_its_f_stem_filtered

ps1_its_f_root <- subset_samples(ps1_its_f, type == "root")
sampledf_ps1_its_f_root <- data.frame(ps1_its_f_root@sam_data)
ps1_its_f_root
ps1_its_f_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_root) > 0, ps1_its_f_root)
sampledf_ps1_its_f_root_filtered <- data.frame(ps1_its_f_root_filtered@sam_data)
ps1_its_f_root_filtered

ps1_its_f_root_associated_soils <- subset_samples(ps1_its_f, type == "root associated soils")
sampledf_ps1_its_f_root_associated_soils <- data.frame(ps1_its_f_root_associated_soils@sam_data)
ps1_its_f_root_associated_soils
ps1_its_f_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_root_associated_soils) > 0, ps1_its_f_root_associated_soils)
sampledf_ps1_its_f_root_associated_soils_filtered <- data.frame(ps1_its_f_root_associated_soils_filtered@sam_data)
ps1_its_f_root_associated_soils_filtered

ps1_its_f_soil <- subset_samples(ps1_its_f, type == "soil")
sampledf_ps1_its_f_soil <- data.frame(ps1_its_f_soil@sam_data)
ps1_its_f_soil
ps1_its_f_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_soil) > 0, ps1_its_f_soil)
sampledf_ps1_its_f_soil_filtered <- data.frame(ps1_its_f_soil_filtered@sam_data)
ps1_its_f_soil_filtered

## Q1. Were bacterial and fungal microbiomes significantly different according to the variation of plant accessions? 
# dataset
ps1_ssu_g_leaf_filtered
ps1_ssu_g_stem_filtered
ps1_ssu_g_root_filtered
ps1_ssu_g_root_associated_soils_filtered
ps1_ssu_g_soil_filtered
ps1_its_g_leaf_filtered
ps1_its_g_stem_filtered
ps1_its_g_root_filtered
ps1_its_g_root_associated_soils_filtered
ps1_its_g_soil_filtered

# Create the function for following analysis
analysis_location <- function (dataset_x) {
  
normalization_dataset_x = transform_sample_counts(dataset_x, function(x) x/sum(x))

# 1. Remove data points with missing metadata
phy_dataset_x <- normalization_dataset_x %>%
  subset_samples(
      !is.na(genotype) &
      !is.na(type2) &
      !is.na(loc.time) &
      !is.na(id))

# 2. Differences depending on input variables
bray_dataset_x <- phyloseq::distance(physeq = phy_dataset_x, method = "bray")
cap_ord_dataset_x <- ordinate(physeq = phy_dataset_x, method = "CAP", distance = bray_dataset_x, formula = ~ type2 + genotype)

# 3. ANOVA
obj_phy_dataset_x <- t(phy_dataset_x)
anova(cap_ord_dataset_x, by="terms")
vif.cca(cap_ord_dataset_x)

# 4. CAP plot
cap_plot_dataset_x <- plot_ordination(physeq = phy_dataset_x, 
                                ordination = cap_ord_dataset_x,
                                color = "genotype", 
                                axes = c(1,2)) +  xlim(-2, 2) + ylim(-2, 2) + 
                                geom_point(aes(colour = genotype), 
                                           alpha = 0.8, size = 3) + geom_point(size = 3) +
                                scale_color_manual(
                                  values = c("#00AFBB","#FC4E07", "#E7B800",
                                             "#36BB07", "#A04000",
                                             "#757575")) + theme_bw()
cap_plot_dataset_x_2 <- cap_plot_dataset_x
cap_plot_dataset_x_2$data$genotype<-factor(cap_plot_dataset_x$data$genotype, levels = c("WT", "L1", "L2", "L3", "L4"))
print(cap_plot_dataset_x_2)

cap_plot_dataset_x_3 <- plot_ordination(physeq = phy_dataset_x, 
                                ordination = cap_ord_dataset_x,
                                color = "type2", 
                                axes = c(1,2)) +  xlim(-2, 2) + ylim(-2, 2) + 
                                geom_point(aes(colour = type2), 
                                           alpha = 0.8, size = 3) + geom_point(size = 3) +
                                scale_color_manual(
                                  values = c("#00FF00","#006600", "#D29E79",
                                             "#D35400", "#6E2C00",
                                             "#757575")) + theme_bw()
cap_plot_dataset_x_4 <- cap_plot_dataset_x_3
cap_plot_dataset_x_4$data$type2<-factor(cap_plot_dataset_x_3$data$type2, levels = c("LF", "ST", "RT", "RS", "SO"))
print(cap_plot_dataset_x_4)

# 5. Create distance matrix for PERMANOVA
dist_phy_dataset_x <-vegdist(decostand(obj_phy_dataset_x@otu_table,"hellinger"),"bray")

# 6. transform in data-frame
sampledf_dist_phy_dataset_x <- data.frame(sample_data(obj_phy_dataset_x))

# 7. PERMANOVA (pairwise_adonis)
adonis_table_dataset_x <- adonis(obj_phy_dataset_x@otu_table ~ genotype, data = sampledf_dist_phy_dataset_x, permutations = 999)
adonis_table_dataset_x
adonis_result_dataset_x <- subset(adonis_table_dataset_x$aov.tab, select = c("R2", "Pr(>F)"))
print(adonis_result_dataset_x)

pairwiseadonis_result_dataset_x <- pairwise.adonis(obj_phy_dataset_x@otu_table, sample_data(obj_phy_dataset_x)$genotype, p.adjust.m = "fdr", perm = 999)
print(pairwiseadonis_result_dataset_x)

adonis_table_dataset_x_2 <- adonis(obj_phy_dataset_x@otu_table ~ type2, data = sampledf_dist_phy_dataset_x, permutations = 999)
adonis_table_dataset_x_2
adonis_result_dataset_x_2 <- subset(adonis_table_dataset_x_2$aov.tab, select = c("R2", "Pr(>F)"))
print(adonis_result_dataset_x_2)

pairwiseadonis_result_dataset_x_2 <- pairwise.adonis(obj_phy_dataset_x@otu_table, sample_data(obj_phy_dataset_x)$type2, p.adjust.m = "fdr", perm = 999)
print(pairwiseadonis_result_dataset_x_2)

adonis_table_dataset_x_3 <- adonis(obj_phy_dataset_x@otu_table ~ genotype*type, data = sampledf_dist_phy_dataset_x, permutations = 999)
adonis_table_dataset_x_3
adonis_result_dataset_x_3 <- subset(adonis_table_dataset_x_3$aov.tab, select = c("R2", "Pr(>F)"))
print(adonis_result_dataset_x_3)

# 8. Betadiversity / dispersion / homogeneity between groups
bet_factor_dataset_x <-betadisper(dist_phy_dataset_x, sampledf_dist_phy_dataset_x$genotype)
betdisper_result_dataset_x <- permutest(bet_factor_dataset_x, pairwise = TRUE, permutations = 999, group=bet_factor_dataset_x$group)
print(betdisper_result_dataset_x)

bet_factor_dataset_x_2 <-betadisper(dist_phy_dataset_x, sampledf_dist_phy_dataset_x$type2)
betdisper_result_dataset_x_2 <- permutest(bet_factor_dataset_x_2, pairwise = TRUE, permutations = 999, group=bet_factor_dataset_x_2$group)
print(betdisper_result_dataset_x_2)

# 9. Alpha diversity plot
diversity_dataset_dataset_x <- prune_taxa(taxa_sums(dataset_x) > 0, dataset_x)
diversity_plot_dataset_x <- plot_richness(diversity_dataset_dataset_x, x="genotype",measures=c("Chao1", "Shannon"), color = "genotype") + theme_bw()
diversity_plot_dataset_x + geom_boxplot(data = diversity_plot_dataset_x$data, aes(x = genotype, y = value, color = NULL), alpha = 0.1)
#write.csv(diversity_plot_dataset_x$plot_env$DF, file="alpha_diversity_ssu.csv")
diversity_plot_dataset_x_data <- diversity_plot_dataset_x$plot_env$DF
ordered_diversity_plot_dataset_x_data <- reorder_levels(diversity_plot_dataset_x_data, genotype, c("WT", "L1", "L2", "L3", "L4"))

res.kruskal_total_shannon_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  kruskal_test(Shannon ~ genotype)
print(res.kruskal_total_shannon_dataset_x)

pwc_total_shannon_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  dunn_test(Shannon ~ genotype, p.adjust.method = "bonferroni") 
print(pwc_total_shannon_dataset_x)

pwc_total_shannon_dataset_x <- pwc_total_shannon_dataset_x %>% add_xy_position()
shannon_plot <- ggboxplot(ordered_diversity_plot_dataset_x_data, x = "genotype", y = "Shannon", color = "genotype", 
          palette = c("#00AFBB","#FC4E07", "#E7B800", "#36BB07", "#7D26CD", "#757575"),
          order = c("WT", "L1", "L2", "L3", "L4"),
          ylab = "Shannon index", xlab = "genotype") +
          stat_pvalue_manual(pwc_total_shannon_dataset_x, y.position = 8.0,
          step.increase = 0.1, hide.ns = TRUE) + 
          labs(subtitle = get_test_label(res.kruskal_total_shannon_dataset_x, detailed = TRUE),
          caption = get_pwc_label(pwc_total_shannon_dataset_x)) + coord_cartesian(ylim = c(0, 15))
print(shannon_plot)

res.kruskal_total_chao1_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  kruskal_test(Chao1 ~ genotype)
print(res.kruskal_total_chao1_dataset_x)

pwc_total_chao1_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  dunn_test(Chao1 ~ genotype, p.adjust.method = "bonferroni") 
print(pwc_total_chao1_dataset_x)

pwc_total_chao1_dataset_x <- pwc_total_chao1_dataset_x %>% add_xy_position()
chao_plot <- ggboxplot(ordered_diversity_plot_dataset_x_data, x = "genotype", y = "Chao1", color = "genotype", 
          palette = c("#00AFBB","#FC4E07", "#E7B800", "#36BB07", "#7D26CD", "#757575"),
          order = c("WT", "L1", "L2", "L3", "L4"),
          ylab = "chao1 index", xlab = "genotype") +
          stat_pvalue_manual(pwc_total_chao1_dataset_x, y.position = 2300,
          step.increase = 0.1, hide.ns = TRUE) + 
          labs(subtitle = get_test_label(res.kruskal_total_chao1_dataset_x, detailed = TRUE),
          caption = get_pwc_label(pwc_total_chao1_dataset_x)) + coord_cartesian(ylim = c(0, 4500))
print(chao_plot)}

analysis_location(ps1_ssu_g_filtered)
analysis_location(ps1_its_g_filtered)
