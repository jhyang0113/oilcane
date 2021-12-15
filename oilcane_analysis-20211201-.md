Untitled
================

## Module and Data load

# Load needed tools

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(vegan)
```

    ## Loading required package: permute

    ## 
    ## Attaching package: 'permute'

    ## The following object is masked from 'package:devtools':
    ## 
    ##     check

    ## Loading required package: lattice

    ## This is vegan 2.5-7

``` r
library(ggplot2)
library(reshape2)
library(corrplot)
```

    ## corrplot 0.90 loaded

``` r
library(wesanderson)
library(phyloseq)
library(pairwiseAdonis)
```

    ## Loading required package: cluster

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(rstatix)
```

    ## 
    ## Attaching package: 'rstatix'

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✓ tibble  3.1.6     ✓ purrr   0.3.4
    ## ✓ tidyr   1.1.4     ✓ stringr 1.4.0
    ## ✓ readr   2.0.2     ✓ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## x rstatix::filter() masks dplyr::filter(), stats::filter()
    ## x dplyr::lag()      masks stats::lag()

``` r
library(knitr)
library(ggrepel)
library(ggpubr)
library(corrplot)
library(ggplot2)
library(GGally)
```

    ## Registered S3 method overwritten by 'GGally':
    ##   method from   
    ##   +.gg   ggplot2

``` r
library(ggpubr)
library(venn)
```

# Setup working directory

``` r
require("knitr")
opts_knit$set(root.dir = "~/Desktop/R_analysis/Oilcane_data/")
```

# Read otu, tax, metadata files - ASV-based (for prevalence-abundance figure and SSU microbial community)

``` r
otu_ssu <- as.matrix(read.csv(file="otu_table_ssu.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
tax_ssu <- as.matrix(read.csv(file="tax_table_ssu.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
metadata_ssu <- read.csv(file="meta_table_ssu.csv", stringsAsFactors=FALSE, header=TRUE)
ps_ssu <- phyloseq(otu_table(otu_ssu, taxa_are_rows=TRUE), tax_table(tax_ssu))
row.names(metadata_ssu) <- metadata_ssu$Samples
metadata_ssu <- sample_data(metadata_ssu)
sample_data(ps_ssu) <- metadata_ssu

otu_its <- as.matrix(read.csv(file="otu_table_its_edited.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
tax_its <- as.matrix(read.csv(file="tax_table_its.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
metadata_its <- read.csv(file="meta_table_its_edited.csv", stringsAsFactors=FALSE, header=TRUE)
ps_its <- phyloseq(otu_table(otu_its, taxa_are_rows=TRUE), tax_table(tax_its))
row.names(metadata_its) <- metadata_its$Samples
metadata_its <- sample_data(metadata_its)
sample_data(ps_its) <- metadata_its
```

# Check the input number of ASVs before pretreatment

``` r
ps_ssu
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 49275 taxa and 238 samples ]
    ## sample_data() Sample Data:       [ 238 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 49275 taxa by 7 taxonomic ranks ]

``` r
ps_its
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 10548 taxa and 232 samples ]
    ## sample_data() Sample Data:       [ 232 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 10548 taxa by 7 taxonomic ranks ]

## Data pretreatment

# Filtering the raw data step 1: read counts per sample

``` r
sample_sum_df_ssu <- data.frame(sum = sample_sums(ps_ssu))
ggplot(sample_sum_df_ssu, aes(x = sum)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
sample_sum_df_its <- data.frame(sum = sample_sums(ps_its))
ggplot(sample_sum_df_its, aes(x = sum)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

# Filtering the raw data step 2: compute feature prevalence in a dataframe

``` r
prevdf_ssu = apply(X = otu_table(ps_ssu),
               MARGIN = ifelse(taxa_are_rows(ps_ssu), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf_its = apply(X = otu_table(ps_its),
               MARGIN = ifelse(taxa_are_rows(ps_its), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
```

# Filtering the raw data step 3: add taxonomy

``` r
prevdf_ssu = data.frame(Prevalence = prevdf_ssu,
                    TotalAbundance = taxa_sums(ps_ssu), tax_table(ps_ssu))
plyr::ddply(prevdf_ssu, "Phylum", function(df1_ssu){cbind(mean(df1_ssu$Prevalence),sum(df1_ssu$Prevalence))}) -> dfprev_ssu
kable(dfprev_ssu)
```

| Phylum                       |         1 |     2 |
|:-----------------------------|----------:|------:|
| Abditibacteriota             |  2.333333 |    91 |
| Acidobacteriota              |  3.187854 | 10708 |
| Actinobacteriota             |  3.928723 | 22213 |
| Aenigmarchaeota              |  1.000000 |     1 |
| AncK6                        |  1.000000 |     1 |
| Armatimonadota               |  2.158076 |  1256 |
| Bacteroidota                 |  2.815886 |  5353 |
| Bdellovibrionota             |  1.637911 |  1443 |
| Calditrichota                |  1.750000 |    35 |
| Campilobacterota             |  1.000000 |     1 |
| Chloroflexi                  |  3.172982 | 18086 |
| Crenarchaeota                | 10.540146 |  1444 |
| Cyanobacteria                |  2.848120 |  1894 |
| Dadabacteria                 |  3.526316 |    67 |
| Deferrisomatota              |  1.695652 |    39 |
| Deinococcota                 |  1.375000 |    22 |
| Dependentiae                 |  1.497942 |   364 |
| Desulfobacterota             |  2.494024 |  1252 |
| DTB120                       |  1.666667 |    10 |
| Elusimicrobiota              |  1.117647 |   133 |
| Entotheonellaeota            |  3.649007 |   551 |
| Euryarchaeota                |  1.571429 |    11 |
| FCPU426                      |  1.133333 |    17 |
| Fibrobacterota               |  2.530303 |   167 |
| Firmicutes                   |  2.977735 |  7757 |
| FW113                        |  1.000000 |     1 |
| GAL15                        |  3.676471 |   125 |
| Gemmatimonadota              |  3.695000 |  5912 |
| Halanaerobiaeota             |  1.000000 |     4 |
| Halobacterota                |  1.500000 |    15 |
| Hydrogenedentes              |  1.905660 |   101 |
| Latescibacterota             |  2.948276 |   855 |
| Margulisbacteria             |  1.000000 |     1 |
| MBNT15                       |  3.228571 |   452 |
| Methylomirabilota            |  5.447552 |  1558 |
| Micrarchaeota                |  2.000000 |     2 |
| Myxococcota                  |  2.345082 |  5722 |
| Nanoarchaeota                |  1.255814 |    54 |
| Nanohaloarchaeota            |  2.000000 |     2 |
| NB1-j                        |  4.356250 |   697 |
| Nitrospinota                 |  3.000000 |     9 |
| Nitrospirota                 |  4.213592 |   868 |
| Patescibacteria              |  1.358757 |   481 |
| PAUC34f                      |  1.000000 |     2 |
| Planctomycetota              |  2.130767 | 13117 |
| Proteobacteria               |  3.855570 | 40630 |
| RCP2-54                      |  7.133333 |   321 |
| SAR324 clade(Marine group B) |  1.428571 |    70 |
| Spirochaetota                |  1.428571 |    50 |
| Sumerlaeota                  |  1.461539 |    95 |
| Sva0485                      |  1.642857 |    23 |
| Thermoplasmatota             |  1.250000 |    35 |
| TX1A-33                      |  1.333333 |     4 |
| Verrucomicrobiota            |  2.134490 |  2952 |
| WPS-2                        |  2.684210 |    51 |
| WS1                          |  1.000000 |     2 |
| WS2                          |  1.976190 |    83 |
| WS4                          |  1.300000 |    13 |
| Zixibacteria                 |  2.333333 |    49 |
| NA                           |  1.321046 |  3333 |

``` r
prevdf_its = data.frame(Prevalence = prevdf_its,
                    TotalAbundance = taxa_sums(ps_its), tax_table(ps_its))
plyr::ddply(prevdf_its, "Phylum", function(df1_its){cbind(mean(df1_its$Prevalence),sum(df1_its$Prevalence))}) -> dfprev_its
kable(dfprev_its)
```

| Phylum                      |        1 |     2 |
|:----------------------------|---------:|------:|
| p\_\_Aphelidiomycota        | 1.461539 |    38 |
| p\_\_Ascomycota             | 4.567918 | 16545 |
| p\_\_Basidiobolomycota      | 1.000000 |     4 |
| p\_\_Basidiomycota          | 2.395833 |  1495 |
| p\_\_Blastocladiomycota     | 1.545454 |    17 |
| p\_\_Calcarisporiellomycota | 1.500000 |     3 |
| p\_\_Chytridiomycota        | 1.977778 |    89 |
| p\_\_Entomophthoromycota    | 1.000000 |     1 |
| p\_\_Glomeromycota          | 1.451613 |   225 |
| p\_\_Kickxellomycota        | 1.000000 |     5 |
| p\_\_Monoblepharomycota     | 4.750000 |    38 |
| p\_\_Mortierellomycota      | 4.477778 |   403 |
| p\_\_Mucoromycota           | 4.764706 |   243 |
| p\_\_Olpidiomycota          | 1.000000 |     1 |
| p\_\_Rozellomycota          | 3.252874 |   283 |
| p\_\_Zoopagomycota          | 4.000000 |    20 |
| NA                          | 3.168646 | 18413 |

# Filtering the raw data step 4: remove singleton

``` r
filterPhyla_ssu = c("Crenarchaeota","Thermoplasmatota", "Nanoarchaeota","Halobacterota", "Euryarchaeota", "Micrarchaeota", "Nanohaloarchaeota", "Aenigmarchaeota", NA)
(ps1_ssu = subset_taxa(ps_ssu, !Phylum %in% filterPhyla_ssu))
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 238 samples ]
    ## sample_data() Sample Data:       [ 238 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
filterPhyla2_ssu <- c("Eukaryota")
ps1_ssu <- subset_taxa(ps1_ssu, !Kingdom %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Phylum %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Class %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Order %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Family %in% filterPhyla2_ssu)
ps1_ssu <- subset_taxa(ps1_ssu, !Genus %in% filterPhyla2_ssu)
ps1_ssu
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 238 samples ]
    ## sample_data() Sample Data:       [ 238 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
filterPhyla_its = c(NA)
(ps1_its = subset_taxa(ps_its, !Phylum %in% filterPhyla_its))
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 4737 taxa and 232 samples ]
    ## sample_data() Sample Data:       [ 232 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 4737 taxa by 7 taxonomic ranks ]

``` r
filterPhyla2_its <- c(NA)
ps1_its <- subset_taxa(ps1_its, !Kingdom %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Phylum %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Class %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Order %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Family %in% filterPhyla2_its)
ps1_its <- subset_taxa(ps1_its, !Genus %in% filterPhyla2_its)
ps1_its
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 232 samples ]
    ## sample_data() Sample Data:       [ 232 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

# Filtering the raw data step 8: read counts per filtered sample

``` r
sample_sum_df_2_ssu <- data.frame(sum_2_ssu = sample_sums(ps1_ssu))
ggplot(sample_sum_df_2_ssu, aes(x = sum_2_ssu)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
sample_sum_df_2_its <- data.frame(sum_2_its = sample_sums(ps1_its))
ggplot(sample_sum_df_2_its, aes(x = sum_2_its)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") +
  xlab("Read counts") +
  theme(axis.title.y = element_blank())
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

# Filtering the raw data step 9: make new data frame with new phyloseq objective

``` r
metadata_ps1_ssu <- sample_data(metadata_ssu)
sample_data(ps1_ssu) <- metadata_ssu
sampledf_ps1_ssu <- data.frame(sample_data(ps1_ssu))

metadata_ps1_its <- sample_data(metadata_its)
sample_data(ps1_its) <- metadata_its
sampledf_ps1_its <- data.frame(sample_data(ps1_its))
```

# Subsampling (each location and compartment)

``` r
# By location
ps1_ssu_g <- subset_samples(ps1_ssu, loc.time == "greenhouse")
sampledf_ps1_ssu_g <- data.frame(ps1_ssu_g@sam_data)
ps1_ssu_g
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 118 samples ]
    ## sample_data() Sample Data:       [ 118 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g) > 0, ps1_ssu_g)
sampledf_ps1_ssu_g_filtered <- data.frame(ps1_ssu_g_filtered@sam_data)
ps1_ssu_g_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 23716 taxa and 118 samples ]
    ## sample_data() Sample Data:       [ 118 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 23716 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f <- subset_samples(ps1_ssu, loc.time == "field")
sampledf_ps1_ssu_f <- data.frame(ps1_ssu_f@sam_data)
ps1_ssu_f
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 120 samples ]
    ## sample_data() Sample Data:       [ 120 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f) > 0, ps1_ssu_f)
sampledf_ps1_ssu_f_filtered <- data.frame(ps1_ssu_f_filtered@sam_data)
ps1_ssu_f_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 29750 taxa and 120 samples ]
    ## sample_data() Sample Data:       [ 120 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 29750 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g <- subset_samples(ps1_its, loc.time == "greenhouse")
sampledf_ps1_its_g <- data.frame(ps1_its_g@sam_data)
ps1_its_g
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 115 samples ]
    ## sample_data() Sample Data:       [ 115 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g) > 0, ps1_its_g)
sampledf_ps1_its_g_filtered <- data.frame(ps1_its_g_filtered@sam_data)
ps1_its_g_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1772 taxa and 115 samples ]
    ## sample_data() Sample Data:       [ 115 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1772 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f <- subset_samples(ps1_its, loc.time == "field")
sampledf_ps1_its_f <- data.frame(ps1_its_f@sam_data)
ps1_its_f
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 117 samples ]
    ## sample_data() Sample Data:       [ 117 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f) > 0, ps1_its_f)
sampledf_ps1_its_f_filtered <- data.frame(ps1_its_f_filtered@sam_data)
ps1_its_f_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2170 taxa and 117 samples ]
    ## sample_data() Sample Data:       [ 117 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2170 taxa by 7 taxonomic ranks ]

``` r
# by compartment
ps1_ssu_g_leaf <- subset_samples(ps1_ssu_g, type == "leaf")
sampledf_ps1_ssu_g_leaf <- data.frame(ps1_ssu_g_leaf@sam_data)
ps1_ssu_g_leaf
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_leaf) > 0, ps1_ssu_g_leaf)
sampledf_ps1_ssu_g_leaf_filtered <- data.frame(ps1_ssu_g_leaf_filtered@sam_data)
ps1_ssu_g_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 254 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 254 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_stem <- subset_samples(ps1_ssu_g, type == "stem")
sampledf_ps1_ssu_g_stem <- data.frame(ps1_ssu_g_stem@sam_data)
ps1_ssu_g_stem
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_stem) > 0, ps1_ssu_g_stem)
sampledf_ps1_ssu_g_stem_filtered <- data.frame(ps1_ssu_g_stem_filtered@sam_data)
ps1_ssu_g_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 377 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 377 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root <- subset_samples(ps1_ssu_g, type == "root")
sampledf_ps1_ssu_g_root <- data.frame(ps1_ssu_g_root@sam_data)
ps1_ssu_g_root
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_root) > 0, ps1_ssu_g_root)
sampledf_ps1_ssu_g_root_filtered <- data.frame(ps1_ssu_g_root_filtered@sam_data)
ps1_ssu_g_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2493 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2493 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root_associated_soils <- subset_samples(ps1_ssu_g, type == "root associated soils")
sampledf_ps1_ssu_g_root_associated_soils <- data.frame(ps1_ssu_g_root_associated_soils@sam_data)
ps1_ssu_g_root_associated_soils
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_root_associated_soils) > 0, ps1_ssu_g_root_associated_soils)
sampledf_ps1_ssu_g_root_associated_soils_filtered <- data.frame(ps1_ssu_g_root_associated_soils_filtered@sam_data)
ps1_ssu_g_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 14910 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 14910 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_soil <- subset_samples(ps1_ssu_g, type == "soil")
sampledf_ps1_ssu_g_soil <- data.frame(ps1_ssu_g_soil@sam_data)
ps1_ssu_g_soil
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_g_soil) > 0, ps1_ssu_g_soil)
sampledf_ps1_ssu_g_soil_filtered <- data.frame(ps1_ssu_g_soil_filtered@sam_data)
ps1_ssu_g_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 11835 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 11835 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_leaf <- subset_samples(ps1_ssu_f, type == "leaf")
sampledf_ps1_ssu_f_leaf <- data.frame(ps1_ssu_f_leaf@sam_data)
ps1_ssu_f_leaf
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_leaf) > 0, ps1_ssu_f_leaf)
sampledf_ps1_ssu_f_leaf_filtered <- data.frame(ps1_ssu_f_leaf_filtered@sam_data)
ps1_ssu_f_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 307 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 307 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_stem <- subset_samples(ps1_ssu_f, type == "stem")
sampledf_ps1_ssu_f_stem <- data.frame(ps1_ssu_f_stem@sam_data)
ps1_ssu_f_stem
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_stem) > 0, ps1_ssu_f_stem)
sampledf_ps1_ssu_f_stem_filtered <- data.frame(ps1_ssu_f_stem_filtered@sam_data)
ps1_ssu_f_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 600 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 600 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root <- subset_samples(ps1_ssu_f, type == "root")
sampledf_ps1_ssu_f_root <- data.frame(ps1_ssu_f_root@sam_data)
ps1_ssu_f_root
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_root) > 0, ps1_ssu_f_root)
sampledf_ps1_ssu_f_root_filtered <- data.frame(ps1_ssu_f_root_filtered@sam_data)
ps1_ssu_f_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3101 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3101 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root_associated_soils <- subset_samples(ps1_ssu_f, type == "root associated soils")
sampledf_ps1_ssu_f_root_associated_soils <- data.frame(ps1_ssu_f_root_associated_soils@sam_data)
ps1_ssu_f_root_associated_soils
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_root_associated_soils) > 0, ps1_ssu_f_root_associated_soils)
sampledf_ps1_ssu_f_root_associated_soils_filtered <- data.frame(ps1_ssu_f_root_associated_soils_filtered@sam_data)
ps1_ssu_f_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 18361 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 18361 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_soil <- subset_samples(ps1_ssu_f, type == "soil")
sampledf_ps1_ssu_f_soil <- data.frame(ps1_ssu_f_soil@sam_data)
ps1_ssu_f_soil
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 46524 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 46524 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_ssu_f_soil) > 0, ps1_ssu_f_soil)
sampledf_ps1_ssu_f_soil_filtered <- data.frame(ps1_ssu_f_soil_filtered@sam_data)
ps1_ssu_f_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 12209 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 12209 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_leaf <- subset_samples(ps1_its_g, type == "leaf")
sampledf_ps1_its_g_leaf <- data.frame(ps1_its_g_leaf@sam_data)
ps1_its_g_leaf
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_leaf) > 0, ps1_its_g_leaf)
sampledf_ps1_its_g_leaf_filtered <- data.frame(ps1_its_g_leaf_filtered@sam_data)
ps1_its_g_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 84 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 84 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_stem <- subset_samples(ps1_its_g, type == "stem")
sampledf_ps1_its_g_stem <- data.frame(ps1_its_g_stem@sam_data)
ps1_its_g_stem
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_stem) > 0, ps1_its_g_stem)
sampledf_ps1_its_g_stem_filtered <- data.frame(ps1_its_g_stem_filtered@sam_data)
ps1_its_g_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 95 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 95 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root <- subset_samples(ps1_its_g, type == "root")
sampledf_ps1_its_g_root <- data.frame(ps1_its_g_root@sam_data)
ps1_its_g_root
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_root) > 0, ps1_its_g_root)
sampledf_ps1_its_g_root_filtered <- data.frame(ps1_its_g_root_filtered@sam_data)
ps1_its_g_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 180 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 180 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root_associated_soils <- subset_samples(ps1_its_g, type == "root associated soils")
sampledf_ps1_its_g_root_associated_soils <- data.frame(ps1_its_g_root_associated_soils@sam_data)
ps1_its_g_root_associated_soils
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_root_associated_soils) > 0, ps1_its_g_root_associated_soils)
sampledf_ps1_its_g_root_associated_soils_filtered <- data.frame(ps1_its_g_root_associated_soils_filtered@sam_data)
ps1_its_g_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1121 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1121 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_soil <- subset_samples(ps1_its_g, type == "soil")
sampledf_ps1_its_g_soil <- data.frame(ps1_its_g_soil@sam_data)
ps1_its_g_soil
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_g_soil) > 0, ps1_its_g_soil)
sampledf_ps1_its_g_soil_filtered <- data.frame(ps1_its_g_soil_filtered@sam_data)
ps1_its_g_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 952 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 952 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_leaf <- subset_samples(ps1_its_f, type == "leaf")
sampledf_ps1_its_f_leaf <- data.frame(ps1_its_f_leaf@sam_data)
ps1_its_f_leaf
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 21 samples ]
    ## sample_data() Sample Data:       [ 21 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_leaf_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_leaf) > 0, ps1_its_f_leaf)
sampledf_ps1_its_f_leaf_filtered <- data.frame(ps1_its_f_leaf_filtered@sam_data)
ps1_its_f_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 81 taxa and 21 samples ]
    ## sample_data() Sample Data:       [ 21 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 81 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_stem <- subset_samples(ps1_its_f, type == "stem")
sampledf_ps1_its_f_stem <- data.frame(ps1_its_f_stem@sam_data)
ps1_its_f_stem
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_stem_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_stem) > 0, ps1_its_f_stem)
sampledf_ps1_its_f_stem_filtered <- data.frame(ps1_its_f_stem_filtered@sam_data)
ps1_its_f_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 190 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 190 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root <- subset_samples(ps1_its_f, type == "root")
sampledf_ps1_its_f_root <- data.frame(ps1_its_f_root@sam_data)
ps1_its_f_root
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_root) > 0, ps1_its_f_root)
sampledf_ps1_its_f_root_filtered <- data.frame(ps1_its_f_root_filtered@sam_data)
ps1_its_f_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 467 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 467 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root_associated_soils <- subset_samples(ps1_its_f, type == "root associated soils")
sampledf_ps1_its_f_root_associated_soils <- data.frame(ps1_its_f_root_associated_soils@sam_data)
ps1_its_f_root_associated_soils
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root_associated_soils_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_root_associated_soils) > 0, ps1_its_f_root_associated_soils)
sampledf_ps1_its_f_root_associated_soils_filtered <- data.frame(ps1_its_f_root_associated_soils_filtered@sam_data)
ps1_its_f_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1322 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1322 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_soil <- subset_samples(ps1_its_f, type == "soil")
sampledf_ps1_its_f_soil <- data.frame(ps1_its_f_soil@sam_data)
ps1_its_f_soil
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3343 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3343 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_soil_filtered <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps1_its_f_soil) > 0, ps1_its_f_soil)
sampledf_ps1_its_f_soil_filtered <- data.frame(ps1_its_f_soil_filtered@sam_data)
ps1_its_f_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 556 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 556 taxa by 7 taxonomic ranks ]

## Q1. Were bacterial and fungal microbiomes significantly different according to the genotype variation?

# dataset

``` r
ps1_ssu_g_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 254 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 254 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 377 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 377 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 2493 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 2493 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 14910 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 14910 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_g_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 11835 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 11835 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 307 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 307 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 600 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 600 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 3101 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 3101 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 18361 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 18361 taxa by 7 taxonomic ranks ]

``` r
ps1_ssu_f_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 12209 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 12209 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 84 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 84 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 95 taxa and 23 samples ]
    ## sample_data() Sample Data:       [ 23 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 95 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 180 taxa and 22 samples ]
    ## sample_data() Sample Data:       [ 22 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 180 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1121 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1121 taxa by 7 taxonomic ranks ]

``` r
ps1_its_g_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 952 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 952 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_leaf_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 81 taxa and 21 samples ]
    ## sample_data() Sample Data:       [ 21 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 81 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_stem_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 190 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 190 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 467 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 467 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_root_associated_soils_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1322 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1322 taxa by 7 taxonomic ranks ]

``` r
ps1_its_f_soil_filtered
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 556 taxa and 24 samples ]
    ## sample_data() Sample Data:       [ 24 samples by 8 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 556 taxa by 7 taxonomic ranks ]

``` r
analysis_f <- function (dataset_x) {
# 1. Remove data points with missing metadata
phy_dataset_x <- dataset_x %>%
  subset_samples(
      !is.na(genotype) &
      !is.na(type) &
      !is.na(loc.time) &
      !is.na(id))

# 2. Differences depending on input variables
bray_dataset_x <- phyloseq::distance(physeq = phy_dataset_x, method = "bray")
cap_ord_dataset_x <- ordinate(physeq = phy_dataset_x, method = "CAP", distance = bray_dataset_x, formula = ~ genotype)

# 3. ANOVA
obj_phy_dataset_x <- t(phy_dataset_x)
anova(cap_ord_dataset_x, by="terms")
vif.cca(cap_ord_dataset_x)

# 4. CAP plot
cap_plot_dataset_x <- plot_ordination(physeq = phy_dataset_x, 
                                ordination = cap_ord_dataset_x,
                                color = "genotype", 
                                axes = c(1,2)) + xlim(-2.3, 2.2) + ylim(-2.2, 2.5) + 
                                geom_point(aes(colour = genotype), 
                                           alpha = 0.8, size = 3) +
                                geom_point(size = 3) + 
                                annotate("text", x = -0.8, y = 2.5, label = bquote("R[PERMANOVA]^2 == 0.015"), 
                                         parse = TRUE) +
                                annotate("text", x = 0.9, y = 2.48, label = bquote("p[PERMANOVA] == 0.018"), 
                                         parse = TRUE) +
                                scale_color_manual(
                                  values = c("#00AFBB","#FC4E07", "#E7B800",
                                             "#36BB07", "#7D26CD",
                                             "#757575")) + theme_bw()
cap_plot_dataset_x_2 <- cap_plot_dataset_x
cap_plot_dataset_x_2$data$genotype<-factor(cap_plot_dataset_x$data$genotype, levels = c("WT", "L1", "L2", "L3", "L4", "L5"))
print(cap_plot_dataset_x_2)

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

# 8. Betadiversity / dispersion / homogeneity between groups
bet_factor_dataset_x <-betadisper(dist_phy_dataset_x, sampledf_dist_phy_dataset_x$genotype)
betdisper_result_dataset_x <- permutest(bet_factor_dataset_x, pairwise = TRUE, permutations = 999, group=bet_factor_dataset_x$group)
print(betdisper_result_dataset_x)

# 9. Alpha diversity plot
diversity_dataset_dataset_x <- prune_taxa(taxa_sums(dataset_x) > 0, dataset_x)
diversity_plot_dataset_x <- plot_richness(diversity_dataset_dataset_x, x="genotype",measures=c("Chao1", "Shannon"), color = "genotype") + theme_bw()
diversity_plot_dataset_x + geom_boxplot(data = diversity_plot_dataset_x$data, aes(x = genotype, y = value, color = NULL), alpha = 0.1)
#write.csv(diversity_plot_dataset_x$plot_env$DF, file="alpha_diversity_ssu.csv")
diversity_plot_dataset_x_data <- diversity_plot_dataset_x$plot_env$DF
ordered_diversity_plot_dataset_x_data <- reorder_levels(diversity_plot_dataset_x_data, genotype, c("WT", "L1", "L2", "L3", "L4", "L5"))

res.kruskal_total_shannon_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  kruskal_test(Shannon ~ genotype)
print(res.kruskal_total_shannon_dataset_x)

pwc_total_shannon_dataset_x <- ordered_diversity_plot_dataset_x_data %>% 
  dunn_test(Shannon ~ genotype, p.adjust.method = "bonferroni") 
print(pwc_total_shannon_dataset_x)

pwc_total_shannon_dataset_x <- pwc_total_shannon_dataset_x %>% add_xy_position()
shannon_plot <- ggboxplot(ordered_diversity_plot_dataset_x_data, x = "genotype", y = "Shannon", color = "genotype", 
          palette = c("#00AFBB","#FC4E07", "#E7B800", "#36BB07", "#7D26CD", "#757575"),
          order = c("WT", "L1", "L2", "L3", "L4", "L5"),
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
          order = c("WT", "L1", "L2", "L3", "L4", "L5"),
          ylab = "chao1 index", xlab = "genotype") +
          stat_pvalue_manual(pwc_total_chao1_dataset_x, y.position = 2300,
          step.increase = 0.1, hide.ns = TRUE) + 
          labs(subtitle = get_test_label(res.kruskal_total_chao1_dataset_x, detailed = TRUE),
          caption = get_pwc_label(pwc_total_chao1_dataset_x)) + coord_cartesian(ylim = c(0, 4500))
print(chao_plot)

ps_normalization_dataset_x = transform_sample_counts(dataset_x, function(x) x/sum(x))
phyla_plot <- phyloseq::plot_bar(ps_normalization_dataset_x, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ genotype, scales = "free") + theme(panel.background = element_blank() +
  theme(legend.position = "bottom"), axis.text.x=element_blank(), axis.ticks.x=element_blank())
print(phyla_plot)

#venn diagramming
wt_sub_ssu <- subset_samples(dataset_x, genotype == "WT") %>% filter_taxa(function(x) sum(x) > 0, T)
l1_sub_ssu <- subset_samples(dataset_x, genotype == "L1") %>% filter_taxa(function(x) sum(x) > 0, T)
l2_sub_ssu <- subset_samples(dataset_x, genotype == "L2") %>% filter_taxa(function(x) sum(x) > 0, T)
l3_sub_ssu <- subset_samples(dataset_x, genotype == "L3") %>% filter_taxa(function(x) sum(x) > 0, T)
l4_sub_ssu <- subset_samples(dataset_x, genotype == "L4") %>% filter_taxa(function(x) sum(x) > 0, T)
l5_sub_ssu <- subset_samples(dataset_x, genotype == "L5") %>% filter_taxa(function(x) sum(x) > 0, T)
wt_asvs_ssu <- taxa_names(wt_sub_ssu)
l1_asvs_ssu <- taxa_names(l1_sub_ssu)
l2_asvs_ssu <- taxa_names(l2_sub_ssu)
l3_asvs_ssu <- taxa_names(l3_sub_ssu)
l4_asvs_ssu <- taxa_names(l4_sub_ssu)
l5_asvs_ssu <- taxa_names(l5_sub_ssu)
genotype_venn_ssu <- venn(list("WT" = wt_asvs_ssu, "L1" = l1_asvs_ssu, "L2" = l2_asvs_ssu, "L3" = l3_asvs_ssu, "L4" = l4_asvs_ssu, "L5"= l5_asvs_ssu))
#write.csv(l1_asvs_ssu,file="l1_asvs_ssu.csv")
#write.csv(l2_asvs_ssu,file="l2_asvs_ssu.csv")
#write.csv(l3_asvs_ssu,file="l3_asvs_ssu.csv")
#write.csv(l4_asvs_ssu,file="l4_asvs_ssu.csv")
#write.csv(l5_asvs_ssu,file="l5_asvs_ssu.csv")
#write.csv(wt_asvs_ssu,file="wt_asvs_ssu.csv")
}
```

``` r
analysis_f(ps1_ssu_g_leaf_filtered)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ##                R2 Pr(>F)  
    ## genotype  0.39651  0.041 *
    ## Residuals 0.60349         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df   SumsOfSqs    F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.095315868  3.6497420 0.37822172   0.115  0.2875000    
    ## 2  L1 vs L3  1 0.074089124  1.7736549 0.22816229   0.169  0.2895000    
    ## 3  L1 vs L4  1 0.059194073  1.0142823 0.14460243   0.403  0.5495455    
    ## 4  L1 vs L5  1 0.076725734  3.0225880 0.37675972   0.144  0.2895000    
    ## 5  L1 vs WT  1 0.003204112  0.1717821 0.02783347   0.859  0.9110000    
    ## 6  L2 vs L3  1 0.022107672  0.5192182 0.07964424   0.546  0.6825000    
    ## 7  L2 vs L4  1 0.017093881  0.2889075 0.04593922   0.911  0.9110000    
    ## 8  L2 vs L5  1 0.283415991 10.7549132 0.68263868   0.033  0.1650000    
    ## 9  L2 vs WT  1 0.102787831  5.2822918 0.46819316   0.030  0.1650000    
    ## 10 L3 vs L4  1 0.024957043  0.3335453 0.05266329   0.640  0.7384615    
    ## 11 L3 vs L5  1 0.220728786  4.8899032 0.49443388   0.026  0.1650000    
    ## 12 L3 vs WT  1 0.082417884  2.3470734 0.28118519   0.176  0.2895000    
    ## 13 L4 vs L5  1 0.209542772  3.2214573 0.39183531   0.096  0.2875000    
    ## 14 L4 vs WT  1 0.067871806  1.3127066 0.17951036   0.193  0.2895000    
    ## 15 L5 vs WT  1 0.077339112  4.4458202 0.47066534   0.093  0.2875000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     5 0.038102 0.0076204 2.8963    999  0.047 *
    ## Residuals 17 0.044729 0.0026311                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.3870000 0.4460000 0.2790000 0.0680000 0.192
    ## L2 0.3859514           0.8770000 0.5890000 0.1530000 0.023
    ## L3 0.4592197 0.8712403           0.5390000 0.1290000 0.039
    ## L4 0.2731515 0.6075788 0.5327854           0.5860000 0.049
    ## L5 0.0749717 0.1641981 0.1306237 0.6020433           0.006
    ## WT 0.2046310 0.0248829 0.0316000 0.0463204 0.0032668

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df       p method        
    ## * <chr>   <int>     <dbl> <int>   <dbl> <chr>         
    ## 1 Shannon    23      15.5     5 0.00857 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic        p  p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>    <dbl>  <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     1.20  0.231    1      ns          
    ##  2 Shannon WT     L2         4     4     1.62  0.106    1      ns          
    ##  3 Shannon WT     L3         4     4     2.81  0.00488  0.0732 ns          
    ##  4 Shannon WT     L4         4     4     0.938 0.348    1      ns          
    ##  5 Shannon WT     L5         4     3     3.36  0.000773 0.0116 *           
    ##  6 Shannon L1     L2         4     4     0.417 0.677    1      ns          
    ##  7 Shannon L1     L3         4     4     1.62  0.106    1      ns          
    ##  8 Shannon L1     L4         4     4    -0.261 0.794    1      ns          
    ##  9 Shannon L1     L5         4     3     2.25  0.0243   0.365  ns          
    ## 10 Shannon L2     L3         4     4     1.20  0.231    1      ns          
    ## 11 Shannon L2     L4         4     4    -0.678 0.498    1      ns          
    ## 12 Shannon L2     L5         4     3     1.87  0.0620   0.930  ns          
    ## 13 Shannon L3     L4         4     4    -1.88  0.0606   0.909  ns          
    ## 14 Shannon L3     L5         4     3     0.756 0.450    1      ns          
    ## 15 Shannon L4     L5         4     3     2.49  0.0126   0.190  ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df      p method        
    ## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Chao1    23      9.46     5 0.0921 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     1.02  0.309   1     ns          
    ##  2 Chao1 WT     L2         4     4     0.339 0.734   1     ns          
    ##  3 Chao1 WT     L3         4     4     1.20  0.230   1     ns          
    ##  4 Chao1 WT     L4         4     4    -0.313 0.754   1     ns          
    ##  5 Chao1 WT     L5         4     3     2.42  0.0157  0.235 ns          
    ##  6 Chao1 L1     L2         4     4    -0.679 0.497   1     ns          
    ##  7 Chao1 L1     L3         4     4     0.183 0.855   1     ns          
    ##  8 Chao1 L1     L4         4     4    -1.33  0.183   1     ns          
    ##  9 Chao1 L1     L5         4     3     1.47  0.140   1     ns          
    ## 10 Chao1 L2     L3         4     4     0.861 0.389   1     ns          
    ## 11 Chao1 L2     L4         4     4    -0.653 0.514   1     ns          
    ## 12 Chao1 L2     L5         4     3     2.10  0.0355  0.533 ns          
    ## 13 Chao1 L3     L4         4     4    -1.51  0.130   1     ns          
    ## 14 Chao1 L3     L5         4     3     1.31  0.192   1     ns          
    ## 15 Chao1 L4     L5         4     3     2.71  0.00680 0.102 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-14-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-14-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-14-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_g_stem_filtered)
```

    ##                R2 Pr(>F)
    ## genotype  0.22028  0.453
    ## Residuals 0.77972       
    ## Total     1.00000

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df   SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.026212244 0.5479463 0.08368217   0.777     0.8510    
    ## 2  L1 vs L3  1 0.043183926 0.6488470 0.09758790   0.666     0.8510    
    ## 3  L1 vs L4  1 0.040661354 0.9932641 0.14203154   0.477     0.7215    
    ## 4  L1 vs L5  1 0.023211645 0.5046013 0.09166900   0.804     0.8510    
    ## 5  L1 vs WT  1 0.039501681 0.8149921 0.11958812   0.404     0.7215    
    ## 6  L2 vs L3  1 0.053245174 1.1536700 0.16126967   0.304     0.7215    
    ## 7  L2 vs L4  1 0.031369065 1.5275830 0.20293141   0.366     0.7215    
    ## 8  L2 vs L5  1 0.007007064 0.3256439 0.06114640   0.851     0.8510    
    ## 9  L2 vs WT  1 0.021799692 0.7767078 0.11461433   0.470     0.7215    
    ## 10 L3 vs L4  1 0.010575825 0.2694291 0.04297506   0.794     0.8510    
    ## 11 L3 vs L5  1 0.028413984 0.6460850 0.11443061   0.481     0.7215    
    ## 12 L3 vs WT  1 0.093703128 2.0028709 0.25026905   0.201     0.7215    
    ## 13 L4 vs L5  1 0.010579615 0.7992221 0.13781539   0.375     0.7215    
    ## 14 L4 vs WT  1 0.072239697 3.4129031 0.36257710   0.077     0.7215    
    ## 15 L5 vs WT  1 0.024955339 1.1203084 0.18304771   0.364     0.7215    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.079012 0.015803 1.5172    999  0.194
    ## Residuals 17 0.177064 0.010415                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.892000 0.129000 0.124000 0.185000 0.506
    ## L2 0.815117          0.183000 0.178000 0.251000 0.659
    ## L3 0.167403 0.207384          0.892000 0.765000 0.013
    ## L4 0.175352 0.218336 0.854182          0.558000 0.013
    ## L5 0.209724 0.243653 0.647715 0.456223          0.012
    ## WT 0.396183 0.527658 0.034946 0.029843 0.017824

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df      p method        
    ## * <chr>   <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Shannon    23      15.0     5 0.0102 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4   -0.0521 0.958   1     ns          
    ##  2 Shannon WT     L2         4     4    0.469  0.639   1     ns          
    ##  3 Shannon WT     L3         4     4   -2.14   0.0326  0.489 ns          
    ##  4 Shannon WT     L4         4     4   -1.93   0.0538  0.806 ns          
    ##  5 Shannon WT     L5         4     3   -2.16   0.0311  0.467 ns          
    ##  6 Shannon L1     L2         4     4    0.521  0.602   1     ns          
    ##  7 Shannon L1     L3         4     4   -2.09   0.0371  0.556 ns          
    ##  8 Shannon L1     L4         4     4   -1.88   0.0606  0.909 ns          
    ##  9 Shannon L1     L5         4     3   -2.11   0.0351  0.526 ns          
    ## 10 Shannon L2     L3         4     4   -2.61   0.00915 0.137 ns          
    ## 11 Shannon L2     L4         4     4   -2.40   0.0165  0.247 ns          
    ## 12 Shannon L2     L5         4     3   -2.59   0.00960 0.144 ns          
    ## 13 Shannon L3     L4         4     4    0.209  0.835   1     ns          
    ## 14 Shannon L3     L5         4     3   -0.177  0.860   1     ns          
    ## 15 Shannon L4     L5         4     3   -0.370  0.711   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    23      2.84     5 0.725 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    0.235  0.814     1 ns          
    ##  2 Chao1 WT     L2         4     4    0.522  0.601     1 ns          
    ##  3 Chao1 WT     L3         4     4    0.444  0.657     1 ns          
    ##  4 Chao1 WT     L4         4     4   -0.0783 0.938     1 ns          
    ##  5 Chao1 WT     L5         4     3   -1.02   0.310     1 ns          
    ##  6 Chao1 L1     L2         4     4    0.287  0.774     1 ns          
    ##  7 Chao1 L1     L3         4     4    0.209  0.835     1 ns          
    ##  8 Chao1 L1     L4         4     4   -0.313  0.754     1 ns          
    ##  9 Chao1 L1     L5         4     3   -1.23   0.218     1 ns          
    ## 10 Chao1 L2     L3         4     4   -0.0783 0.938     1 ns          
    ## 11 Chao1 L2     L4         4     4   -0.601  0.548     1 ns          
    ## 12 Chao1 L2     L5         4     3   -1.50   0.134     1 ns          
    ## 13 Chao1 L3     L4         4     4   -0.522  0.601     1 ns          
    ## 14 Chao1 L3     L5         4     3   -1.43   0.154     1 ns          
    ## 15 Chao1 L4     L5         4     3   -0.943  0.346     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_g_root_filtered)
```

    ##                R2 Pr(>F)    
    ## genotype  0.37689  0.001 ***
    ## Residuals 0.62311           
    ## Total     1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.3015539 1.6581941 0.2165255   0.038 0.06150000    
    ## 2  L1 vs L3  1 0.5939397 2.3212481 0.2789543   0.029 0.06150000    
    ## 3  L1 vs L4  1 0.3282434 1.6789089 0.2186390   0.062 0.08454545    
    ## 4  L1 vs L5  1 0.9716502 4.2382354 0.4139615   0.024 0.06150000    
    ## 5  L1 vs WT  1 0.2032540 1.0792018 0.1524468   0.289 0.30964286    
    ## 6  L2 vs L3  1 0.3796129 1.6116145 0.2117310   0.027 0.06150000    
    ## 7  L2 vs L4  1 0.1667751 0.9519817 0.1369367   0.560 0.56000000    
    ## 8  L2 vs L5  1 0.7378677 3.5315568 0.3705121   0.039 0.06150000    
    ## 9  L2 vs WT  1 0.2062848 1.2277784 0.1698694   0.221 0.25500000    
    ## 10 L3 vs L4  1 0.4064005 1.6308121 0.2137141   0.025 0.06150000    
    ## 11 L3 vs L5  1 0.4655387 1.6453063 0.2152048   0.036 0.06150000    
    ## 12 L3 vs WT  1 0.5153129 2.1291396 0.2619145   0.030 0.06150000    
    ## 13 L4 vs L5  1 0.6658030 2.9911810 0.3326794   0.030 0.06150000    
    ## 14 L4 vs WT  1 0.2528473 1.3918113 0.1882910   0.114 0.14250000    
    ## 15 L5 vs WT  1 0.8431848 3.9142153 0.3948084   0.041 0.06150000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
    ## Groups     5 0.048286 0.0096571 4.3886    999  0.009 **
    ## Residuals 18 0.039609 0.0022005                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.5950000 0.1810000 0.9750000 0.0730000 0.603
    ## L2 0.6232538           0.0320000 0.4230000 0.0070000 0.987
    ## L3 0.1566753 0.0197348           0.0570000 0.6740000 0.026
    ## L4 0.9789926 0.4398846 0.0437843           0.0090000 0.387
    ## L5 0.0818091 0.0029080 0.6992491 0.0060327           0.004
    ## WT 0.6071786 0.9816215 0.0163457 0.4021270 0.0019041

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df       p method        
    ## * <chr>   <int>     <dbl> <int>   <dbl> <chr>         
    ## 1 Shannon    24      16.1     5 0.00651 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -0.45 0.653   1      ns          
    ##  2 Shannon WT     L2         4     4     -0.9  0.368   1      ns          
    ##  3 Shannon WT     L3         4     4     -2.5  0.0124  0.186  ns          
    ##  4 Shannon WT     L4         4     4     -0.5  0.617   1      ns          
    ##  5 Shannon WT     L5         4     4     -3.15 0.00163 0.0245 *           
    ##  6 Shannon L1     L2         4     4     -0.45 0.653   1      ns          
    ##  7 Shannon L1     L3         4     4     -2.05 0.0404  0.605  ns          
    ##  8 Shannon L1     L4         4     4     -0.05 0.960   1      ns          
    ##  9 Shannon L1     L5         4     4     -2.7  0.00693 0.104  ns          
    ## 10 Shannon L2     L3         4     4     -1.6  0.110   1      ns          
    ## 11 Shannon L2     L4         4     4      0.4  0.689   1      ns          
    ## 12 Shannon L2     L5         4     4     -2.25 0.0244  0.367  ns          
    ## 13 Shannon L3     L4         4     4      2    0.0455  0.683  ns          
    ## 14 Shannon L3     L5         4     4     -0.65 0.516   1      ns          
    ## 15 Shannon L4     L5         4     4     -2.65 0.00805 0.121  ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df       p method        
    ## * <chr> <int>     <dbl> <int>   <dbl> <chr>         
    ## 1 Chao1    24      16.4     5 0.00584 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic        p   p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>    <dbl>   <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    -0.550 0.582    1       ns          
    ##  2 Chao1 WT     L2         4     4    -0.975 0.329    1       ns          
    ##  3 Chao1 WT     L3         4     4    -2.43  0.0153   0.229   ns          
    ##  4 Chao1 WT     L4         4     4    -1.43  0.154    1       ns          
    ##  5 Chao1 WT     L5         4     4    -3.48  0.000509 0.00764 **          
    ##  6 Chao1 L1     L2         4     4    -0.425 0.671    1       ns          
    ##  7 Chao1 L1     L3         4     4    -1.88  0.0607   0.911   ns          
    ##  8 Chao1 L1     L4         4     4    -0.875 0.381    1       ns          
    ##  9 Chao1 L1     L5         4     4    -2.93  0.00344  0.0516  ns          
    ## 10 Chao1 L2     L3         4     4    -1.45  0.147    1       ns          
    ## 11 Chao1 L2     L4         4     4    -0.450 0.653    1       ns          
    ## 12 Chao1 L2     L5         4     4    -2.50  0.0124   0.186   ns          
    ## 13 Chao1 L3     L4         4     4     1.00  0.317    1       ns          
    ## 14 Chao1 L3     L5         4     4    -1.05  0.294    1       ns          
    ## 15 Chao1 L4     L5         4     4    -2.05  0.0403   0.605   ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-16-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-16-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-16-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_g_root_associated_soils_filtered)
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

    ##                R2 Pr(>F)   
    ## genotype  0.32323  0.007 **
    ## Residuals 0.67677          
    ## Total     1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.1528234 1.0297117 0.1464799   0.423  0.4880769    
    ## 2  L1 vs L3  1 0.5121188 3.1674533 0.3455107   0.029  0.1125000    
    ## 3  L1 vs L4  1 0.1528015 0.8033054 0.1180758   0.730  0.7300000    
    ## 4  L1 vs L5  1 0.2410214 1.1431091 0.1600296   0.289  0.4131818    
    ## 5  L1 vs WT  1 0.1712713 1.0592952 0.1500568   0.303  0.4131818    
    ## 6  L2 vs L3  1 0.4077496 3.2358703 0.3503590   0.027  0.1125000    
    ## 7  L2 vs L4  1 0.2075929 1.3432642 0.1829247   0.166  0.2766667    
    ## 8  L2 vs L5  1 0.3043705 1.7375231 0.2245581   0.071  0.1775000    
    ## 9  L2 vs WT  1 0.1211996 0.9618106 0.1381552   0.398  0.4880769    
    ## 10 L3 vs L4  1 0.5143148 3.0648384 0.3381018   0.025  0.1125000    
    ## 11 L3 vs L5  1 0.3495461 1.8549195 0.2361475   0.101  0.1912500    
    ## 12 L3 vs WT  1 0.4166981 2.9918086 0.3327260   0.030  0.1125000    
    ## 13 L4 vs L5  1 0.1922447 0.8860140 0.1286686   0.482  0.5164286    
    ## 14 L4 vs WT  1 0.2509121 1.4951794 0.1994855   0.102  0.1912500    
    ## 15 L5 vs WT  1 0.3438585 1.8247119 0.2331986   0.041  0.1230000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     5 0.036710 0.0073421 2.7233    999  0.039 *
    ## Residuals 18 0.048529 0.0026961                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.0810000 0.0790000 0.7200000 0.5950000 0.214
    ## L2 0.0917669           0.7380000 0.0060000 0.0870000 0.344
    ## L3 0.0750967 0.7294999           0.0010000 0.0980000 0.314
    ## L4 0.6938777 0.0131123 0.0015548           0.7530000 0.020
    ## L5 0.5944983 0.0996290 0.1005386 0.7217120           0.183
    ## WT 0.2222610 0.3490084 0.3112410 0.0288311 0.1852616      
    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      4.57     5 0.471 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      0.55 0.582  1     ns          
    ##  2 Shannon WT     L2         4     4      0.8  0.424  1     ns          
    ##  3 Shannon WT     L3         4     4      1.4  0.162  1     ns          
    ##  4 Shannon WT     L4         4     4      1.4  0.162  1     ns          
    ##  5 Shannon WT     L5         4     4      1.85 0.0643 0.965 ns          
    ##  6 Shannon L1     L2         4     4      0.25 0.803  1     ns          
    ##  7 Shannon L1     L3         4     4      0.85 0.395  1     ns          
    ##  8 Shannon L1     L4         4     4      0.85 0.395  1     ns          
    ##  9 Shannon L1     L5         4     4      1.3  0.194  1     ns          
    ## 10 Shannon L2     L3         4     4      0.6  0.549  1     ns          
    ## 11 Shannon L2     L4         4     4      0.6  0.549  1     ns          
    ## 12 Shannon L2     L5         4     4      1.05 0.294  1     ns          
    ## 13 Shannon L3     L4         4     4      0    1      1     ns          
    ## 14 Shannon L3     L5         4     4      0.45 0.653  1     ns          
    ## 15 Shannon L4     L5         4     4      0.45 0.653  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      6.34     5 0.275 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4      1.4  0.162  1     ns          
    ##  2 Chao1 WT     L2         4     4      0.85 0.395  1     ns          
    ##  3 Chao1 WT     L3         4     4      0.65 0.516  1     ns          
    ##  4 Chao1 WT     L4         4     4      2    0.0455 0.683 ns          
    ##  5 Chao1 WT     L5         4     4      2    0.0455 0.683 ns          
    ##  6 Chao1 L1     L2         4     4     -0.55 0.582  1     ns          
    ##  7 Chao1 L1     L3         4     4     -0.75 0.453  1     ns          
    ##  8 Chao1 L1     L4         4     4      0.6  0.549  1     ns          
    ##  9 Chao1 L1     L5         4     4      0.6  0.549  1     ns          
    ## 10 Chao1 L2     L3         4     4     -0.2  0.841  1     ns          
    ## 11 Chao1 L2     L4         4     4      1.15 0.250  1     ns          
    ## 12 Chao1 L2     L5         4     4      1.15 0.250  1     ns          
    ## 13 Chao1 L3     L4         4     4      1.35 0.177  1     ns          
    ## 14 Chao1 L3     L5         4     4      1.35 0.177  1     ns          
    ## 15 Chao1 L4     L5         4     4      0    1      1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-17-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-17-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_g_soil_filtered)
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

    ##                R2 Pr(>F)  
    ## genotype  0.30101  0.025 *
    ## Residuals 0.69899         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.12124125 1.0294434 0.14644736   0.335  0.5025000    
    ## 2  L1 vs L3  1 0.19436962 1.0063647 0.14363579   0.454  0.5238462    
    ## 3  L1 vs L4  1 0.09446161 0.8159360 0.11971005   0.748  0.8014286    
    ## 4  L1 vs L5  1 0.18942634 1.1427050 0.15998211   0.335  0.5025000    
    ## 5  L1 vs WT  1 0.35034408 2.6586607 0.30705219   0.021  0.1200000    
    ## 6  L2 vs L3  1 0.27162799 1.5044733 0.20047686   0.317  0.5025000    
    ## 7  L2 vs L4  1 0.10720312 1.0390172 0.14760828   0.416  0.5200000    
    ## 8  L2 vs L5  1 0.25532782 1.6668843 0.21741352   0.173  0.4325000    
    ## 9  L2 vs WT  1 0.22881169 1.9198636 0.24241120   0.024  0.1200000    
    ## 10 L3 vs L4  1 0.20306203 1.1373210 0.15934844   0.384  0.5200000    
    ## 11 L3 vs L5  1 0.09034860 0.3953235 0.06181447   0.946  0.9460000    
    ## 12 L3 vs WT  1 0.48521739 2.4940759 0.29362533   0.030  0.1200000    
    ## 13 L4 vs L5  1 0.18373936 1.2154168 0.16844721   0.297  0.5025000    
    ## 14 L4 vs WT  1 0.33388930 2.8494080 0.32198854   0.032  0.1200000    
    ## 15 L5 vs WT  1 0.48369035 2.8932700 0.32533253   0.041  0.1230000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)   
    ## Groups     5 0.043385 0.0086769 4.7947    999  0.004 **
    ## Residuals 18 0.032574 0.0018097                        
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.370000 0.063000 0.208000 0.052000 0.470
    ## L2 0.350744          0.030000 0.663000 0.015000 0.125
    ## L3 0.068214 0.032898          0.019000 0.445000 0.112
    ## L4 0.216522 0.618592 0.026387          0.011000 0.067
    ## L5 0.065651 0.014912 0.410014 0.011499          0.178
    ## WT 0.482087 0.122454 0.123808 0.080190 0.188142      
    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      8.03     5 0.155 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      2.65 0.00805 0.121 ns          
    ##  2 Shannon WT     L2         4     4      1.9  0.0574  0.861 ns          
    ##  3 Shannon WT     L3         4     4      1.65 0.0989  1     ns          
    ##  4 Shannon WT     L4         4     4      2    0.0455  0.683 ns          
    ##  5 Shannon WT     L5         4     4      2    0.0455  0.683 ns          
    ##  6 Shannon L1     L2         4     4     -0.75 0.453   1     ns          
    ##  7 Shannon L1     L3         4     4     -1    0.317   1     ns          
    ##  8 Shannon L1     L4         4     4     -0.65 0.516   1     ns          
    ##  9 Shannon L1     L5         4     4     -0.65 0.516   1     ns          
    ## 10 Shannon L2     L3         4     4     -0.25 0.803   1     ns          
    ## 11 Shannon L2     L4         4     4      0.1  0.920   1     ns          
    ## 12 Shannon L2     L5         4     4      0.1  0.920   1     ns          
    ## 13 Shannon L3     L4         4     4      0.35 0.726   1     ns          
    ## 14 Shannon L3     L5         4     4      0.35 0.726   1     ns          
    ## 15 Shannon L4     L5         4     4      0    1       1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      12.8     5 0.025 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     2.50  0.0124  0.186  ns          
    ##  2 Chao1 WT     L2         4     4     1.50  0.134   1      ns          
    ##  3 Chao1 WT     L3         4     4     2.80  0.00510 0.0765 ns          
    ##  4 Chao1 WT     L4         4     4     2.08  0.0379  0.569  ns          
    ##  5 Chao1 WT     L5         4     4     3.13  0.00177 0.0266 *           
    ##  6 Chao1 L1     L2         4     4    -1.00  0.317   1      ns          
    ##  7 Chao1 L1     L3         4     4     0.300 0.764   1      ns          
    ##  8 Chao1 L1     L4         4     4    -0.425 0.671   1      ns          
    ##  9 Chao1 L1     L5         4     4     0.625 0.532   1      ns          
    ## 10 Chao1 L2     L3         4     4     1.30  0.194   1      ns          
    ## 11 Chao1 L2     L4         4     4     0.575 0.565   1      ns          
    ## 12 Chao1 L2     L5         4     4     1.63  0.104   1      ns          
    ## 13 Chao1 L3     L4         4     4    -0.725 0.468   1      ns          
    ## 14 Chao1 L3     L5         4     4     0.325 0.745   1      ns          
    ## 15 Chao1 L4     L5         4     4     1.05  0.294   1      ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-18-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-18-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-18-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_f_leaf_filtered)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ##                R2 Pr(>F)
    ## genotype  0.27573  0.103
    ## Residuals 0.72427       
    ## Total     1.00000       
    ##       pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.20476449 1.5463094 0.20490936   0.220  0.4125000    
    ## 2  L1 vs L3  1 0.21263886 1.5469491 0.20497675   0.210  0.4125000    
    ## 3  L1 vs L4  1 0.15929150 1.0178442 0.14503659   0.595  0.6865385    
    ## 4  L1 vs L5  1 0.18494203 1.2619525 0.17377592   0.272  0.4533333    
    ## 5  L1 vs WT  1 0.20403940 1.6043216 0.21097498   0.166  0.4125000    
    ## 6  L2 vs L3  1 0.01823303 0.2358257 0.03781788   1.000  1.0000000    
    ## 7  L2 vs L4  1 0.20003186 2.0759302 0.25705153   0.061  0.3937500    
    ## 8  L2 vs L5  1 0.08190320 0.9478321 0.13642127   0.550  0.6865385    
    ## 9  L2 vs WT  1 0.04858533 0.7247227 0.10776990   0.566  0.6865385    
    ## 10 L3 vs L4  1 0.19127356 1.8864535 0.23920175   0.103  0.3937500    
    ## 11 L3 vs L5  1 0.12145357 1.3281367 0.18123798   0.304  0.4560000    
    ## 12 L3 vs WT  1 0.04131214 0.5731795 0.08719973   0.835  0.8946429    
    ## 13 L4 vs L5  1 0.20997613 1.9004330 0.24054795   0.105  0.3937500    
    ## 14 L4 vs WT  1 0.13509803 1.4826812 0.19814838   0.140  0.4125000    
    ## 15 L5 vs WT  1 0.14493736 1.7855856 0.22934506   0.096  0.3937500    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.023947 0.0047893 1.2148    999  0.383
    ## Residuals 18 0.070966 0.0039425                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.312000 0.021000 0.285000 0.740000 0.822
    ## L2 0.300238          0.605000 0.685000 0.347000 0.393
    ## L3 0.015771 0.602618          0.109000 0.155000 0.029
    ## L4 0.282678 0.678191 0.110061          0.409000 0.383
    ## L5 0.716208 0.341731 0.157692 0.413689          0.637
    ## WT 0.809082 0.358959 0.021619 0.381402 0.644614

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df      p method        
    ## * <chr>   <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Shannon    24      11.2     5 0.0468 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4       1.1 0.271   1     ns          
    ##  2 Shannon WT     L2         4     4      -0.8 0.424   1     ns          
    ##  3 Shannon WT     L3         4     4      -0.7 0.484   1     ns          
    ##  4 Shannon WT     L4         4     4       1.8 0.0719  1     ns          
    ##  5 Shannon WT     L5         4     4       1   0.317   1     ns          
    ##  6 Shannon L1     L2         4     4      -1.9 0.0574  0.861 ns          
    ##  7 Shannon L1     L3         4     4      -1.8 0.0719  1     ns          
    ##  8 Shannon L1     L4         4     4       0.7 0.484   1     ns          
    ##  9 Shannon L1     L5         4     4      -0.1 0.920   1     ns          
    ## 10 Shannon L2     L3         4     4       0.1 0.920   1     ns          
    ## 11 Shannon L2     L4         4     4       2.6 0.00932 0.140 ns          
    ## 12 Shannon L2     L5         4     4       1.8 0.0719  1     ns          
    ## 13 Shannon L3     L4         4     4       2.5 0.0124  0.186 ns          
    ## 14 Shannon L3     L5         4     4       1.7 0.0891  1     ns          
    ## 15 Shannon L4     L5         4     4      -0.8 0.424   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      2.33     5 0.802 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    -0.452 0.651     1 ns          
    ##  2 Chao1 WT     L2         4     4    -1.15  0.248     1 ns          
    ##  3 Chao1 WT     L3         4     4    -0.954 0.340     1 ns          
    ##  4 Chao1 WT     L4         4     4    -1.31  0.192     1 ns          
    ##  5 Chao1 WT     L5         4     4    -0.803 0.422     1 ns          
    ##  6 Chao1 L1     L2         4     4    -0.703 0.482     1 ns          
    ##  7 Chao1 L1     L3         4     4    -0.502 0.616     1 ns          
    ##  8 Chao1 L1     L4         4     4    -0.853 0.393     1 ns          
    ##  9 Chao1 L1     L5         4     4    -0.351 0.725     1 ns          
    ## 10 Chao1 L2     L3         4     4     0.201 0.841     1 ns          
    ## 11 Chao1 L2     L4         4     4    -0.151 0.880     1 ns          
    ## 12 Chao1 L2     L5         4     4     0.351 0.725     1 ns          
    ## 13 Chao1 L3     L4         4     4    -0.351 0.725     1 ns          
    ## 14 Chao1 L3     L5         4     4     0.151 0.880     1 ns          
    ## 15 Chao1 L4     L5         4     4     0.502 0.616     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-19-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_f_stem_filtered)
```

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ##                R2 Pr(>F)
    ## genotype  0.26801  0.113
    ## Residuals 0.73199       
    ## Total     1.00000       
    ##       pairs Df   SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.016977642 0.5251298 0.08047807   0.688    0.97300    
    ## 2  L1 vs L3  1 0.232995411 2.6635871 0.30744622   0.082    0.41000    
    ## 3  L1 vs L4  1 0.089978671 0.6773000 0.10143321   0.913    0.97300    
    ## 4  L1 vs L5  1 0.114992443 0.7379904 0.10952678   0.972    0.97300    
    ## 5  L1 vs WT  1 0.007944038 0.2263227 0.03634933   0.973    0.97300    
    ## 6  L2 vs L3  1 0.317420853 4.7596177 0.44235937   0.041    0.41000    
    ## 7  L2 vs L4  1 0.102710543 0.9165251 0.13251236   0.724    0.97300    
    ## 8  L2 vs L5  1 0.130366875 0.9654341 0.13860358   0.647    0.97300    
    ## 9  L2 vs WT  1 0.013542377 0.9459207 0.13618363   0.428    0.97300    
    ## 10 L3 vs L4  1 0.311935870 1.8655436 0.23717923   0.111    0.41625    
    ## 11 L3 vs L5  1 0.272508097 1.4329078 0.19277890   0.205    0.61500    
    ## 12 L3 vs WT  1 0.255345277 3.6761201 0.37991675   0.056    0.41000    
    ## 13 L4 vs L5  1 0.119605442 0.5077641 0.07802436   0.948    0.97300    
    ## 14 L4 vs WT  1 0.104794282 0.9125615 0.13201496   0.810    0.97300    
    ## 15 L5 vs WT  1 0.132158463 0.9590280 0.13781063   0.455    0.97300    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.08909 0.017819 0.5947    999  0.737
    ## Residuals 18 0.53932 0.029962                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##         L1      L2      L3      L4      L5    WT
    ## L1         0.90200 0.62600 0.12900 0.63700 0.736
    ## L2 0.81998         0.70100 0.09900 0.65700 0.621
    ## L3 0.49523 0.56756         0.21000 0.92000 0.486
    ## L4 0.17602 0.15825 0.23053         0.37400 0.325
    ## L5 0.52245 0.56267 0.82293 0.34582         0.577
    ## WT 0.61690 0.49936 0.38385 0.31566 0.45527

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      2.46     5 0.783 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      0.15 0.881     1 ns          
    ##  2 Shannon WT     L2         4     4     -0.2  0.841     1 ns          
    ##  3 Shannon WT     L3         4     4     -0.3  0.764     1 ns          
    ##  4 Shannon WT     L4         4     4     -1.25 0.211     1 ns          
    ##  5 Shannon WT     L5         4     4     -0.5  0.617     1 ns          
    ##  6 Shannon L1     L2         4     4     -0.35 0.726     1 ns          
    ##  7 Shannon L1     L3         4     4     -0.45 0.653     1 ns          
    ##  8 Shannon L1     L4         4     4     -1.4  0.162     1 ns          
    ##  9 Shannon L1     L5         4     4     -0.65 0.516     1 ns          
    ## 10 Shannon L2     L3         4     4     -0.1  0.920     1 ns          
    ## 11 Shannon L2     L4         4     4     -1.05 0.294     1 ns          
    ## 12 Shannon L2     L5         4     4     -0.3  0.764     1 ns          
    ## 13 Shannon L3     L4         4     4     -0.95 0.342     1 ns          
    ## 14 Shannon L3     L5         4     4     -0.2  0.841     1 ns          
    ## 15 Shannon L4     L5         4     4      0.75 0.453     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-20-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      8.31     5  0.14 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    0.150  0.881   1      ns          
    ##  2 Chao1 WT     L2         4     4    0.100  0.920   1      ns          
    ##  3 Chao1 WT     L3         4     4    1.15   0.249   1      ns          
    ##  4 Chao1 WT     L4         4     4   -1.65   0.0984  1      ns          
    ##  5 Chao1 WT     L5         4     4    0.250  0.802   1      ns          
    ##  6 Chao1 L1     L2         4     4   -0.0501 0.960   1      ns          
    ##  7 Chao1 L1     L3         4     4    1.00   0.316   1      ns          
    ##  8 Chao1 L1     L4         4     4   -1.80   0.0714  1      ns          
    ##  9 Chao1 L1     L5         4     4    0.100  0.920   1      ns          
    ## 10 Chao1 L2     L3         4     4    1.05   0.293   1      ns          
    ## 11 Chao1 L2     L4         4     4   -1.75   0.0796  1      ns          
    ## 12 Chao1 L2     L5         4     4    0.150  0.881   1      ns          
    ## 13 Chao1 L3     L4         4     4   -2.80   0.00503 0.0755 ns          
    ## 14 Chao1 L3     L5         4     4   -0.902  0.367   1      ns          
    ## 15 Chao1 L4     L5         4     4    1.90   0.0570  0.855  ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-20-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-20-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-20-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_f_root_filtered)
```

    ##                R2 Pr(>F)  
    ## genotype  0.25016  0.044 *
    ## Residuals 0.74984         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.1920416 0.7981851 0.1174115   0.916  0.9160000    
    ## 2  L1 vs L3  1 0.2442521 1.0293695 0.1464384   0.391  0.6516667    
    ## 3  L1 vs L4  1 0.2379813 1.0216817 0.1455038   0.495  0.6761538    
    ## 4  L1 vs L5  1 0.3153431 1.3474687 0.1833922   0.077  0.2460000    
    ## 5  L1 vs WT  1 0.2221933 0.8856615 0.1286240   0.633  0.6782143    
    ## 6  L2 vs L3  1 0.2174793 0.9515836 0.1368873   0.574  0.6761538    
    ## 7  L2 vs L4  1 0.2374445 1.0591104 0.1500345   0.332  0.6225000    
    ## 8  L2 vs L5  1 0.2845479 1.2630423 0.1738999   0.162  0.3535714    
    ## 9  L2 vs WT  1 0.2303080 0.9511364 0.1368318   0.569  0.6761538    
    ## 10 L3 vs L4  1 0.3079258 1.3941006 0.1885423   0.025  0.2460000    
    ## 11 L3 vs L5  1 0.4384869 1.9754060 0.2476872   0.036  0.2460000    
    ## 12 L3 vs WT  1 0.2305881 0.9655099 0.1386130   0.586  0.6761538    
    ## 13 L4 vs L5  1 0.2585458 1.1880561 0.1652820   0.165  0.3535714    
    ## 14 L4 vs WT  1 0.3481445 1.4847960 0.1983749   0.082  0.2460000    
    ## 15 L5 vs WT  1 0.4205534 1.7852719 0.2293140   0.066  0.2460000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.004154 0.00083079 0.5396    999   0.76
    ## Residuals 18 0.027714 0.00153967                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.952000 0.912000 0.389000 0.559000 0.974
    ## L2 0.960346          0.857000 0.233000 0.503000 0.961
    ## L3 0.912210 0.821320          0.114000 0.368000 0.811
    ## L4 0.403459 0.220660 0.122166          0.740000 0.013
    ## L5 0.573785 0.479506 0.345512 0.722421          0.331
    ## WT 0.977922 0.961388 0.788942 0.010729 0.315286

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      6.47     5 0.263 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      1.35 0.177  1     ns          
    ##  2 Shannon WT     L2         4     4      1.25 0.211  1     ns          
    ##  3 Shannon WT     L3         4     4      0.4  0.689  1     ns          
    ##  4 Shannon WT     L4         4     4      2.05 0.0404 0.605 ns          
    ##  5 Shannon WT     L5         4     4      1.85 0.0643 0.965 ns          
    ##  6 Shannon L1     L2         4     4     -0.1  0.920  1     ns          
    ##  7 Shannon L1     L3         4     4     -0.95 0.342  1     ns          
    ##  8 Shannon L1     L4         4     4      0.7  0.484  1     ns          
    ##  9 Shannon L1     L5         4     4      0.5  0.617  1     ns          
    ## 10 Shannon L2     L3         4     4     -0.85 0.395  1     ns          
    ## 11 Shannon L2     L4         4     4      0.8  0.424  1     ns          
    ## 12 Shannon L2     L5         4     4      0.6  0.549  1     ns          
    ## 13 Shannon L3     L4         4     4      1.65 0.0989 1     ns          
    ## 14 Shannon L3     L5         4     4      1.45 0.147  1     ns          
    ## 15 Shannon L4     L5         4     4     -0.2  0.841  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      8.14     5 0.149 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    1.15   0.250  1     ns          
    ##  2 Chao1 WT     L2         4     4    0.775  0.438  1     ns          
    ##  3 Chao1 WT     L3         4     4    0.200  0.841  1     ns          
    ##  4 Chao1 WT     L4         4     4    2.13   0.0335 0.503 ns          
    ##  5 Chao1 WT     L5         4     4    2.05   0.0403 0.605 ns          
    ##  6 Chao1 L1     L2         4     4   -0.375  0.708  1     ns          
    ##  7 Chao1 L1     L3         4     4   -0.950  0.342  1     ns          
    ##  8 Chao1 L1     L4         4     4    0.975  0.329  1     ns          
    ##  9 Chao1 L1     L5         4     4    0.900  0.368  1     ns          
    ## 10 Chao1 L2     L3         4     4   -0.575  0.565  1     ns          
    ## 11 Chao1 L2     L4         4     4    1.35   0.177  1     ns          
    ## 12 Chao1 L2     L5         4     4    1.28   0.202  1     ns          
    ## 13 Chao1 L3     L4         4     4    1.93   0.0542 0.813 ns          
    ## 14 Chao1 L3     L5         4     4    1.85   0.0643 0.964 ns          
    ## 15 Chao1 L4     L5         4     4   -0.0750 0.940  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-21-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-21-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_f_root_associated_soils_filtered)
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

    ##                R2 Pr(>F)   
    ## genotype  0.28926   0.01 **
    ## Residuals 0.71074          
    ## Total     1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.1750592 0.8997805 0.13040712   0.481  0.5700000    
    ## 2  L1 vs L3  1 0.2524964 1.3844445 0.18748119   0.169  0.3168750    
    ## 3  L1 vs L4  1 0.1957800 0.9976979 0.14257516   0.369  0.5535000    
    ## 4  L1 vs L5  1 0.2571296 1.4274315 0.19218373   0.064  0.1675000    
    ## 5  L1 vs WT  1 0.3300238 1.3743853 0.18637286   0.092  0.1971429    
    ## 6  L2 vs L3  1 0.1588754 0.8966835 0.13001662   0.494  0.5700000    
    ## 7  L2 vs L4  1 0.1234730 0.6463475 0.09724853   0.956  0.9560000    
    ## 8  L2 vs L5  1 0.1767170 1.0101890 0.14410296   0.355  0.5535000    
    ## 9  L2 vs WT  1 0.5104425 2.1727915 0.26585672   0.024  0.1537500    
    ## 10 L3 vs L4  1 0.1510574 0.8445793 0.12339390   0.837  0.8967857    
    ## 11 L3 vs L5  1 0.2184563 1.3422163 0.18280806   0.067  0.1675000    
    ## 12 L3 vs WT  1 0.5762198 2.5868682 0.30125863   0.030  0.1537500    
    ## 13 L4 vs L5  1 0.1781923 1.0089672 0.14395377   0.446  0.5700000    
    ## 14 L4 vs WT  1 0.5111378 2.1603568 0.26473803   0.032  0.1537500    
    ## 15 L5 vs WT  1 0.5357739 2.4297971 0.28823909   0.041  0.1537500    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.022134 0.0044268 1.7523    999  0.154
    ## Residuals 18 0.045473 0.0025263                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.699000 0.425000 0.865000 0.276000 0.008
    ## L2 0.714889          0.692000 0.852000 0.630000 0.075
    ## L3 0.453831 0.703272          0.565000 0.977000 0.060
    ## L4 0.865695 0.858462 0.580594          0.477000 0.063
    ## L5 0.293005 0.653152 0.982466 0.494298          0.021
    ## WT 0.003895 0.063500 0.060001 0.054393 0.014078      
    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      6.41     5 0.268 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      0.35 0.726  1     ns          
    ##  2 Shannon WT     L2         4     4     -0.15 0.881  1     ns          
    ##  3 Shannon WT     L3         4     4      0.7  0.484  1     ns          
    ##  4 Shannon WT     L4         4     4      0.95 0.342  1     ns          
    ##  5 Shannon WT     L5         4     4      2.05 0.0404 0.605 ns          
    ##  6 Shannon L1     L2         4     4     -0.5  0.617  1     ns          
    ##  7 Shannon L1     L3         4     4      0.35 0.726  1     ns          
    ##  8 Shannon L1     L4         4     4      0.6  0.549  1     ns          
    ##  9 Shannon L1     L5         4     4      1.7  0.0891 1     ns          
    ## 10 Shannon L2     L3         4     4      0.85 0.395  1     ns          
    ## 11 Shannon L2     L4         4     4      1.1  0.271  1     ns          
    ## 12 Shannon L2     L5         4     4      2.2  0.0278 0.417 ns          
    ## 13 Shannon L3     L4         4     4      0.25 0.803  1     ns          
    ## 14 Shannon L3     L5         4     4      1.35 0.177  1     ns          
    ## 15 Shannon L4     L5         4     4      1.1  0.271  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-22-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      8.63     5 0.125 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4      1.15 0.250   1      ns          
    ##  2 Chao1 WT     L2         4     4      0.8  0.424   1      ns          
    ##  3 Chao1 WT     L3         4     4      1.3  0.194   1      ns          
    ##  4 Chao1 WT     L4         4     4      1.8  0.0719  1      ns          
    ##  5 Chao1 WT     L5         4     4      2.75 0.00596 0.0894 ns          
    ##  6 Chao1 L1     L2         4     4     -0.35 0.726   1      ns          
    ##  7 Chao1 L1     L3         4     4      0.15 0.881   1      ns          
    ##  8 Chao1 L1     L4         4     4      0.65 0.516   1      ns          
    ##  9 Chao1 L1     L5         4     4      1.6  0.110   1      ns          
    ## 10 Chao1 L2     L3         4     4      0.5  0.617   1      ns          
    ## 11 Chao1 L2     L4         4     4      1    0.317   1      ns          
    ## 12 Chao1 L2     L5         4     4      1.95 0.0512  0.768  ns          
    ## 13 Chao1 L3     L4         4     4      0.5  0.617   1      ns          
    ## 14 Chao1 L3     L5         4     4      1.45 0.147   1      ns          
    ## 15 Chao1 L4     L5         4     4      0.95 0.342   1      ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-22-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-22-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-22-5.png)<!-- -->

``` r
analysis_f(ps1_ssu_f_soil_filtered)
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

    ##                R2 Pr(>F)
    ## genotype  0.21682  0.541
    ## Residuals 0.78318       
    ## Total     1.00000       
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.3733542 0.9630581 0.1383096   0.830  0.9203571    
    ## 2  L1 vs L3  1 0.3804045 1.0089000 0.1439456   0.389  0.8754545    
    ## 3  L1 vs L4  1 0.3678560 0.9771701 0.1400525   0.591  0.8754545    
    ## 4  L1 vs L5  1 0.3997888 1.0469205 0.1485643   0.196  0.8754545    
    ## 5  L1 vs WT  1 0.3565308 0.9750888 0.1397959   0.624  0.8754545    
    ## 6  L2 vs L3  1 0.3634037 0.9804633 0.1404582   0.610  0.8754545    
    ## 7  L2 vs L4  1 0.3637993 0.9831184 0.1407850   0.465  0.8754545    
    ## 8  L2 vs L5  1 0.3777168 1.0059915 0.1435902   0.339  0.8754545    
    ## 9  L2 vs WT  1 0.3459407 0.9629920 0.1383015   0.855  0.9203571    
    ## 10 L3 vs L4  1 0.3517941 0.9787844 0.1402514   0.642  0.8754545    
    ## 11 L3 vs L5  1 0.3968615 1.0877678 0.1534711   0.054  0.5175000    
    ## 12 L3 vs WT  1 0.3318536 0.9519379 0.1369313   0.859  0.9203571    
    ## 13 L4 vs L5  1 0.3603105 0.9892071 0.1415335   0.511  0.8754545    
    ## 14 L4 vs WT  1 0.3262595 0.9375007 0.1351352   0.928  0.9280000    
    ## 15 L5 vs WT  1 0.3881723 1.0982979 0.1547269   0.069  0.5175000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.003523 0.00070451 0.2825    999  0.946
    ## Residuals 18 0.044891 0.00249395                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.564000 0.289000 0.734000 0.454000 0.031
    ## L2 0.526470          0.761000 0.953000 0.956000 0.431
    ## L3 0.289042 0.723183          0.884000 0.755000 0.679
    ## L4 0.694360 0.949724 0.876510          0.985000 0.723
    ## L5 0.418824 0.949425 0.749390 0.977879          0.362
    ## WT 0.033549 0.382258 0.649460 0.684313 0.364516      
    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      3.37     5 0.643 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -1.15 0.250     1 ns          
    ##  2 Shannon WT     L2         4     4     -0.9  0.368     1 ns          
    ##  3 Shannon WT     L3         4     4      0    1         1 ns          
    ##  4 Shannon WT     L4         4     4      0.3  0.764     1 ns          
    ##  5 Shannon WT     L5         4     4     -0.65 0.516     1 ns          
    ##  6 Shannon L1     L2         4     4      0.25 0.803     1 ns          
    ##  7 Shannon L1     L3         4     4      1.15 0.250     1 ns          
    ##  8 Shannon L1     L4         4     4      1.45 0.147     1 ns          
    ##  9 Shannon L1     L5         4     4      0.5  0.617     1 ns          
    ## 10 Shannon L2     L3         4     4      0.9  0.368     1 ns          
    ## 11 Shannon L2     L4         4     4      1.2  0.230     1 ns          
    ## 12 Shannon L2     L5         4     4      0.25 0.803     1 ns          
    ## 13 Shannon L3     L4         4     4      0.3  0.764     1 ns          
    ## 14 Shannon L3     L5         4     4     -0.65 0.516     1 ns          
    ## 15 Shannon L4     L5         4     4     -0.95 0.342     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      3.38     5 0.642 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     -1.2  0.230     1 ns          
    ##  2 Chao1 WT     L2         4     4     -0.85 0.395     1 ns          
    ##  3 Chao1 WT     L3         4     4      0.05 0.960     1 ns          
    ##  4 Chao1 WT     L4         4     4      0.25 0.803     1 ns          
    ##  5 Chao1 WT     L5         4     4     -0.65 0.516     1 ns          
    ##  6 Chao1 L1     L2         4     4      0.35 0.726     1 ns          
    ##  7 Chao1 L1     L3         4     4      1.25 0.211     1 ns          
    ##  8 Chao1 L1     L4         4     4      1.45 0.147     1 ns          
    ##  9 Chao1 L1     L5         4     4      0.55 0.582     1 ns          
    ## 10 Chao1 L2     L3         4     4      0.9  0.368     1 ns          
    ## 11 Chao1 L2     L4         4     4      1.1  0.271     1 ns          
    ## 12 Chao1 L2     L5         4     4      0.2  0.841     1 ns          
    ## 13 Chao1 L3     L4         4     4      0.2  0.841     1 ns          
    ## 14 Chao1 L3     L5         4     4     -0.7  0.484     1 ns          
    ## 15 Chao1 L4     L5         4     4     -0.9  0.368     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-23-5.png)<!-- -->

``` r
analysis_f(ps1_its_g_leaf_filtered)
```

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ## Warning: Removed 2 rows containing missing values (geom_point).

    ##               R2 Pr(>F)
    ## genotype  0.2051  0.753
    ## Residuals 0.7949       
    ## Total     1.0000

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.08273896 0.3144260 0.05916463   0.940  1.0000000    
    ## 2  L1 vs L3  1 0.25847464 0.9630253 0.13830559   0.297  0.8745000    
    ## 3  L1 vs L4  1 0.45812009 1.2170586 0.16863638   0.312  0.8745000    
    ## 4  L1 vs L5  1 0.33321054 1.1390851 0.18554640   0.363  0.8745000    
    ## 5  L1 vs WT  1 0.12559350 0.4049082 0.06321843   1.000  1.0000000    
    ## 6  L2 vs L3  1 0.24056330 0.8682992 0.14796437   0.506  0.8745000    
    ## 7  L2 vs L4  1 0.30804463 0.7574772 0.13156408   0.583  0.8745000    
    ## 8  L2 vs L5  1 0.24911567 0.8052303 0.16757372   0.700  0.9545455    
    ## 9  L2 vs WT  1 0.14067056 0.4299405 0.07917960   1.000  1.0000000    
    ## 10 L3 vs L4  1 0.40903265 1.0541913 0.14944184   0.488  0.8745000    
    ## 11 L3 vs L5  1 0.17294591 0.5643840 0.10142794   0.944  1.0000000    
    ## 12 L3 vs WT  1 0.29796559 0.9260262 0.13370238   0.372  0.8745000    
    ## 13 L4 vs L5  1 0.40434936 0.9272926 0.15644455   0.567  0.8745000    
    ## 14 L4 vs WT  1 0.36672491 0.8532749 0.12450615   0.514  0.8745000    
    ## 15 L5 vs WT  1 0.31484349 0.8829831 0.15009104   0.571  0.8745000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)  
    ## Groups     5 0.25380 0.050760 2.344    999  0.085 .
    ## Residuals 16 0.34649 0.021656                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.6160000 0.8300000 0.0170000 0.4320000 0.401
    ## L2 0.5930029           0.3860000 0.0220000 0.8570000 0.692
    ## L3 0.8326575 0.4033990           0.0030000 0.1960000 0.282
    ## L4 0.0117554 0.0260991 0.0018297           0.0030000 0.239
    ## L5 0.4355654 0.8655759 0.2083351 0.0008254           0.749
    ## WT 0.3984899 0.7014967 0.2949052 0.2419961 0.7515285

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    22      2.56     5 0.767 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4    -0.408 0.683     1 ns          
    ##  2 Shannon WT     L2         4     3     0.437 0.662     1 ns          
    ##  3 Shannon WT     L3         4     4    -0.408 0.683     1 ns          
    ##  4 Shannon WT     L4         4     4    -0.980 0.327     1 ns          
    ##  5 Shannon WT     L5         4     3     0.303 0.762     1 ns          
    ##  6 Shannon L1     L2         4     3     0.815 0.415     1 ns          
    ##  7 Shannon L1     L3         4     4     0     1         1 ns          
    ##  8 Shannon L1     L4         4     4    -0.572 0.567     1 ns          
    ##  9 Shannon L1     L5         4     3     0.681 0.496     1 ns          
    ## 10 Shannon L2     L3         3     4    -0.815 0.415     1 ns          
    ## 11 Shannon L2     L4         3     4    -1.34  0.179     1 ns          
    ## 12 Shannon L2     L5         3     3    -0.126 0.900     1 ns          
    ## 13 Shannon L3     L4         4     4    -0.572 0.567     1 ns          
    ## 14 Shannon L3     L5         4     3     0.681 0.496     1 ns          
    ## 15 Shannon L4     L5         4     3     1.21  0.226     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-24-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    22      4.03     5 0.546 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4  -0.494   0.621  1     ns          
    ##  2 Chao1 WT     L2         4     3   0.00847 0.993  1     ns          
    ##  3 Chao1 WT     L3         4     4   0.467   0.641  1     ns          
    ##  4 Chao1 WT     L4         4     4  -1.37    0.170  1     ns          
    ##  5 Chao1 WT     L5         4     3   0.0424  0.966  1     ns          
    ##  6 Chao1 L1     L2         4     3   0.466   0.641  1     ns          
    ##  7 Chao1 L1     L3         4     4   0.961   0.337  1     ns          
    ##  8 Chao1 L1     L4         4     4  -0.879   0.380  1     ns          
    ##  9 Chao1 L1     L5         4     3   0.500   0.617  1     ns          
    ## 10 Chao1 L2     L3         3     4   0.424   0.672  1     ns          
    ## 11 Chao1 L2     L4         3     4  -1.28    0.201  1     ns          
    ## 12 Chao1 L2     L5         3     3   0.0317  0.975  1     ns          
    ## 13 Chao1 L3     L4         4     4  -1.84    0.0658 0.987 ns          
    ## 14 Chao1 L3     L5         4     3  -0.390   0.697  1     ns          
    ## 15 Chao1 L4     L5         4     3   1.31    0.189  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-24-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-24-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-24-5.png)<!-- -->

``` r
analysis_f(ps1_its_g_stem_filtered)
```

    ##                R2 Pr(>F)
    ## genotype  0.24485  0.241
    ## Residuals 0.75515       
    ## Total     1.00000

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.4049447 0.9694466 0.13909951   0.547  0.7459091    
    ## 2  L1 vs L3  1 0.5111959 1.2123845 0.16809759   0.237  0.4283333    
    ## 3  L1 vs L4  1 0.3558904 0.8359560 0.12228809   0.636  0.7950000    
    ## 4  L1 vs L5  1 0.3896624 0.8382825 0.14358375   0.506  0.7459091    
    ## 5  L1 vs WT  1 0.6537635 1.6093582 0.21149723   0.046  0.3000000    
    ## 6  L2 vs L3  1 0.2383428 0.5989088 0.09075876   0.838  0.9669231    
    ## 7  L2 vs L4  1 0.2154886 0.5359807 0.08200463   1.000  1.0000000    
    ## 8  L2 vs L5  1 0.6241920 1.4302739 0.22242815   0.049  0.3000000    
    ## 9  L2 vs WT  1 0.5259976 1.3750025 0.18644095   0.136  0.3400000    
    ## 10 L3 vs L4  1 0.2299749 0.5664638 0.08626619   0.958  1.0000000    
    ## 11 L3 vs L5  1 0.5688036 1.2893951 0.20501099   0.112  0.3360000    
    ## 12 L3 vs WT  1 0.4233687 1.0954453 0.15438711   0.257  0.4283333    
    ## 13 L4 vs L5  1 0.6049782 1.3563310 0.21338269   0.248  0.4283333    
    ## 14 L4 vs WT  1 0.5403862 1.3836025 0.18738855   0.080  0.3000000    
    ## 15 L5 vs WT  1 0.6438064 1.5233072 0.23351763   0.063  0.3000000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.011429 0.0022858 0.1255    999  0.983
    ## Residuals 17 0.309530 0.0182076                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##         L1      L2      L3      L4      L5    WT
    ## L1         0.82600 0.92300 0.99700 0.47800 0.845
    ## L2 0.83174         0.71200 0.88200 0.46100 0.969
    ## L3 0.92654 0.72201         0.95500 0.32000 0.744
    ## L4 0.99846 0.88103 0.94499         0.63900 0.914
    ## L5 0.47750 0.48246 0.31732 0.64594         0.484
    ## WT 0.85822 0.96794 0.75210 0.89991 0.46161

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    23      4.36     5 0.499 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4   -1.25   0.210     1 ns          
    ##  2 Shannon WT     L2         4     4   -0.522  0.602     1 ns          
    ##  3 Shannon WT     L3         4     4   -0.157  0.876     1 ns          
    ##  4 Shannon WT     L4         4     4   -1.30   0.192     1 ns          
    ##  5 Shannon WT     L5         4     3   -1.56   0.118     1 ns          
    ##  6 Shannon L1     L2         4     4    0.731  0.465     1 ns          
    ##  7 Shannon L1     L3         4     4    1.10   0.273     1 ns          
    ##  8 Shannon L1     L4         4     4   -0.0522 0.958     1 ns          
    ##  9 Shannon L1     L5         4     3   -0.403  0.687     1 ns          
    ## 10 Shannon L2     L3         4     4    0.365  0.715     1 ns          
    ## 11 Shannon L2     L4         4     4   -0.783  0.434     1 ns          
    ## 12 Shannon L2     L5         4     3   -1.08   0.281     1 ns          
    ## 13 Shannon L3     L4         4     4   -1.15   0.251     1 ns          
    ## 14 Shannon L3     L5         4     3   -1.42   0.156     1 ns          
    ## 15 Shannon L4     L5         4     3   -0.354  0.723     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-25-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    23      3.86     5  0.57 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4   -0.739  0.460      1 ns          
    ##  2 Chao1 WT     L2         4     4   -0.580  0.562      1 ns          
    ##  3 Chao1 WT     L3         4     4   -0.501  0.616      1 ns          
    ##  4 Chao1 WT     L4         4     4   -1.32   0.187      1 ns          
    ##  5 Chao1 WT     L5         4     3   -1.74   0.0815     1 ns          
    ##  6 Chao1 L1     L2         4     4    0.158  0.874      1 ns          
    ##  7 Chao1 L1     L3         4     4    0.237  0.812      1 ns          
    ##  8 Chao1 L1     L4         4     4   -0.580  0.562      1 ns          
    ##  9 Chao1 L1     L5         4     3   -1.06   0.290      1 ns          
    ## 10 Chao1 L2     L3         4     4    0.0791 0.937      1 ns          
    ## 11 Chao1 L2     L4         4     4   -0.739  0.460      1 ns          
    ## 12 Chao1 L2     L5         4     3   -1.20   0.228      1 ns          
    ## 13 Chao1 L3     L4         4     4   -0.818  0.414      1 ns          
    ## 14 Chao1 L3     L5         4     3   -1.28   0.201      1 ns          
    ## 15 Chao1 L4     L5         4     3   -0.521  0.602      1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-25-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-25-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-25-5.png)<!-- -->

``` r
analysis_f(ps1_its_g_root_filtered)
```

    ##                R2 Pr(>F)
    ## genotype  0.23765  0.522
    ## Residuals 0.76235       
    ## Total     1.00000

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.5335644 1.1494200 0.1607711   0.293  0.8162500    
    ## 2  L1 vs L3  1 0.4994795 1.0212873 0.2033915   0.400  0.8162500    
    ## 3  L1 vs L4  1 0.4978088 1.0698883 0.1513303   0.360  0.8162500    
    ## 4  L1 vs L5  1 0.5154235 1.0466512 0.1485317   0.163  0.8162500    
    ## 5  L1 vs WT  1 0.4276056 0.8840936 0.1284256   0.738  0.8515385    
    ## 6  L2 vs L3  1 0.4247375 0.9304589 0.1887165   0.600  0.8162500    
    ## 7  L2 vs L4  1 0.4433919 0.9996076 0.1428091   0.477  0.8162500    
    ## 8  L2 vs L5  1 0.5453354 1.1584995 0.1618355   0.216  0.8162500    
    ## 9  L2 vs WT  1 0.4237253 0.9172714 0.1326060   0.653  0.8162500    
    ## 10 L3 vs L4  1 0.4463686 0.9743643 0.1958771   0.600  0.8162500    
    ## 11 L3 vs L5  1 0.4479874 0.8980364 0.1833462   1.000  1.0000000    
    ## 12 L3 vs WT  1 0.4589659 0.9450056 0.1911030   0.600  0.8162500    
    ## 13 L4 vs L5  1 0.4732401 1.0030249 0.1432274   0.416  0.8162500    
    ## 14 L4 vs WT  1 0.4191363 0.9052066 0.1310904   0.860  0.9214286    
    ## 15 L5 vs WT  1 0.4926273 1.0049762 0.1434660   0.198  0.8162500    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.021115 0.0042230 0.8149    999  0.592
    ## Residuals 16 0.082918 0.0051824                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##            L1         L2         L3         L4         L5    WT
    ## L1            0.85100000 0.02800000 0.85100000 0.18100000 0.410
    ## L2 0.85604940            0.49600000 0.95300000 0.53700000 0.647
    ## L3 0.03999197 0.50701324            0.28600000 0.00100000 0.023
    ## L4 0.86406068 0.95363719 0.29121659            0.41100000 0.517
    ## L5 0.17811989 0.54424443 0.00004648 0.40701415            0.811
    ## WT 0.41840457 0.61268750 0.02145380 0.51718334 0.78680471

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    22      3.69     5 0.595 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4    -1.03  0.301     1 ns          
    ##  2 Shannon WT     L2         4     4     0.218 0.828     1 ns          
    ##  3 Shannon WT     L3         4     2    -0.178 0.859     1 ns          
    ##  4 Shannon WT     L4         4     4     0.544 0.586     1 ns          
    ##  5 Shannon WT     L5         4     4    -0.817 0.414     1 ns          
    ##  6 Shannon L1     L2         4     4     1.25  0.210     1 ns          
    ##  7 Shannon L1     L3         4     2     0.667 0.505     1 ns          
    ##  8 Shannon L1     L4         4     4     1.58  0.114     1 ns          
    ##  9 Shannon L1     L5         4     4     0.218 0.828     1 ns          
    ## 10 Shannon L2     L3         4     2    -0.356 0.722     1 ns          
    ## 11 Shannon L2     L4         4     4     0.327 0.744     1 ns          
    ## 12 Shannon L2     L5         4     4    -1.03  0.301     1 ns          
    ## 13 Shannon L3     L4         2     4     0.622 0.534     1 ns          
    ## 14 Shannon L3     L5         2     4    -0.489 0.625     1 ns          
    ## 15 Shannon L4     L5         4     4    -1.36  0.173     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-26-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    22      4.32     5 0.504 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    -0.873 0.383      1 ns          
    ##  2 Chao1 WT     L2         4     4     0.873 0.383      1 ns          
    ##  3 Chao1 WT     L3         4     2     0.423 0.672      1 ns          
    ##  4 Chao1 WT     L4         4     4     0.737 0.461      1 ns          
    ##  5 Chao1 WT     L5         4     4    -0.246 0.806      1 ns          
    ##  6 Chao1 L1     L2         4     4     1.75  0.0808     1 ns          
    ##  7 Chao1 L1     L3         4     2     1.14  0.256      1 ns          
    ##  8 Chao1 L1     L4         4     4     1.61  0.107      1 ns          
    ##  9 Chao1 L1     L5         4     4     0.628 0.530      1 ns          
    ## 10 Chao1 L2     L3         4     2    -0.290 0.772      1 ns          
    ## 11 Chao1 L2     L4         4     4    -0.136 0.891      1 ns          
    ## 12 Chao1 L2     L5         4     4    -1.12  0.263      1 ns          
    ## 13 Chao1 L3     L4         2     4     0.178 0.859      1 ns          
    ## 14 Chao1 L3     L5         2     4    -0.624 0.533      1 ns          
    ## 15 Chao1 L4     L5         4     4    -0.982 0.326      1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-26-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-26-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-26-5.png)<!-- -->

``` r
analysis_f(ps1_its_g_root_associated_soils_filtered)
```

    ##                R2 Pr(>F)  
    ## genotype  0.31672  0.014 *
    ## Residuals 0.68328         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.0816277 0.5673524 0.08638982   0.906  0.9060000    
    ## 2  L1 vs L3  1 0.5410933 3.2192097 0.34918499   0.029  0.1500000    
    ## 3  L1 vs L4  1 0.1230119 0.7758914 0.11450764   0.726  0.7778571    
    ## 4  L1 vs L5  1 0.2653120 1.3706257 0.18595785   0.166  0.3112500    
    ## 5  L1 vs WT  1 0.1212285 0.7379106 0.10951624   0.703  0.7778571    
    ## 6  L2 vs L3  1 0.4408978 2.8535786 0.32230793   0.027  0.1500000    
    ## 7  L2 vs L4  1 0.1685948 1.1629868 0.16236060   0.289  0.4335000    
    ## 8  L2 vs L5  1 0.3070363 1.7058111 0.22136684   0.136  0.3085714    
    ## 9  L2 vs WT  1 0.1006024 0.6675207 0.10011528   0.719  0.7778571    
    ## 10 L3 vs L4  1 0.4842254 2.8622756 0.32297299   0.045  0.1687500    
    ## 11 L3 vs L5  1 0.4162746 2.0385414 0.25359593   0.144  0.3085714    
    ## 12 L3 vs WT  1 0.4327863 2.4742181 0.29197008   0.030  0.1500000    
    ## 13 L4 vs L5  1 0.2291900 1.1773727 0.16403951   0.343  0.4677273    
    ## 14 L4 vs WT  1 0.1995718 1.2067582 0.16744813   0.235  0.3916667    
    ## 15 L5 vs WT  1 0.3726332 1.8593941 0.23658237   0.095  0.2850000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.055934 0.0111868 1.9453    999  0.113
    ## Residuals 18 0.103511 0.0057506                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.157000 0.211000 0.410000 0.337000 0.402
    ## L2 0.164016          0.809000 0.039000 0.094000 0.726
    ## L3 0.223798 0.763110          0.054000 0.139000 0.929
    ## L4 0.381018 0.048245 0.060012          0.549000 0.143
    ## L5 0.322702 0.117358 0.140709 0.525336          0.175
    ## WT 0.395368 0.723018 0.905034 0.156054 0.181366

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-27-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      2.97     5 0.705 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -0.75 0.453     1 ns          
    ##  2 Shannon WT     L2         4     4     -0.2  0.841     1 ns          
    ##  3 Shannon WT     L3         4     4      0.7  0.484     1 ns          
    ##  4 Shannon WT     L4         4     4      0.55 0.582     1 ns          
    ##  5 Shannon WT     L5         4     4     -0.3  0.764     1 ns          
    ##  6 Shannon L1     L2         4     4      0.55 0.582     1 ns          
    ##  7 Shannon L1     L3         4     4      1.45 0.147     1 ns          
    ##  8 Shannon L1     L4         4     4      1.3  0.194     1 ns          
    ##  9 Shannon L1     L5         4     4      0.45 0.653     1 ns          
    ## 10 Shannon L2     L3         4     4      0.9  0.368     1 ns          
    ## 11 Shannon L2     L4         4     4      0.75 0.453     1 ns          
    ## 12 Shannon L2     L5         4     4     -0.1  0.920     1 ns          
    ## 13 Shannon L3     L4         4     4     -0.15 0.881     1 ns          
    ## 14 Shannon L3     L5         4     4     -1    0.317     1 ns          
    ## 15 Shannon L4     L5         4     4     -0.85 0.395     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-27-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      6.03     5 0.304 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    -0.750 0.453  1     ns          
    ##  2 Chao1 WT     L2         4     4    -0.125 0.900  1     ns          
    ##  3 Chao1 WT     L3         4     4     1.25  0.211  1     ns          
    ##  4 Chao1 WT     L4         4     4    -0.575 0.565  1     ns          
    ##  5 Chao1 WT     L5         4     4    -0.850 0.395  1     ns          
    ##  6 Chao1 L1     L2         4     4     0.625 0.532  1     ns          
    ##  7 Chao1 L1     L3         4     4     2.00  0.0454 0.681 ns          
    ##  8 Chao1 L1     L4         4     4     0.175 0.861  1     ns          
    ##  9 Chao1 L1     L5         4     4    -0.100 0.920  1     ns          
    ## 10 Chao1 L2     L3         4     4     1.38  0.169  1     ns          
    ## 11 Chao1 L2     L4         4     4    -0.450 0.653  1     ns          
    ## 12 Chao1 L2     L5         4     4    -0.725 0.468  1     ns          
    ## 13 Chao1 L3     L4         4     4    -1.83  0.0679 1     ns          
    ## 14 Chao1 L3     L5         4     4    -2.10  0.0356 0.535 ns          
    ## 15 Chao1 L4     L5         4     4    -0.275 0.783  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-27-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-27-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-27-5.png)<!-- -->

``` r
analysis_f(ps1_its_g_soil_filtered)
```

    ##                R2 Pr(>F)  
    ## genotype  0.28943  0.079 .
    ## Residuals 0.71057         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.07503769 0.6677369 0.10014446   0.815   0.905000    
    ## 2  L1 vs L3  1 0.16402191 0.8309556 0.12164558   0.573   0.716250    
    ## 3  L1 vs L4  1 0.07326880 0.6242008 0.09423036   0.873   0.905000    
    ## 4  L1 vs L5  1 0.14834069 0.8292855 0.12143078   0.560   0.716250    
    ## 5  L1 vs WT  1 0.23874162 1.8058672 0.23134742   0.107   0.321000    
    ## 6  L2 vs L3  1 0.24697405 1.5272483 0.20289596   0.268   0.556875    
    ## 7  L2 vs L4  1 0.09003669 1.1020073 0.15516842   0.394   0.591000    
    ## 8  L2 vs L5  1 0.27160892 1.8967105 0.24018995   0.209   0.522500    
    ## 9  L2 vs WT  1 0.15901338 1.6473704 0.21541658   0.072   0.277500    
    ## 10 L3 vs L4  1 0.21040745 1.2620723 0.17378956   0.356   0.591000    
    ## 11 L3 vs L5  1 0.04338807 0.1901206 0.03071356   0.905   0.905000    
    ## 12 L3 vs WT  1 0.50238742 2.7673806 0.31564508   0.050   0.277500    
    ## 13 L4 vs L5  1 0.21691629 1.4636328 0.19610193   0.297   0.556875    
    ## 14 L4 vs WT  1 0.28803365 2.8369417 0.32103207   0.032   0.277500    
    ## 15 L5 vs WT  1 0.51338206 3.1490583 0.34419480   0.074   0.277500    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)  
    ## Groups     5 0.084866 0.016973 3.1767    999  0.035 *
    ## Residuals 18 0.096173 0.005343                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.323000 0.097000 0.433000 0.191000 0.985
    ## L2 0.299229          0.026000 0.592000 0.048000 0.423
    ## L3 0.087944 0.019147          0.023000 0.652000 0.108
    ## L4 0.406042 0.579201 0.019127          0.044000 0.523
    ## L5 0.193295 0.043749 0.631842 0.047933          0.240
    ## WT 0.986457 0.397002 0.118424 0.528199 0.236457

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24     0.700     5 0.983 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -0.65 0.516     1 ns          
    ##  2 Shannon WT     L2         4     4     -0.2  0.841     1 ns          
    ##  3 Shannon WT     L3         4     4     -0.55 0.582     1 ns          
    ##  4 Shannon WT     L4         4     4     -0.05 0.960     1 ns          
    ##  5 Shannon WT     L5         4     4     -0.35 0.726     1 ns          
    ##  6 Shannon L1     L2         4     4      0.45 0.653     1 ns          
    ##  7 Shannon L1     L3         4     4      0.1  0.920     1 ns          
    ##  8 Shannon L1     L4         4     4      0.6  0.549     1 ns          
    ##  9 Shannon L1     L5         4     4      0.3  0.764     1 ns          
    ## 10 Shannon L2     L3         4     4     -0.35 0.726     1 ns          
    ## 11 Shannon L2     L4         4     4      0.15 0.881     1 ns          
    ## 12 Shannon L2     L5         4     4     -0.15 0.881     1 ns          
    ## 13 Shannon L3     L4         4     4      0.5  0.617     1 ns          
    ## 14 Shannon L3     L5         4     4      0.2  0.841     1 ns          
    ## 15 Shannon L4     L5         4     4     -0.3  0.764     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      4.19     5 0.522 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4      1.3  0.194  1     ns          
    ##  2 Chao1 WT     L2         4     4      1.2  0.230  1     ns          
    ##  3 Chao1 WT     L3         4     4      0.95 0.342  1     ns          
    ##  4 Chao1 WT     L4         4     4      2    0.0455 0.683 ns          
    ##  5 Chao1 WT     L5         4     4      1.15 0.250  1     ns          
    ##  6 Chao1 L1     L2         4     4     -0.1  0.920  1     ns          
    ##  7 Chao1 L1     L3         4     4     -0.35 0.726  1     ns          
    ##  8 Chao1 L1     L4         4     4      0.7  0.484  1     ns          
    ##  9 Chao1 L1     L5         4     4     -0.15 0.881  1     ns          
    ## 10 Chao1 L2     L3         4     4     -0.25 0.803  1     ns          
    ## 11 Chao1 L2     L4         4     4      0.8  0.424  1     ns          
    ## 12 Chao1 L2     L5         4     4     -0.05 0.960  1     ns          
    ## 13 Chao1 L3     L4         4     4      1.05 0.294  1     ns          
    ## 14 Chao1 L3     L5         4     4      0.2  0.841  1     ns          
    ## 15 Chao1 L4     L5         4     4     -0.85 0.395  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-28-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-28-5.png)<!-- -->

``` r
analysis_f(ps1_its_f_leaf_filtered)
```

    ##                R2 Pr(>F)
    ## genotype  0.23816  0.677
    ## Residuals 0.76184       
    ## Total     1.00000

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ## 'nperm' >= set of all permutations: complete enumeration.

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.5316324 1.1825969 0.1912783   0.152      0.570    
    ## 2  L1 vs L3  1 0.6560060 1.4720312 0.1970055   0.033      0.275    
    ## 3  L1 vs L4  1 0.6029468 1.3388779 0.2112169   0.055      0.275    
    ## 4  L1 vs L5  1 0.3611829 0.8013458 0.1381310   0.710      1.000    
    ## 5  L1 vs WT  1 0.5880380 1.2807776 0.1759122   0.052      0.275    
    ## 6  L2 vs L3  1 0.3479059 0.7216906 0.1261324   0.822      1.000    
    ## 7  L2 vs L4  1 0.3413118 0.6866889 0.1465190   1.000      1.000    
    ## 8  L2 vs L5  1 0.4280222 0.8603138 0.1770079   0.800      1.000    
    ## 9  L2 vs WT  1 0.3599996 0.7225346 0.1262613   0.973      1.000    
    ## 10 L3 vs L4  1 0.3498350 0.7245039 0.1265619   0.871      1.000    
    ## 11 L3 vs L5  1 0.4627071 0.9575008 0.1607219   0.539      1.000    
    ## 12 L3 vs WT  1 0.3410043 0.7013242 0.1046546   0.851      1.000    
    ## 13 L4 vs L5  1 0.4942724 0.9915049 0.1986385   0.600      1.000    
    ## 14 L4 vs WT  1 0.3910597 0.7836298 0.1354910   0.465      1.000    
    ## 15 L5 vs WT  1 0.4122536 0.8254659 0.1416995   0.889      1.000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.011373 0.0022746 0.5736    999  0.742
    ## Residuals 15 0.059485 0.0039657                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.7560000 0.7840000 0.6190000 0.5620000 0.154
    ## L2 0.7555431           0.9800000 0.8480000 0.8340000 0.191
    ## L3 0.7778139 0.9793445           0.8440000 0.8300000 0.300
    ## L4 0.6225660 0.8547512 0.8690580           0.8210000 0.004
    ## L5 0.6035757 0.8250129 0.8471324 0.8384910           0.009
    ## WT 0.1688874 0.1861042 0.2920471 0.0059053 0.0063456

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    21      8.00     5 0.156 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4    1.77   0.0772  1     ns          
    ##  2 Shannon WT     L2         4     3    0.457  0.647   1     ns          
    ##  3 Shannon WT     L3         4     4   -0.883  0.377   1     ns          
    ##  4 Shannon WT     L4         4     3    0.422  0.673   1     ns          
    ##  5 Shannon WT     L5         4     3    0.985  0.325   1     ns          
    ##  6 Shannon L1     L2         4     3   -1.18   0.239   1     ns          
    ##  7 Shannon L1     L3         4     4   -2.65   0.00804 0.121 ns          
    ##  8 Shannon L1     L4         4     3   -1.21   0.225   1     ns          
    ##  9 Shannon L1     L5         4     3   -0.651  0.515   1     ns          
    ## 10 Shannon L2     L3         3     4   -1.28   0.202   1     ns          
    ## 11 Shannon L2     L4         3     3   -0.0329 0.974   1     ns          
    ## 12 Shannon L2     L5         3     3    0.494  0.622   1     ns          
    ## 13 Shannon L3     L4         4     3    1.24   0.215   1     ns          
    ## 14 Shannon L3     L5         4     3    1.80   0.0714  1     ns          
    ## 15 Shannon L4     L5         3     3    0.527  0.599   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-29-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    21      4.71     5 0.452 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    1.56   0.119  1     ns          
    ##  2 Chao1 WT     L2         4     3    0.330  0.742  1     ns          
    ##  3 Chao1 WT     L3         4     4   -0.433  0.665  1     ns          
    ##  4 Chao1 WT     L4         4     3    0.294  0.769  1     ns          
    ##  5 Chao1 WT     L5         4     3    0.793  0.428  1     ns          
    ##  6 Chao1 L1     L2         4     3   -1.11   0.265  1     ns          
    ##  7 Chao1 L1     L3         4     4   -1.99   0.0462 0.694 ns          
    ##  8 Chao1 L1     L4         4     3   -1.15   0.250  1     ns          
    ##  9 Chao1 L1     L5         4     3   -0.651  0.515  1     ns          
    ## 10 Chao1 L2     L3         3     4   -0.731  0.465  1     ns          
    ## 11 Chao1 L2     L4         3     3   -0.0334 0.973  1     ns          
    ## 12 Chao1 L2     L5         3     3    0.434  0.665  1     ns          
    ## 13 Chao1 L3     L4         4     3    0.695  0.487  1     ns          
    ## 14 Chao1 L3     L5         4     3    1.19   0.232  1     ns          
    ## 15 Chao1 L4     L5         3     3    0.467  0.641  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-29-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-29-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-29-5.png)<!-- -->

``` r
analysis_f(ps1_its_f_stem_filtered)
```

    ##                R2 Pr(>F)  
    ## genotype  0.30479  0.014 *
    ## Residuals 0.69521         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.6933138 2.6896909 0.30952665   0.094  0.2164286    
    ## 2  L1 vs L3  1 0.2117942 0.7233483 0.10758751   0.773  0.8282143    
    ## 3  L1 vs L4  1 0.5237497 1.5878229 0.20925936   0.174  0.2900000    
    ## 4  L1 vs L5  1 0.3513580 1.1444455 0.16018675   0.164  0.2900000    
    ## 5  L1 vs WT  1 0.2654252 1.6561665 0.21631798   0.048  0.1800000    
    ## 6  L2 vs L3  1 0.6436935 2.0142906 0.25133736   0.088  0.2164286    
    ## 7  L2 vs L4  1 0.4282654 1.2009001 0.16677083   0.282  0.3712500    
    ## 8  L2 vs L5  1 0.5971372 1.7890255 0.22968541   0.101  0.2164286    
    ## 9  L2 vs WT  1 0.8449242 4.5175589 0.42952542   0.024  0.1800000    
    ## 10 L3 vs L4  1 0.4313421 1.1013451 0.15508965   0.297  0.3712500    
    ## 11 L3 vs L5  1 0.2327846 0.6311813 0.09518383   0.919  0.9190000    
    ## 12 L3 vs WT  1 0.2830036 1.2744404 0.17519428   0.292  0.3712500    
    ## 13 L4 vs L5  1 0.4482221 1.1043632 0.15544858   0.435  0.5019231    
    ## 14 L4 vs WT  1 0.5465489 2.1092667 0.26010572   0.032  0.1800000    
    ## 15 L5 vs WT  1 0.4880067 2.0654130 0.25608274   0.041  0.1800000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.14692 0.029383 1.5262    999  0.213
    ## Residuals 18 0.34655 0.019253                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##          L1       L2       L3       L4       L5    WT
    ## L1          0.741000 0.578000 0.084000 0.366000 0.281
    ## L2 0.723443          0.786000 0.152000 0.543000 0.185
    ## L3 0.595137 0.791892          0.345000 0.782000 0.273
    ## L4 0.075023 0.148350 0.332499          0.467000 0.020
    ## L5 0.371706 0.544925 0.777416 0.480844          0.127
    ## WT 0.277538 0.216596 0.256680 0.020145 0.138510

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      7.51     5 0.185 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -1.5  0.134  1     ns          
    ##  2 Shannon WT     L2         4     4     -1.95 0.0512 0.768 ns          
    ##  3 Shannon WT     L3         4     4     -2.35 0.0188 0.282 ns          
    ##  4 Shannon WT     L4         4     4     -2.3  0.0214 0.322 ns          
    ##  5 Shannon WT     L5         4     4     -1.5  0.134  1     ns          
    ##  6 Shannon L1     L2         4     4     -0.45 0.653  1     ns          
    ##  7 Shannon L1     L3         4     4     -0.85 0.395  1     ns          
    ##  8 Shannon L1     L4         4     4     -0.8  0.424  1     ns          
    ##  9 Shannon L1     L5         4     4      0    1      1     ns          
    ## 10 Shannon L2     L3         4     4     -0.4  0.689  1     ns          
    ## 11 Shannon L2     L4         4     4     -0.35 0.726  1     ns          
    ## 12 Shannon L2     L5         4     4      0.45 0.653  1     ns          
    ## 13 Shannon L3     L4         4     4      0.05 0.960  1     ns          
    ## 14 Shannon L3     L5         4     4      0.85 0.395  1     ns          
    ## 15 Shannon L4     L5         4     4      0.8  0.424  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-30-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df      p method        
    ## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Chao1    24      10.1     5 0.0726 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4    -1.28  0.201   1     ns          
    ##  2 Chao1 WT     L2         4     4    -2.71  0.00672 0.101 ns          
    ##  3 Chao1 WT     L3         4     4    -1.08  0.280   1     ns          
    ##  4 Chao1 WT     L4         4     4    -2.53  0.0112  0.169 ns          
    ##  5 Chao1 WT     L5         4     4    -1.73  0.0833  1     ns          
    ##  6 Chao1 L1     L2         4     4    -1.43  0.153   1     ns          
    ##  7 Chao1 L1     L3         4     4     0.201 0.841   1     ns          
    ##  8 Chao1 L1     L4         4     4    -1.25  0.210   1     ns          
    ##  9 Chao1 L1     L5         4     4    -0.452 0.651   1     ns          
    ## 10 Chao1 L2     L3         4     4     1.63  0.103   1     ns          
    ## 11 Chao1 L2     L4         4     4     0.176 0.861   1     ns          
    ## 12 Chao1 L2     L5         4     4     0.979 0.328   1     ns          
    ## 13 Chao1 L3     L4         4     4    -1.46  0.145   1     ns          
    ## 14 Chao1 L3     L5         4     4    -0.653 0.514   1     ns          
    ## 15 Chao1 L4     L5         4     4     0.803 0.422   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-30-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-30-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-30-5.png)<!-- -->

``` r
analysis_f(ps1_its_f_root_filtered)
```

    ##                R2 Pr(>F)   
    ## genotype  0.31265   0.01 **
    ## Residuals 0.68735          
    ## Total     1.00000          
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.1379457 0.7269169 0.1080609   0.666  0.6660000    
    ## 2  L1 vs L3  1 0.4037263 2.3654792 0.2827667   0.029  0.1025000    
    ## 3  L1 vs L4  1 0.2083547 1.0095451 0.1440243   0.463  0.5342308    
    ## 4  L1 vs L5  1 0.2534287 1.3994027 0.1891237   0.129  0.2150000    
    ## 5  L1 vs WT  1 0.4561952 2.7893054 0.3173522   0.021  0.1025000    
    ## 6  L2 vs L3  1 0.2394695 1.5302006 0.2032085   0.054  0.1157143    
    ## 7  L2 vs L4  1 0.1759052 0.9151912 0.1323450   0.555  0.5946429    
    ## 8  L2 vs L5  1 0.2302015 1.3791213 0.1868951   0.160  0.2400000    
    ## 9  L2 vs WT  1 0.3002886 2.0103295 0.2509671   0.024  0.1025000    
    ## 10 L3 vs L4  1 0.2855274 1.6493793 0.2156226   0.080  0.1500000    
    ## 11 L3 vs L5  1 0.3367680 2.2781525 0.2752006   0.036  0.1025000    
    ## 12 L3 vs WT  1 0.1423992 1.0930341 0.1540997   0.398  0.4975000    
    ## 13 L4 vs L5  1 0.2147014 1.1698085 0.1631576   0.281  0.3831818    
    ## 14 L4 vs WT  1 0.3594626 2.1655763 0.2652080   0.032  0.1025000    
    ## 15 L5 vs WT  1 0.3787817 2.6920770 0.3097162   0.041  0.1025000    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.013010 0.0026021 0.6811    999  0.635
    ## Residuals 18 0.068771 0.0038206                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##         L1      L2      L3      L4      L5    WT
    ## L1         0.98500 0.50000 0.82400 0.17900 0.613
    ## L2 0.98626         0.51000 0.85000 0.20300 0.658
    ## L3 0.50516 0.52749         0.58600 0.60900 0.737
    ## L4 0.82336 0.84742 0.60613         0.21600 0.769
    ## L5 0.15911 0.18247 0.59511 0.20435         0.283
    ## WT 0.61528 0.65068 0.74037 0.77835 0.27085

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df      p method        
    ## * <chr>   <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Shannon    24      10.2     5 0.0706 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      2.6  0.00932 0.140 ns          
    ##  2 Shannon WT     L2         4     4      1.55 0.121   1     ns          
    ##  3 Shannon WT     L3         4     4      0.65 0.516   1     ns          
    ##  4 Shannon WT     L4         4     4      2.4  0.0164  0.246 ns          
    ##  5 Shannon WT     L5         4     4      1.8  0.0719  1     ns          
    ##  6 Shannon L1     L2         4     4     -1.05 0.294   1     ns          
    ##  7 Shannon L1     L3         4     4     -1.95 0.0512  0.768 ns          
    ##  8 Shannon L1     L4         4     4     -0.2  0.841   1     ns          
    ##  9 Shannon L1     L5         4     4     -0.8  0.424   1     ns          
    ## 10 Shannon L2     L3         4     4     -0.9  0.368   1     ns          
    ## 11 Shannon L2     L4         4     4      0.85 0.395   1     ns          
    ## 12 Shannon L2     L5         4     4      0.25 0.803   1     ns          
    ## 13 Shannon L3     L4         4     4      1.75 0.0801  1     ns          
    ## 14 Shannon L3     L5         4     4      1.15 0.250   1     ns          
    ## 15 Shannon L4     L5         4     4     -0.6  0.549   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-31-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df      p method        
    ## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Chao1    24      9.24     5 0.0998 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     2.28  0.0227  0.340 ns          
    ##  2 Chao1 WT     L2         4     4     2.65  0.00795 0.119 ns          
    ##  3 Chao1 WT     L3         4     4     2.38  0.0174  0.261 ns          
    ##  4 Chao1 WT     L4         4     4     1.85  0.0639  0.959 ns          
    ##  5 Chao1 WT     L5         4     4     1.50  0.133   1     ns          
    ##  6 Chao1 L1     L2         4     4     0.376 0.707   1     ns          
    ##  7 Chao1 L1     L3         4     4     0.100 0.920   1     ns          
    ##  8 Chao1 L1     L4         4     4    -0.426 0.670   1     ns          
    ##  9 Chao1 L1     L5         4     4    -0.776 0.438   1     ns          
    ## 10 Chao1 L2     L3         4     4    -0.275 0.783   1     ns          
    ## 11 Chao1 L2     L4         4     4    -0.801 0.423   1     ns          
    ## 12 Chao1 L2     L5         4     4    -1.15  0.249   1     ns          
    ## 13 Chao1 L3     L4         4     4    -0.526 0.599   1     ns          
    ## 14 Chao1 L3     L5         4     4    -0.876 0.381   1     ns          
    ## 15 Chao1 L4     L5         4     4    -0.351 0.726   1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-31-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-31-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-31-5.png)<!-- -->

``` r
analysis_f(ps1_its_f_root_associated_soils_filtered)
```

    ##                R2 Pr(>F)  
    ## genotype  0.26359  0.065 .
    ## Residuals 0.73641         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.1758381 0.9645092 0.1384892   0.485  0.5989286    
    ## 2  L1 vs L3  1 0.2245259 1.2664647 0.1742890   0.175  0.5250000    
    ## 3  L1 vs L4  1 0.1620217 0.9720018 0.1394150   0.491  0.5989286    
    ## 4  L1 vs L5  1 0.1894568 1.1151814 0.1567327   0.376  0.5989286    
    ## 5  L1 vs WT  1 0.2829924 1.4353139 0.1930401   0.256  0.5850000    
    ## 6  L2 vs L3  1 0.2138286 1.0351063 0.1471344   0.378  0.5989286    
    ## 7  L2 vs L4  1 0.1866277 0.9522811 0.1369739   0.419  0.5989286    
    ## 8  L2 vs L5  1 0.1887214 0.9474928 0.1363791   0.559  0.5989286    
    ## 9  L2 vs WT  1 0.4118105 1.8185083 0.2325902   0.109  0.4087500    
    ## 10 L3 vs L4  1 0.2369345 1.2407753 0.1713594   0.273  0.5850000    
    ## 11 L3 vs L5  1 0.1794172 0.9240838 0.1334594   0.615  0.6150000    
    ## 12 L3 vs WT  1 0.4113402 1.8576341 0.2364114   0.106  0.4087500    
    ## 13 L4 vs L5  1 0.1719940 0.9369906 0.1350716   0.523  0.5989286    
    ## 14 L4 vs WT  1 0.4042691 1.9174629 0.2421815   0.032  0.4087500    
    ## 15 L5 vs WT  1 0.3440800 1.6075835 0.2113133   0.090  0.4087500    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.015358 0.0030715 0.9163    999  0.485
    ## Residuals 18 0.060340 0.0033522                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.8600000 0.5900000 0.5120000 0.4910000 0.012
    ## L2 0.8747768           0.6630000 0.6420000 0.5880000 0.054
    ## L3 0.5878485 0.6693774           0.9340000 0.9590000 0.363
    ## L4 0.5284959 0.6573042 0.9302325           0.9660000 0.182
    ## L5 0.4844019 0.6147613 0.9630370 0.9583715           0.189
    ## WT 0.0074637 0.0397309 0.3478530 0.1712188 0.1864584

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      8.05     5 0.154 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic      p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>  <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4      1.05 0.294  1     ns          
    ##  2 Shannon WT     L2         4     4      1.5  0.134  1     ns          
    ##  3 Shannon WT     L3         4     4      1.7  0.0891 1     ns          
    ##  4 Shannon WT     L4         4     4      2.4  0.0164 0.246 ns          
    ##  5 Shannon WT     L5         4     4      2.35 0.0188 0.282 ns          
    ##  6 Shannon L1     L2         4     4      0.45 0.653  1     ns          
    ##  7 Shannon L1     L3         4     4      0.65 0.516  1     ns          
    ##  8 Shannon L1     L4         4     4      1.35 0.177  1     ns          
    ##  9 Shannon L1     L5         4     4      1.3  0.194  1     ns          
    ## 10 Shannon L2     L3         4     4      0.2  0.841  1     ns          
    ## 11 Shannon L2     L4         4     4      0.9  0.368  1     ns          
    ## 12 Shannon L2     L5         4     4      0.85 0.395  1     ns          
    ## 13 Shannon L3     L4         4     4      0.7  0.484  1     ns          
    ## 14 Shannon L3     L5         4     4      0.65 0.516  1     ns          
    ## 15 Shannon L4     L5         4     4     -0.05 0.960  1     ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df      p method        
    ## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Chao1    24      9.96     5 0.0764 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p  p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl>  <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     1.90  0.0573  0.860  ns          
    ##  2 Chao1 WT     L2         4     4     2.90  0.00372 0.0557 ns          
    ##  3 Chao1 WT     L3         4     4     2.20  0.0277  0.416  ns          
    ##  4 Chao1 WT     L4         4     4     1.70  0.0890  1      ns          
    ##  5 Chao1 WT     L5         4     4     2.40  0.0163  0.245  ns          
    ##  6 Chao1 L1     L2         4     4     1.00  0.317   1      ns          
    ##  7 Chao1 L1     L3         4     4     0.300 0.764   1      ns          
    ##  8 Chao1 L1     L4         4     4    -0.200 0.841   1      ns          
    ##  9 Chao1 L1     L5         4     4     0.500 0.617   1      ns          
    ## 10 Chao1 L2     L3         4     4    -0.700 0.484   1      ns          
    ## 11 Chao1 L2     L4         4     4    -1.20  0.230   1      ns          
    ## 12 Chao1 L2     L5         4     4    -0.500 0.617   1      ns          
    ## 13 Chao1 L3     L4         4     4    -0.500 0.617   1      ns          
    ## 14 Chao1 L3     L5         4     4     0.200 0.841   1      ns          
    ## 15 Chao1 L4     L5         4     4     0.700 0.484   1      ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-32-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-32-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-32-5.png)<!-- -->

``` r
analysis_f(ps1_its_f_soil_filtered)
```

    ##                R2 Pr(>F)
    ## genotype  0.20214  0.838
    ## Residuals 0.79786       
    ## Total     1.00000       
    ##       pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.3371861 0.8399887 0.1228056   0.732  0.9621429    
    ## 2  L1 vs L3  1 0.3678616 0.9780774 0.1401643   0.437  0.9621429    
    ## 3  L1 vs L4  1 0.3732895 0.9434561 0.1358770   0.681  0.9621429    
    ## 4  L1 vs L5  1 0.3952730 1.0259765 0.1460262   0.431  0.9621429    
    ## 5  L1 vs WT  1 0.3613888 0.8874086 0.1288451   0.747  0.9621429    
    ## 6  L2 vs L3  1 0.3841450 0.9306005 0.1342742   0.474  0.9621429    
    ## 7  L2 vs L4  1 0.4037676 0.9338960 0.1346856   0.570  0.9621429    
    ## 8  L2 vs L5  1 0.3673264 0.8705429 0.1267066   0.693  0.9621429    
    ## 9  L2 vs WT  1 0.3516152 0.7920575 0.1166153   1.000  1.0000000    
    ## 10 L3 vs L4  1 0.4454128 1.0942810 0.1542483   0.321  0.9621429    
    ## 11 L3 vs L5  1 0.3562795 0.8982433 0.1302133   0.766  0.9621429    
    ## 12 L3 vs WT  1 0.4092023 0.9775128 0.1400947   0.510  0.9621429    
    ## 13 L4 vs L5  1 0.3436729 0.8257493 0.1209756   0.840  0.9621429    
    ## 14 L4 vs WT  1 0.3829364 0.8739435 0.1271386   0.898  0.9621429    
    ## 15 L5 vs WT  1 0.3583699 0.8377550 0.1225190   0.889  0.9621429    
    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
    ## Groups     5 0.003892 0.0007785 0.1395    999  0.976
    ## Residuals 18 0.100474 0.0055819                     
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##         L1      L2      L3      L4      L5    WT
    ## L1         0.14900 0.72000 0.76200 0.62300 0.472
    ## L2 0.14487         0.65100 0.79300 0.66100 0.956
    ## L3 0.74567 0.62818         0.95900 0.95000 0.757
    ## L4 0.76658 0.78269 0.95010         0.99600 0.809
    ## L5 0.63593 0.65061 0.94397 0.99106         0.742
    ## WT 0.48298 0.95666 0.70852 0.79805 0.73780

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df     p method        
    ## * <chr>   <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Shannon    24      2.80     5 0.731 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     -0.55 0.582     1 ns          
    ##  2 Shannon WT     L2         4     4      0.1  0.920     1 ns          
    ##  3 Shannon WT     L3         4     4     -0.35 0.726     1 ns          
    ##  4 Shannon WT     L4         4     4      0.55 0.582     1 ns          
    ##  5 Shannon WT     L5         4     4      0.85 0.395     1 ns          
    ##  6 Shannon L1     L2         4     4      0.65 0.516     1 ns          
    ##  7 Shannon L1     L3         4     4      0.2  0.841     1 ns          
    ##  8 Shannon L1     L4         4     4      1.1  0.271     1 ns          
    ##  9 Shannon L1     L5         4     4      1.4  0.162     1 ns          
    ## 10 Shannon L2     L3         4     4     -0.45 0.653     1 ns          
    ## 11 Shannon L2     L4         4     4      0.45 0.653     1 ns          
    ## 12 Shannon L2     L5         4     4      0.75 0.453     1 ns          
    ## 13 Shannon L3     L4         4     4      0.9  0.368     1 ns          
    ## 14 Shannon L3     L5         4     4      1.2  0.230     1 ns          
    ## 15 Shannon L4     L5         4     4      0.3  0.764     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-33-2.png)<!-- -->

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df     p method        
    ## * <chr> <int>     <dbl> <int> <dbl> <chr>         
    ## 1 Chao1    24      1.93     5 0.859 Kruskal-Wallis
    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic     p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4   -0.450  0.652     1 ns          
    ##  2 Chao1 WT     L2         4     4    0.200  0.841     1 ns          
    ##  3 Chao1 WT     L3         4     4   -0.0250 0.980     1 ns          
    ##  4 Chao1 WT     L4         4     4    0.776  0.438     1 ns          
    ##  5 Chao1 WT     L5         4     4    0.550  0.582     1 ns          
    ##  6 Chao1 L1     L2         4     4    0.650  0.515     1 ns          
    ##  7 Chao1 L1     L3         4     4    0.425  0.671     1 ns          
    ##  8 Chao1 L1     L4         4     4    1.23   0.220     1 ns          
    ##  9 Chao1 L1     L5         4     4    1.00   0.317     1 ns          
    ## 10 Chao1 L2     L3         4     4   -0.225  0.822     1 ns          
    ## 11 Chao1 L2     L4         4     4    0.575  0.565     1 ns          
    ## 12 Chao1 L2     L5         4     4    0.350  0.726     1 ns          
    ## 13 Chao1 L3     L4         4     4    0.801  0.423     1 ns          
    ## 14 Chao1 L3     L5         4     4    0.575  0.565     1 ns          
    ## 15 Chao1 L4     L5         4     4   -0.225  0.822     1 ns

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-33-3.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-33-4.png)<!-- -->![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-33-5.png)<!-- -->

# Example

``` r
# 1. Remove data points with missing metadata
phy_na_ps1_ssu_g_leaf_filtered <- ps1_ssu_g_leaf_filtered %>%
  subset_samples(
      !is.na(genotype) &
      !is.na(type) &
      !is.na(loc.time) &
      !is.na(id))

# 2. Differences depending on input variables
ps1_bray_ps1_ssu_g_leaf_filtered <- phyloseq::distance(physeq = phy_na_ps1_ssu_g_leaf_filtered, method = "bray")
cap_ord_ps1_ssu_g_leaf_filtered <- ordinate(physeq = phy_na_ps1_ssu_g_leaf_filtered, method = "CAP", 
                        distance = ps1_bray_ps1_ssu_g_leaf_filtered,
                        formula = ~ genotype)

# 3. ANOVA
obj_phy_na_ps1_ssu_g_leaf_filtered <- t(phy_na_ps1_ssu_g_leaf_filtered)
anova(cap_ord_ps1_ssu_g_leaf_filtered, by="terms")
```

    ## Permutation test for capscale under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: capscale(formula = distance ~ genotype, data = data)
    ##          Df SumOfSqs     F Pr(>F)  
    ## genotype  5  0.46043 2.203   0.04 *
    ## Residual 17  0.71059               
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
vif.cca(cap_ord_ps1_ssu_g_leaf_filtered)
```

    ## genotypeL2 genotypeL3 genotypeL4 genotypeL5 genotypeWT 
    ##   1.652174   1.652174   1.652174   1.521739   1.652174

``` r
# 4. CAP plot
cap_plot_ps1_ssu_g_leaf_filtered <- plot_ordination(physeq = phy_na_ps1_ssu_g_leaf_filtered, 
                                ordination = cap_ord_ps1_ssu_g_leaf_filtered,
                                color = "genotype", 
                                axes = c(1,2)) + xlim(-2.3, 2.2) + ylim(-2.2, 2.5) + 
                                geom_point(aes(colour = genotype), 
                                           alpha = 0.8, size = 3) +
                                geom_point(size = 3) + 
                                annotate("text", x = -0.8, y = 2.5, label = bquote("R[PERMANOVA]^2 == 0.015"), 
                                         parse = TRUE) +
                                annotate("text", x = 0.9, y = 2.48, label = bquote("p[PERMANOVA] == 0.018"), 
                                         parse = TRUE) +
                                scale_color_manual(
                                  values = c("#00AFBB","#FC4E07", "#E7B800",
                                             "#36BB07", "#7D26CD",
                                             "#757575")) + theme_bw()
cap_plot_ps1_ssu_g_leaf_filtered_2 <- cap_plot_ps1_ssu_g_leaf_filtered
cap_plot_ps1_ssu_g_leaf_filtered_2$data$genotype<-factor(cap_plot_ps1_ssu_g_leaf_filtered$data$genotype, levels = c("WT", "L1", "L2", "L3", "L4", "L5"))
print(cap_plot_ps1_ssu_g_leaf_filtered_2)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_point).

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
# 5. Create distance matrix for PERMANOVA
dist_phy_na_ps1_ssu_g_leaf_filtered <-vegdist(decostand(obj_phy_na_ps1_ssu_g_leaf_filtered@otu_table,"hellinger"),"bray")

# 6. transform in data-frame
sampledf_dist_phy_na_ps1_ssu_g_leaf_filtered <- data.frame(sample_data(obj_phy_na_ps1_ssu_g_leaf_filtered))

# 7. PERMANOVA (pairwise_adonis)
adonis_table_ps1_ssu_g_leaf_filtered <- adonis(obj_phy_na_ps1_ssu_g_leaf_filtered@otu_table ~ genotype, data = sampledf_dist_phy_na_ps1_ssu_g_leaf_filtered, permutations = 999)
adonis_result_ps1_ssu_g_leaf_filtered <- subset(adonis_table_ps1_ssu_g_leaf_filtered$aov.tab, select = c("R2", "Pr(>F)"))
adonis_result_ps1_ssu_g_leaf_filtered
```

    ##                R2 Pr(>F)  
    ## genotype  0.39651  0.054 .
    ## Residuals 0.60349         
    ## Total     1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
pairwise.adonis(obj_phy_na_ps1_ssu_g_leaf_filtered@otu_table,
                sample_data(obj_phy_na_ps1_ssu_g_leaf_filtered)$genotype, p.adjust.m = "fdr",
                perm = 999)
```

    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.
    ## Set of permutations < 'minperm'. Generating entire set.

    ##       pairs Df   SumsOfSqs    F.Model         R2 p.value p.adjusted sig
    ## 1  L1 vs L2  1 0.095315868  3.6497420 0.37822172   0.124  0.2833333    
    ## 2  L1 vs L3  1 0.074089124  1.7736549 0.22816229   0.168  0.2833333    
    ## 3  L1 vs L4  1 0.059194073  1.0142823 0.14460243   0.419  0.5713636    
    ## 4  L1 vs L5  1 0.076725734  3.0225880 0.37675972   0.134  0.2833333    
    ## 5  L1 vs WT  1 0.003204112  0.1717821 0.02783347   0.821  0.8796429    
    ## 6  L2 vs L3  1 0.022107672  0.5192182 0.07964424   0.579  0.7237500    
    ## 7  L2 vs L4  1 0.017093881  0.2889075 0.04593922   0.908  0.9080000    
    ## 8  L2 vs L5  1 0.283415991 10.7549132 0.68263868   0.022  0.1300000    
    ## 9  L2 vs WT  1 0.102787831  5.2822918 0.46819316   0.022  0.1300000    
    ## 10 L3 vs L4  1 0.024957043  0.3335453 0.05266329   0.639  0.7373077    
    ## 11 L3 vs L5  1 0.220728786  4.8899032 0.49443388   0.026  0.1300000    
    ## 12 L3 vs WT  1 0.082417884  2.3470734 0.28118519   0.170  0.2833333    
    ## 13 L4 vs L5  1 0.209542772  3.2214573 0.39183531   0.086  0.2790000    
    ## 14 L4 vs WT  1 0.067871806  1.3127066 0.17951036   0.201  0.3015000    
    ## 15 L5 vs WT  1 0.077339112  4.4458202 0.47066534   0.093  0.2790000

``` r
# 8. Betadiversity / dispersion / homogeneity between groups
bet_factor_ps_ps1_ssu_g_leaf_filtered <-betadisper(dist_phy_na_ps1_ssu_g_leaf_filtered,
                               sampledf_dist_phy_na_ps1_ssu_g_leaf_filtered$genotype)
permutest(bet_factor_ps_ps1_ssu_g_leaf_filtered, pairwise = TRUE, permutations = 999,
          group=bet_factor_ps_ps1_ssu_g_leaf_filtered$group)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
    ## Groups     5 0.038102 0.0076204 2.8963    999  0.039 *
    ## Residuals 17 0.044729 0.0026311                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##           L1        L2        L3        L4        L5    WT
    ## L1           0.3660000 0.4430000 0.2670000 0.0820000 0.192
    ## L2 0.3859514           0.8810000 0.5950000 0.1400000 0.026
    ## L3 0.4592197 0.8712403           0.5250000 0.1220000 0.045
    ## L4 0.2731515 0.6075788 0.5327854           0.5740000 0.055
    ## L5 0.0749717 0.1641981 0.1306237 0.6020433           0.006
    ## WT 0.2046310 0.0248829 0.0316000 0.0463204 0.0032668

``` r
# 4. Alpha diversity plot
GP_ps1_ssu_g_leaf_filtered <- prune_taxa(taxa_sums(ps1_ssu_g_leaf_filtered) > 0, ps1_ssu_g_leaf_filtered)
p_ps1_ssu_g_leaf_filtered<-plot_richness(GP_ps1_ssu_g_leaf_filtered, x="genotype",measures=c("Chao1", "Shannon"), color = "genotype") + theme_bw()
```

    ## Warning in estimate_richness(physeq, split = TRUE, measures = measures): The data you have provided does not have
    ## any singletons. This is highly suspicious. Results of richness
    ## estimates (for example) are probably unreliable, or wrong, if you have already
    ## trimmed low-abundance taxa from the data.
    ## 
    ## We recommended that you find the un-trimmed data and retry.

``` r
p_ps1_ssu_g_leaf_filtered + geom_boxplot(data = p_ps1_ssu_g_leaf_filtered$data, aes(x = genotype, y = value, color = NULL), alpha = 0.1)
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-2.png)<!-- -->

``` r
#write.csv(p_ssu$plot_env$DF, file="alpha_diversity_ssu.csv")

p_ps1_ssu_g_leaf_filtered_data <- p_ps1_ssu_g_leaf_filtered$plot_env$DF
ordered_p_ps1_ssu_g_leaf_filtered_data <- reorder_levels(p_ps1_ssu_g_leaf_filtered_data, genotype, c("WT", "L1", "L2", "L3", "L4", "L5"))

res.kruskal_total_shannon_ps1_ssu_g_leaf_filtered <- ordered_p_ps1_ssu_g_leaf_filtered_data %>% 
  kruskal_test(Shannon ~ genotype)
res.kruskal_total_shannon_ps1_ssu_g_leaf_filtered
```

    ## # A tibble: 1 × 6
    ##   .y.         n statistic    df       p method        
    ## * <chr>   <int>     <dbl> <int>   <dbl> <chr>         
    ## 1 Shannon    23      15.5     5 0.00857 Kruskal-Wallis

``` r
pwc_total_shannon_ps1_ssu_g_leaf_filtered <- ordered_p_ps1_ssu_g_leaf_filtered_data %>% 
  dunn_test(Shannon ~ genotype, p.adjust.method = "bonferroni") 
pwc_total_shannon_ps1_ssu_g_leaf_filtered
```

    ## # A tibble: 15 × 9
    ##    .y.     group1 group2    n1    n2 statistic        p  p.adj p.adj.signif
    ##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>    <dbl>  <dbl> <chr>       
    ##  1 Shannon WT     L1         4     4     1.20  0.231    1      ns          
    ##  2 Shannon WT     L2         4     4     1.62  0.106    1      ns          
    ##  3 Shannon WT     L3         4     4     2.81  0.00488  0.0732 ns          
    ##  4 Shannon WT     L4         4     4     0.938 0.348    1      ns          
    ##  5 Shannon WT     L5         4     3     3.36  0.000773 0.0116 *           
    ##  6 Shannon L1     L2         4     4     0.417 0.677    1      ns          
    ##  7 Shannon L1     L3         4     4     1.62  0.106    1      ns          
    ##  8 Shannon L1     L4         4     4    -0.261 0.794    1      ns          
    ##  9 Shannon L1     L5         4     3     2.25  0.0243   0.365  ns          
    ## 10 Shannon L2     L3         4     4     1.20  0.231    1      ns          
    ## 11 Shannon L2     L4         4     4    -0.678 0.498    1      ns          
    ## 12 Shannon L2     L5         4     3     1.87  0.0620   0.930  ns          
    ## 13 Shannon L3     L4         4     4    -1.88  0.0606   0.909  ns          
    ## 14 Shannon L3     L5         4     3     0.756 0.450    1      ns          
    ## 15 Shannon L4     L5         4     3     2.49  0.0126   0.190  ns

``` r
pwc_total_shannon_ps1_ssu_g_leaf_filtered <- pwc_total_shannon_ps1_ssu_g_leaf_filtered %>% add_xy_position()
ggboxplot(ordered_p_ps1_ssu_g_leaf_filtered_data, x = "genotype", y = "Shannon", color = "genotype", 
          palette = c("#00AFBB","#FC4E07", "#E7B800", "#36BB07", "#7D26CD", "#757575"),
          order = c("WT", "L1", "L2", "L3", "L4", "L5"),
          ylab = "Shannon index", xlab = "genotype") +
          stat_pvalue_manual(pwc_total_shannon_ps1_ssu_g_leaf_filtered, y.position = 8.0,
          step.increase = 0.1, hide.ns = TRUE) + 
          labs(subtitle = get_test_label(res.kruskal_total_shannon_ps1_ssu_g_leaf_filtered, detailed = TRUE),
          caption = get_pwc_label(pwc_total_shannon_ps1_ssu_g_leaf_filtered)) + coord_cartesian(ylim = c(0, 15))
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-3.png)<!-- -->

``` r
res.kruskal_total_chao1_ps1_ssu_g_leaf_filtered <- ordered_p_ps1_ssu_g_leaf_filtered_data %>% 
  kruskal_test(Chao1 ~ genotype)
res.kruskal_total_chao1_ps1_ssu_g_leaf_filtered
```

    ## # A tibble: 1 × 6
    ##   .y.       n statistic    df      p method        
    ## * <chr> <int>     <dbl> <int>  <dbl> <chr>         
    ## 1 Chao1    23      9.46     5 0.0921 Kruskal-Wallis

``` r
pwc_total_chao1_ps1_ssu_g_leaf_filtered <- ordered_p_ps1_ssu_g_leaf_filtered_data %>% 
  dunn_test(Chao1 ~ genotype, p.adjust.method = "bonferroni") 
pwc_total_chao1_ps1_ssu_g_leaf_filtered
```

    ## # A tibble: 15 × 9
    ##    .y.   group1 group2    n1    n2 statistic       p p.adj p.adj.signif
    ##  * <chr> <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl> <chr>       
    ##  1 Chao1 WT     L1         4     4     1.02  0.309   1     ns          
    ##  2 Chao1 WT     L2         4     4     0.339 0.734   1     ns          
    ##  3 Chao1 WT     L3         4     4     1.20  0.230   1     ns          
    ##  4 Chao1 WT     L4         4     4    -0.313 0.754   1     ns          
    ##  5 Chao1 WT     L5         4     3     2.42  0.0157  0.235 ns          
    ##  6 Chao1 L1     L2         4     4    -0.679 0.497   1     ns          
    ##  7 Chao1 L1     L3         4     4     0.183 0.855   1     ns          
    ##  8 Chao1 L1     L4         4     4    -1.33  0.183   1     ns          
    ##  9 Chao1 L1     L5         4     3     1.47  0.140   1     ns          
    ## 10 Chao1 L2     L3         4     4     0.861 0.389   1     ns          
    ## 11 Chao1 L2     L4         4     4    -0.653 0.514   1     ns          
    ## 12 Chao1 L2     L5         4     3     2.10  0.0355  0.533 ns          
    ## 13 Chao1 L3     L4         4     4    -1.51  0.130   1     ns          
    ## 14 Chao1 L3     L5         4     3     1.31  0.192   1     ns          
    ## 15 Chao1 L4     L5         4     3     2.71  0.00680 0.102 ns

``` r
pwc_total_chao1_ps1_ssu_g_leaf_filtered <- pwc_total_chao1_ps1_ssu_g_leaf_filtered %>% add_xy_position()
ggboxplot(ordered_p_ps1_ssu_g_leaf_filtered_data, x = "genotype", y = "Chao1", color = "genotype", 
          palette = c("#00AFBB","#FC4E07", "#E7B800", "#36BB07", "#7D26CD", "#757575"),
          order = c("WT", "L1", "L2", "L3", "L4", "L5"),
          ylab = "chao1 index", xlab = "genotype") +
          stat_pvalue_manual(pwc_total_chao1_ps1_ssu_g_leaf_filtered, y.position = 2300,
          step.increase = 0.1, hide.ns = TRUE) + 
          labs(subtitle = get_test_label(res.kruskal_total_chao1_ps1_ssu_g_leaf_filtered, detailed = TRUE),
          caption = get_pwc_label(pwc_total_chao1_ps1_ssu_g_leaf_filtered)) + coord_cartesian(ylim = c(0, 4500))
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-4.png)<!-- -->

``` r
ps_normalization_ps1_ssu_g_leaf_filtered = transform_sample_counts(ps1_ssu_g_leaf_filtered, function(x) x/sum(x))
phyloseq::plot_bar(ps_normalization_ps1_ssu_g_leaf_filtered, fill = "Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ genotype, scales = "free") + theme(panel.background = element_blank() +
  theme(legend.position = "bottom"), axis.text.x=element_blank(), axis.ticks.x=element_blank())
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-5.png)<!-- -->

``` r
#venn diagramming
wt_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "WT") %>% filter_taxa(function(x) sum(x) > 0, T)
l1_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "L1") %>% filter_taxa(function(x) sum(x) > 0, T)
l2_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "L2") %>% filter_taxa(function(x) sum(x) > 0, T)
l3_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "L3") %>% filter_taxa(function(x) sum(x) > 0, T)
l4_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "L4") %>% filter_taxa(function(x) sum(x) > 0, T)
l5_sub_ssu <- subset_samples(ps1_ssu_g_leaf_filtered, genotype == "L5") %>% filter_taxa(function(x) sum(x) > 0, T)
wt_asvs_ssu <- taxa_names(wt_sub_ssu)
l1_asvs_ssu <- taxa_names(l1_sub_ssu)
l2_asvs_ssu <- taxa_names(l2_sub_ssu)
l3_asvs_ssu <- taxa_names(l3_sub_ssu)
l4_asvs_ssu <- taxa_names(l4_sub_ssu)
l5_asvs_ssu <- taxa_names(l5_sub_ssu)
genotype_venn_ssu <- venn(list("WT" = wt_asvs_ssu, "L1" = l1_asvs_ssu, "L2" = l2_asvs_ssu, "L3" = l3_asvs_ssu, "L4" = l4_asvs_ssu, "L5"= l5_asvs_ssu))
```

![](oilcane_analysis-20211201-_files/figure-gfm/unnamed-chunk-34-6.png)<!-- -->

``` r
#write.csv(l1_asvs_ssu,file="l1_asvs_ssu.csv")
#write.csv(l2_asvs_ssu,file="l2_asvs_ssu.csv")
#write.csv(l3_asvs_ssu,file="l3_asvs_ssu.csv")
#write.csv(l4_asvs_ssu,file="l4_asvs_ssu.csv")
#write.csv(l5_asvs_ssu,file="l5_asvs_ssu.csv")
#write.csv(wt_asvs_ssu,file="wt_asvs_ssu.csv")
```
