# Load library
library(devtools)
library(vegan)
library(gplots)
library(pheatmap)
library(ggplot2)
library(pvclust)

# Setup working directory
require("knitr")
opts_knit$set(root.dir = "~/Desktop/R_analysis/Oilcane_data/")


# Read TAG yield
TAG <- as.matrix(read.csv(file="TAG_yield.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1))
TAG_a <- t(TAG)

# Create heatmap
heatmap.2(TAG_a, scale = "row", Rowv = NA)

# Statistical analysis
pvclust(TAG_a, method.hclust = "ward.D2",
        method.dist = "euclidean", nboot = 999)
res.pv <- pvclust(TAG_a, method.dist="euclidean", 
                  method.hclust="ward.D2", nboot = 999)
plot(res.pv, hang = -1, cex = 0.5)
pvrect(res.pv)
clusters <- pvpick(res.pv)
clusters
