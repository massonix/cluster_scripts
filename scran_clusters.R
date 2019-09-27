# Load packages
library(scater)
library(scran)
library(Seurat)
library(ggpubr)
library(tidyverse)

# Load data
sce <- readRDS("../data/R_objects/cll_sce_3_donors.rds")

# Normalize
clusters <- quickCluster(
  sce, 
  use.ranks = FALSE, 
  block = factor(sce$donor), 
  BPPARAM = MulticoreParam(4),
  block.BPPARAM = MulticoreParam(4)
)
sce <- computeSumFactors(sce, cluster = clusters, BPPARAM =  MulticoreParam(4))
summary(sizeFactors(sce))
sce <- normalize(sce)

# Save
saveRDS(sce, "../data/R_objects/cll_sce_3_donors_normalized.rds")