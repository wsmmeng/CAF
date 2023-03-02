rm(list=ls())
OV <- readRDS('OV.rds')
#取细胞子集
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
##提取细胞子集
Cells.sub <- subset(OV@meta.data, celltype=="stromal cell")
stromal <- subset(OV, cells=row.names(Cells.sub))
saveRDS(stromal,'stromal.rds')





