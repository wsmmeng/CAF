rm(list=ls())
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
library(ggplot2)
library(patchwork)

stromal<- readRDS('stromal.rds')
Idents(stromal) <- 'celltype'
stromalmarkers <- FindAllMarkers(stromal,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
celltypemarkers <- stromalmarkers
save(celltypemarkers,file='celltypemarkers.Rdata')
