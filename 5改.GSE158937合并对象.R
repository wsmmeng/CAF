library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
#####自动读取cellranger(LINUX)输出的feature barcode matric
rm(list = ls())
OVH1 <- Read10X(data.dir = "GSM4816045/")
OVH2 <- Read10X(data.dir = "GSM4816046/")
OVH3 <- Read10X(data.dir = "GSM4816047/")
#自动读取10X的数据，是一些tsv与mtx文件
OVH1 <- CreateSeuratObject(counts = OVH1)
OVH2 <- CreateSeuratObject(counts = OVH2)
OVH3 <- CreateSeuratObject(counts = OVH3)

OVH1 <- as.data.frame(OVH1[["RNA"]]@counts)
OVH2 <- as.data.frame(OVH2[["RNA"]]@counts)
OVH3 <- as.data.frame(OVH3[["RNA"]]@counts)

##不用intersect
colnames(OVH1) <- paste(colnames(OVH1),'H1',sep = '-')
colnames(OVH2) <- paste(colnames(OVH2),'H2',sep = '-')
colnames(OVH3) <- paste(colnames(OVH3),'H3',sep = '-')
OVH <- CreateSeuratObject(cbind(OVH1,OVH2,OVH3),project = 'OVH',min.cells =5)
saveRDS(OVH,'OVH.rds')
