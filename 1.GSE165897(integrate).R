rm(list=ls())
##合并对象--数据整合成一个seurat对象，看看需不需要harmony
{library(multtest)
  if(!require(multtest))install.packages("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(mindr))install.packages("mindr")
  if(!require(mindr))install.packages("tidyverse")
  }

#读取数据
memory.limit(10000000)
library(data.table)

data <- fread("GSE165897_UMIcounts_HGSOC.tsv",sep = "\t",header = T)
dim(data)
data <- data.frame(data)
rownames(data) <- data[,1]
data <- data[,-1]
cellinfo <- fread("GSE165897_cellInfo_HGSOC.tsv",sep = "\t",header = T)
table(cellinfo$sample)
table(cellinfo$cell_type)
EOC1005I <- cellinfo[which(cellinfo$sample=='EOC1005_interval_Tumor'),]
EOC1005P <- cellinfo[which(cellinfo$sample=='EOC1005_primary_Peritoneum'),]
EOC136I <- cellinfo[which(cellinfo$sample=='EOC136_interval_Omentum'),]
EOC136P <- cellinfo[which(cellinfo$sample=='EOC136_primary_Mesentery'),]
EOC153I <- cellinfo[which(cellinfo$sample=='EOC153_interval_Omentum'),]
EOC153P <- cellinfo[which(cellinfo$sample=='EOC153_primary_Omentum'),]
EOC227I <- cellinfo[which(cellinfo$sample=='EOC227_interval_Omentum'),]
EOC227P <- cellinfo[which(cellinfo$sample=='EOC227_primary_Omentum'),]
EOC3I <- cellinfo[which(cellinfo$sample=='EOC3_interval_Omentum'),]
EOC3P <- cellinfo[which(cellinfo$sample=='EOC3_primary_Peritoneum'),]
EOC349I <- cellinfo[which(cellinfo$sample=='EOC349_interval_Omentum'),]
EOC349P <- cellinfo[which(cellinfo$sample=='EOC349_primary_Peritoneum'),]
EOC372I <- cellinfo[which(cellinfo$sample=='EOC372_interval_Peritoneum'),]
EOC372P <- cellinfo[which(cellinfo$sample=='EOC372_primary_Peritoneum'),]
EOC443I <- cellinfo[which(cellinfo$sample=='EOC443_interval_Omentum'),]
EOC443P <- cellinfo[which(cellinfo$sample=='EOC443_primary_Omentum'),]
EOC540I <- cellinfo[which(cellinfo$sample=='EOC540_interval_Omentum'),]
EOC540P <- cellinfo[which(cellinfo$sample=='EOC540_primary_Omentum'),]
EOC733I <- cellinfo[which(cellinfo$sample=='EOC733_interval_Omentum'),]
EOC733P <- cellinfo[which(cellinfo$sample=='EOC733_primary_Peritoneum'),]
EOC87I <- cellinfo[which(cellinfo$sample=='EOC87_interval_Omentum'),]
EOC87P <- cellinfo[which(cellinfo$sample=='EOC87_primary_Peritoneum'),]

##data
colnames(data)[1:20]
colnames(data) <- gsub('\\.','-',colnames(data))
count1005I <- data[,EOC1005I$cell] 
count1005P <- data[,EOC1005P$cell]
count136I <- data[,EOC136I$cell]
count136P <- data[,EOC136P$cell]
count153I <- data[,EOC153I$cell]
count153P <- data[,EOC153P$cell]
count227I <- data[,EOC227I$cell]
count227P <- data[,EOC227P$cell]
count349I <- data[,EOC349I$cell]
count349P <- data[,EOC349P$cell]
count372I <- data[,EOC372I$cell]
count372P <- data[,EOC372P$cell]
count3I <- data[,EOC3I$cell]
count3P <- data[,EOC3P$cell]
count443I <- data[,EOC443I$cell]
count443P <- data[,EOC443P$cell]
count540I <- data[,EOC540I$cell]
count540P <- data[,EOC540P$cell]
count733I <- data[,EOC733I$cell]
count733P <- data[,EOC733P$cell]
count87I <- data[,EOC87I$cell]
count87P <- data[,EOC87P$cell]
##构建seurat对象
S1005I <- CreateSeuratObject(counts=count1005I)
S1005P <- CreateSeuratObject(counts=count1005P)
S136I <- CreateSeuratObject(counts=count136I)
S136P <- CreateSeuratObject(counts=count136P)
S153I <- CreateSeuratObject(counts=count153I)
S153P <- CreateSeuratObject(counts=count153P)
S227I <- CreateSeuratObject(counts=count227I)
S227P <- CreateSeuratObject(counts=count227P)
S349I <- CreateSeuratObject(counts=count349I)
S349P <- CreateSeuratObject(counts=count349P)
S372I <- CreateSeuratObject(counts=count372I)
S372P <- CreateSeuratObject(counts=count372P)
S3I <- CreateSeuratObject(counts=count3I)
S3P <- CreateSeuratObject(counts=count3P)
S443I <- CreateSeuratObject(counts=count443I)
S443P <- CreateSeuratObject(counts=count443P)
S540I <- CreateSeuratObject(counts=count540I)
S540P <- CreateSeuratObject(counts=count540P)
S733I <- CreateSeuratObject(counts=count733I)
S733P <- CreateSeuratObject(counts=count733P)
S87I <- CreateSeuratObject(counts=count87I)
S87P <- CreateSeuratObject(counts=count87P)

saveRDS(S1005I,'S1005I.rds')
saveRDS(S1005P,'S1005P.rds')
saveRDS(S136I,'S136I.rds')     
saveRDS(S136P,'S136P.rds')  
saveRDS(S153I,'S153I.rds')     
saveRDS(S153P,'S153P.rds')     
saveRDS(S227I,'S227I.rds')     
saveRDS(S227P,'S227P.rds')
saveRDS(S349I,'S349I.rds')     
saveRDS(S349P,'S349P.rds')
saveRDS(S372I,'S372I.rds')     
saveRDS(S372P,'S372P.rds')
saveRDS(S3I,'S3I.rds')     
saveRDS(S3P,'S3P.rds')
saveRDS(S443I,'S443I.rds')     
saveRDS(S443P,'S443P.rds')
saveRDS(S540I,'S540I.rds')     
saveRDS(S540P,'S540P.rds')
saveRDS(S733I,'S733I.rds')     
saveRDS(S733P,'S733P.rds')
saveRDS(S87I,'S87I.rds')     
saveRDS(S87P,'S87P.rds')
