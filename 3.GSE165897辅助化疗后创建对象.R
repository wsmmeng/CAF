rm(list=ls())
{library(multtest)
  if(!require(multtest))install.packages("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(mindr))install.packages("mindr")
  if(!require(mindr))install.packages("tidyverse")}

##读取数据
{S1005I <- readRDS('S1005I.rds')
  S1005P <- readRDS('S1005P.rds')
  S136I <- readRDS('S136I.rds')     
  S136P <- readRDS('S136P.rds')  
  S153I <- readRDS('S153I.rds')     
  S153P <- readRDS('S153P.rds')     
  S227I <- readRDS('S227I.rds')     
  S227P <- readRDS('S227P.rds')
  S349I <- readRDS('S349I.rds')     
  S349P <- readRDS('S349P.rds')
  S372I <- readRDS('S372I.rds')     
  S372P <- readRDS('S372P.rds')
  S3I <- readRDS('S3I.rds')     
  S3P <- readRDS('S3P.rds')
  S443I <- readRDS('S443I.rds')     
  S443P <- readRDS('S443P.rds')
  S540I <- readRDS('S540I.rds')     
  S540P <- readRDS('S540P.rds')
  S733I <- readRDS('S733I.rds')     
  S733P <- readRDS('S733P.rds')
  S87I <- readRDS('S87I.rds')     
  S87P <- readRDS('S87P.rds')
  load('sampleall.Rdata')
}
#####直接merge起来，看变化
##一、先merge化疗前的腹膜组织，提取成纤维细胞
##两个稀疏表达矩阵整合
##取出稀疏表达矩阵
I1 <- as.data.frame(S1005I[["RNA"]]@counts)
I2 <- as.data.frame(S136I[["RNA"]]@counts)
I3 <- as.data.frame(S153I[["RNA"]]@counts)
I4 <- as.data.frame(S227I[["RNA"]]@counts)
I5 <- as.data.frame(S349I[["RNA"]]@counts)
I6 <- as.data.frame(S372I[["RNA"]]@counts)
I7 <- as.data.frame(S3I[["RNA"]]@counts)
I8 <- as.data.frame(S443I[["RNA"]]@counts)
I9 <- as.data.frame(S540I[["RNA"]]@counts)
I10 <- as.data.frame(S733I[["RNA"]]@counts)
I11 <- as.data.frame(S87I[['RNA']]@counts)
{I1gene <- rownames(I1)
  I2gene <- rownames(I2)
  I3gene <- rownames(I3)
  I4gene <- rownames(I4)
  row <- intersect(I1gene,I2gene) %>% intersect(I3gene) %>% intersect(I4gene)
  I1 <- I1[row,]
  I2 <- I2[row,]
  I3 <- I3[row,]
  I4 <- I4[row,]
  }#因为基因数一样，所以不运行这个代码
memory.limit(1000000)
OVI <- CreateSeuratObject(counts = cbind(I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11), project = "OVI", min.cells = 5)
saveRDS(OVI,'OVI.rds')
