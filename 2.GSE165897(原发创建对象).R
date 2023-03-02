rm(list = ls())
{library(multtest)
  if(!require(multtest))install.packages("multtest")
  if(!require(Seurat))install.packages("Seurat")
  if(!require(dplyr))install.packages("dplyr")
  if(!require(mindr))install.packages("mindr")
  if(!require(mindr))install.packages("tidyverse")}
##读取数据
{S1005P <- readRDS('S1005P.rds')
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
P1 <- as.data.frame(S1005P[["RNA"]]@counts)
P2 <- as.data.frame(S136P[["RNA"]]@counts)
P3 <- as.data.frame(S153P[["RNA"]]@counts)
P4 <- as.data.frame(S227P[["RNA"]]@counts)
P5 <- as.data.frame(S349P[["RNA"]]@counts)
P6 <- as.data.frame(S372P[["RNA"]]@counts)
P7 <- as.data.frame(S3P[["RNA"]]@counts)
P8 <- as.data.frame(S443P[["RNA"]]@counts)
P9 <- as.data.frame(S540P[["RNA"]]@counts)
P10 <- as.data.frame(S733P[["RNA"]]@counts)
P11 <- as.data.frame(S87P[['RNA']]@counts)
{P1gene <- rownames(P1)
P2gene <- rownames(P2)
p3gene <- rownames(P3)
p4gene <- rownames(P4)
row <- intersect(p1gene,P2gene) %>% intersect(p3gene) %>% intersect(p4gene)
P1 <- P1[row,]
P2 <- P2[row,]
P3 <- P3[row,]
P4 <- P4[row,]
}#因为基因数一样，所以不运行这个代码
memory.limit(100000)
OVPPP <- CreateSeuratObject(counts = cbind(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11), project = "OVPPP", min.cells = 5)

saveRDS(OVPPP,'OVPPP.rds')
