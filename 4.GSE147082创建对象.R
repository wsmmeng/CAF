##合并对象--数据整合成一个seurat对象，看看需不需要harmony
rm(list=ls())
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")

#读取数据
matrix_data1 <- read.csv("PT-2834.csv", header=T, row.names=1)
dim(matrix_data1)
seurat_obj1 <- CreateSeuratObject(counts = matrix_data1)

matrix_data2 <- read.csv("PT-3232.csv", header=T, row.names=1)
dim(matrix_data2)
seurat_obj2 <- CreateSeuratObject(counts = matrix_data2)

matrix_data3 <- read.csv("PT-3401.csv", header=T, row.names=1)
dim(matrix_data3)
seurat_obj3 <- CreateSeuratObject(counts = matrix_data3)

matrix_data4 <- read.csv("PT-4806.csv", header=T, row.names=1)
dim(matrix_data4)
seurat_obj4 <- CreateSeuratObject(counts = matrix_data4)

matrix_data5 <- read.csv("PT-5150.csv", header=T, row.names=1)
dim(matrix_data5)
seurat_obj5 <- CreateSeuratObject(counts = matrix_data5)

matrix_data6 <- read.csv("PT-6885.csv", header=T, row.names=1)
dim(matrix_data6)
seurat_obj6 <- CreateSeuratObject(counts = matrix_data6)

##多样本数据整合
###########harmony 速度快、内存少################
##参考https://blog.csdn.net/qazplm12_3/article/details/104791765
##参考《“harmony”整合不同平台的单细胞数据之旅》

##取出稀疏表达矩阵
P1 <- as.data.frame(seurat_obj1[["RNA"]]@counts)
P2 <- as.data.frame(seurat_obj2[["RNA"]]@counts)
P3 <- as.data.frame(seurat_obj3[["RNA"]]@counts)
P4 <- as.data.frame(seurat_obj4[["RNA"]]@counts)
P5 <- as.data.frame(seurat_obj5[["RNA"]]@counts)
P6 <- as.data.frame(seurat_obj6[["RNA"]]@counts)

p1gene <- rownames(P1)
P2gene <- rownames(P2)
p3gene <- rownames(P3)
p4gene <- rownames(P4)
p5gene <- rownames(P5)
p6gene <- rownames(P6)
row <- intersect(p1gene,P2gene) %>% intersect(p3gene) %>% intersect(p4gene) %>% intersect(p5gene) %>% intersect(p6gene)
P1 <- P1[row,]
P2 <- P2[row,]
P3 <- P3[row,]
P4 <- P4[row,]
P5 <- P5[row,]
P6 <- P6[row,]
colnames(P1) <- paste(colnames(P1),'p1',sep = '-')
colnames(P2) <- paste(colnames(P2),'p2',sep = '-')
colnames(P3) <- paste(colnames(P3),'p3',sep = '-')
colnames(P4) <- paste(colnames(P4),'p4',sep = '-')
colnames(P5) <- paste(colnames(P5),'p5',sep = '-')
colnames(P6) <- paste(colnames(P6),'p6',sep = '-')
##创建seurat对象
memory.limit(1000000)
##这个地方要改一下，因为线粒体这里是.,不是-
P <- cbind(P1,P2,P3,P4,P5,P6)
rownames(P) <- gsub("\\.","-",rownames(P))
##
OVP <- CreateSeuratObject(counts = P, project = "OVP", min.cells = 5)

saveRDS(OVP,'OVP.rds')
