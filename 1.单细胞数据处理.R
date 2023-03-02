rm(list=ls())
setwd('~/合并/spatial/123467/A_955_OvarianTumor/')
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
library(ggplot2)
library(patchwork)

##读取数据，参考https://cloud.tencent.com/developer/article/1892336
matrix_data <- read.csv("raw_counts.csv",header=T, row.names=1)
dim(matrix_data)
colnames(matrix_data)<- gsub('\\.','-',colnames(matrix_data))##这个地方保证和后面读的image是统一的
OV955 <- CreateSeuratObject(counts = matrix_data,
                            assay = "Spatial", # seurat对象中存放表达矩阵的assay名称
                            slice = "slice1"   # seurat对象中设置的空转样本ID，多样本合并分析时有用
)##先创建一个seurat对象
##再把图读进来
img <- Read10X_Image(image.dir = "~/合并/spatial/123467/A_955_OvarianTumor/")
DefaultAssay(object = img) <- 'Spatial'

img <- img[colnames(x = OV955)]
img##看一下img的sample数，要和seurat对象统一
OV955[['image']] <- img



##数据处理
plot1 <- VlnPlot(OV955, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(OV955, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
OV955 <- SCTransform(OV955, assay = "Spatial", verbose = FALSE)
#SpatialFeaturePlot(OV955, features = c("Hpca", "Ttr"))
OV955 <- RunPCA(OV955, assay = "SCT", verbose = FALSE)
OV955 <- FindNeighbors(OV955, reduction = "pca", dims = 1:30)
OV955 <- FindClusters(OV955, verbose = FALSE)
OV955 <- RunUMAP(OV955, reduction = "pca", dims = 1:30)

##可视化特征基因
genes <- c('POSTN',"SFRP2","PRRX1","CTSK",'SULF1')
p1 <- SpatialFeaturePlot(OV955, features = genes, pt.size.factor = 1)
p2 <- SpatialFeaturePlot(OV955, features = "Ttr", alpha = c(0.1, 1))
p1 + p2
##可视化分群图
p1 <- DimPlot(OV955, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(OV955, label = TRUE, label.size = 3)
p1 + p2
##分面可视化



##提取细胞群可视化                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
X <- subset(OV955, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 |
# image_imagecol < 150))
X <- subset(OV955, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
X <- subset(OV955, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
X <- subset(OV955, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(OV955, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(OV955, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2

##常规注释方法--直接看后面的
de_markers <- FindMarkers(OV955, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = OV955, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
##空间转录组的关键基因
OV955 <- FindSpatiallyVariableFeatures(OV955, assay = "SCT", features = VariableFeatures(OV955)[1:1000], 
                                       selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(OV955, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(OV955, features = top.features, ncol = 3, alpha = c(0.1, 1))

##联合单细胞注释
OV <- readRDS("~/合并/OV.rds")
# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k cells
# this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
OV <- readRDS('OV.rds')
OV <- SCTransform(OV, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
OV955 <- SCTransform(OV955, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(OV955, group.by =  "subclass", label = TRUE)
library(future)
plan()
plan('multiprocess',workers=16)

Idents(OV) <- 'celltype2'
anchors <- FindTransferAnchors(reference = OV, query = OV955, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = OV$subclass, prediction.assay = TRUE, 
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
anchors <- FindTransferAnchors(reference = OV, query = OV955, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = OV$celltype, prediction.assay = TRUE, 
                                  weight.reduction = OV955[["pca"]], dims = 1:30)
OV955[["predictions"]] <- predictions.assay
DefaultAssay(OV955) <- "predictions"
SpatialFeaturePlot(OV955, features = c("", ""), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
##所有细胞亚群
SpatialFeaturePlot(OV955, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", 
                                    "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
##不明白什么意思,和上面应该是一样的
OV955 <- FindSpatiallyVariableFeatures(OV955, assay = "predictions", selection.method = "markvariogram", 
                                        features = rownames(OV955), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(OV955), 4)
SpatialPlot(object = OV955 , features = top.clusters, ncol = 2)

##感觉这个注释不是很科学，用addmodulescorecelltype注释一下
load('~/合并/aftercombinedcelltypemarkers.Rdata')
epimarker <- rownames(subset(markers,cluster=='epithelial tumor cell'))
resistance_features <- list(epimarker)
resistancescore <- AddModuleScore(OV955,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")
colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[8] <- 'resistance_Score'
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(1, 1))

DCmarker <- rownames(subset(markers,cluster=='demosplastic CAFs'))
table(markers$cluster)
resistance_features <- list(DCmarker)
resistancescore <- AddModuleScore(OV955,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")
colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[8] <- 'resistance_Score'
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(1, 1))

ICmarker <- rownames(subset(markers,cluster=='inflammatory CAFs'))
table(markers$cluster)
resistance_features <- list(DCmarker)
resistancescore <- AddModuleScore(OV955,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")
colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[8] <- 'resistance_Score'
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(1, 1))
saveRDS(OV955,'OV955.rds')
