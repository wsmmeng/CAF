rm(list=ls())
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

stromal <- readRDS(stromal,'stromal.rds')

plot2 <- ElbowPlot(stromal, ndims=50, reduction="pca") 
##降维
##感觉用harmony不好,还是用pca
##本数据集用harmony
library(future)
options(future.globals.maxSize=4000000000)
plan('multiprocess',workers=16)
plan()
stromal <- FindNeighbors(stromal,reduction = "harmony",dims = 1:50) 
stromal <- FindClusters(stromal,reduction = "harmony",resolution = 0.3)
table(stromal@meta.data$seurat_clusters)
metadata <- stromal@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

stromal <- RunUMAP(stromal,reduction = "harmony", dims = 1:50)
stromal <- RunTSNE(stromal,reduction = "harmony", dims = 1:50)

plot2 = DimPlot(stromal, reduction = "umap")
plot1 = DimPlot(stromal, reduction = "tsne")
plot1+plot2
saveRDS(stromal,'stromal.rds')
stromal <- readRDS('stromal.rds')
##注释
plan()
diff.wilcox = FindAllMarkers(stromal)##默认使用wilcox方法挑选差异基因，大概4-5min
all.markers = diff.wilcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
#An adjusted P value < 0.05and | log 2 [fold change (FC)] | > 0.5 
#were considered the 2 cutoff criteria for identifying marker genes.
save(diff.wilcox,file = "allmarkerstromalcluster.Rdata")
saveRDS(stromal,'stromal.rds')
##top10基因绘制热图
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(stromal)) 
DoHeatmap(stromal, features = top10, group.by = "seurat_clusters")
saveRDS(stromal, file="stromal.rds")
