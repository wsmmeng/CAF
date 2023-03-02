rm(list=ls())
OV <- readRDS('OV.rds')
#取细胞子集
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)
##提取细胞子集
table(OV@meta.data$celltype)
Cells.sub <- subset(OV@meta.data, celltype=="immune cell")
immune <- subset(OV, cells=row.names(Cells.sub))
saveRDS(immune,'immune.rds')


##进行降低维度聚类
immune <- readRDS(immune,'immune.rds')
plot2 <- ElbowPlot(immune, ndims=50, reduction="pca") 
##降维
##感觉用harmony不好,还是用pca
##本数据集用harmony
library(future)
options(future.globals.maxSize=4000000000)
plan('multiprocess',workers=16)
plan()
immune <- FindNeighbors(immune,reduction = "harmony",dims = 1:50) 
immune <- FindClusters(immune,reduction = "harmony",resolution = 0.3)
table(immune@meta.data$seurat_clusters)
metadata <- immune@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

immune <- RunUMAP(immune,reduction = "harmony", dims = 1:50)
immune <- RunTSNE(immune,reduction = "harmony", dims = 1:50)

plot2 = DimPlot(immune, reduction = "umap")
plot1 = DimPlot(immune, reduction = "tsne")
plot1+plot2
saveRDS(immune,'immune.rds')
immune <- readRDS('immune.rds')
##注释
plan()
diff.wilcox = FindAllMarkers(immune)##默认使用wilcox方法挑选差异基因，大概4-5min
all.markers = diff.wilcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
#An adjusted P value < 0.05and | log 2 [fold change (FC)] | > 0.5 
#were considered the 2 cutoff criteria for identifying marker genes.
save(diff.wilcox,file = "allmarkerimmunecluster.Rdata")
saveRDS(immune,'immune.rds')
##top10基因绘制热图
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10
top10 = CaseMatch(search = as.vector(top10$gene), match = rownames(immune)) 
DoHeatmap(immune, features = top10, group.by = "seurat_clusters")
saveRDS(immune, file="immune.rds")


##进行注释
top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
top30gene <- top20$gene
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top10gene <- top10$gene
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50gene <- top50$gene


NK <- c('NCAM1','KLRD1','NKG7',"GNLY")
FeaturePlot(immune,features = NK)
DimPlot(immune, reduction = "umap")##cluster5
Neu <- c('FCGR3B','CXCL8','MNDA','SELL','CXCR2','FCGR3A')
FeaturePlot(immune,features = Neu)
DimPlot(immune, reduction = "umap")##cluster3,7
Mac <- c('CD14','CD163','CD68','CSF1R')
FeaturePlot(immune,features = Mac)
DimPlot(immune, reduction = "umap")##cluster2

CD8 <- c('CD8A','CD8B','CD3D','IL7R') 
FeaturePlot(immune,features = CD8)
DimPlot(immune, reduction = "umap")##cluster0,1,6
Effectormemory <- c("IL7R","EOMES")
FeaturePlot(immune,features = Effectormemory )
DimPlot(immune, reduction = "umap")##cluster1
Exhauted <- c("LAG3","TIGIT")
FeaturePlot(immune,features = Exhauted  )
DimPlot(immune, reduction = "umap")##cluster6
Effector <- c("PRF1","CCL5")
FeaturePlot(immune,features = Effector )
DimPlot(immune, reduction = "umap")##cluster1


Treg <- c("FOXP3","CTLA4","IL2RA")
FeaturePlot(immune,features = Treg)
DimPlot(immune, reduction = "umap")##cluster4

B <- c("CD19","CD79A",'MS4A1','CD40','IGHM','HLA-DRA','CCR7')
FeaturePlot(immune,features = B)
DimPlot(immune, reduction = "umap")
B2 <- c("IGKC","MZB1")
FeaturePlot(immune,features = B2)
DimPlot(immune, reduction = "umap")##cluster7


##unknown cluster8

celltypemarker <- c('NKG7',"GNLY",##先鉴定天然免疫细胞NK ##cluster5
                    #"CXCL8",'MNDA',
                    'CD14','CD163','CSF1R',##leukocyte##cluster2,3,7##只能分为白细胞了
                    'CD8A','CD8B','IL7R',"GZMB",##再鉴定适应性免疫细胞#cluster0,1,6
                    "FOXP3","CTLA4","IL2RA",##Tregcluster4
)


new.cluster.ids <- c("0"="CD8+ T cells", 
                     "1"="CD8+ T cells", 
                     "2"="leukocytes", 
                     "3"="leukocytes", 
                     "4"="Treg",
                     "5"='NK cells',
                     "6"="CD8+ T cells", 
                     "7"="leukocytes", 
                     "8"="low quality cells")

##命个名
Idents(immune) <- 'seurat_clusters'
immune <- RenameIdents(immune , new.cluster.ids)                        
immune$celltype<- immune@active.ident
##画个图
DimPlot(immune, group.by = "celltype",reduction = 'umap',label = T)

saveRDS(immune,file='immune.rds')



