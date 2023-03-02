rm(list = ls (all = TRUE))

library(monocle)
library(Seurat)
stromal <- readRDS('stromal.rds')
load("allmarkerstromalcluster.Rdata")
library(dplyr)
{diff.wilcox = FindAllMarkers(stromal)##默认使用wilcox方法挑选差异基因，大概4-5min
  all.markers = diff.wilcox %>% select(gene, everything()) %>%
    subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
}

top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
top30gene <- top20$gene
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10gene <- top10$gene
top50 = all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50gene <- top50$gene
##
{myCAF <- c('ACTA2','MYH11','MCAM','TAGLN','MYLK') #myofibroblasts
  deCAF <- c('COL1A1','COL3A1','DCN','PDPN') #desmoplastic fibroblasts
  apCAF <- c('HLA-DPA1','HLA-DPA','HLA-DQA1','SLPI') #antigen-presenting fibroblasts
  iCAF <- c('IL6,CFD','CXCL2','LMNA') #inflammatory fibroblasts
  normal <- 
  low-quality
  pCAF <- c('BIRC5','TOP2A')
  }
VlnPlot(stromal,features = myCAF)#3,5
FeaturePlot(stromal,features = myCAF)
VlnPlot(stromal,features = c('COL1A1','COL3A1','DCN'))#0
VlnPlot(stromal,features = c("CFD", "C3", "CXCL14","CXCL12"))#1
VlnPlot(stromal,features = c("TOP2A",'BIRC5'))#4
VlnPlot(stromal,features = c('FAP'))

VlnPlot(stromal,features = c('ACTA2','MYH11','MCAM','MYLK'))#3,5
VlnPlot(stromal,features = c('COL1A1','COL3A1','DCN'))#0
VlnPlot(stromal,features = c("CFD", "C3", "CXCL14","CXCL12"))#1
VlnPlot(stromal,features = c("TOP2A",'BIRC5'))#4
VlnPlot(stromal,features = c('FAP'))


celltypemarker <- c('CTHRC1',"FN1","COL11A2","COL1A1",##dCAF
                    "RSAD1",'PDPN',##lymphatic endothelial cell
                    'CFD','APOE','C7',"IL6",##inflammatory CAF
                    "LOX",##mesenchymal stem cell
                    "RGS5","MCAM","MYH11","ACTA2",##myofibroblast CAF
                    'ANGPT2','PECAM1','VWF')##vascular endothelial cell
##命名 ##myofibroblasts
new.cluster.ids <- c("0"="dCAF1", 
                     "1"="lymphatic endothelial cell", 
                     "2"="iCAF1", 
                     "3"="iCAF2", 
                     "4"="mesenchymal stem cell",
                     "5"='iCAF3',
                     "6"="myo CAF", 
                     "7"="vascular endothelial cells", 
                     "8"="dCAF2", 
                     "9"="low quality cells", 
                     "10"="low quality cells"
                    
)
new.cluster.ids <- c("0"="demosplastic CAF", 
                     "1"="lymphatic endothelial cell", 
                     "2"="inflammatory CAF", 
                     "3"="inflammatory CAF", 
                     "4"="mesenchymal stem cell",
                     "5"='inflammatory CAF',
                     "6"="myofibroblast CAF", 
                     "7"="vascular endothelial cell", 
                     "8"="demosplastic CAF", 
                     "9"="demosplastic CAF", 
                     "10"="low quality cell"
                     
)
##把成纤维细胞的marker单独提取出来
FeaturePlot(stromal, features = c('FAP'))
FeaturePlot(stromal, features = c('ACTA2'))
FeaturePlot(stromal, features = c('DCN'))
FeaturePlot(stromal, features = c('COL1A1'))
FeaturePlot(stromal, features = c('COL1A2'))
FeaturePlot(stromal,features = c('POSTN','SFRP2','PRRX1','SULF1','CTSK'))
##命个名
Idents(stromal) <- 'seurat_clusters'
stromal <- RenameIdents(stromal , new.cluster.ids)                        
stromal$celltype2 <- stromal@active.ident
stromal$celltype <- stromal@active.ident
##画个图
DimPlot(stromal, group.by = "celltype",reduction = 'umap',label = T)
p = DoHeatmap(stromal, features = top20gene, group.by = "seurat_clusters")
p
p = DoHeatmap(stromal, features = top10gene, group.by = "seurat_clusters")
p
p = DoHeatmap(stromal, features = top10gene, group.by = "celltype")
p
saveRDS(stromal,'stromal.rds')

