##注释
rm(list=ls())
rm(list=ls())
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
library(ggplot2)
library(patchwork)

OV <- readRDS('OV.rds')
load("allmarkerovariancluster.Rdata")
##不要用singleR注释空转，感觉效果不好
markers <- c("FAP","ACTA2","CD8A","CD3D","FCGR3B",'CXCL8',
             "CD44","EPCAM","PECAM1")
DotPlot(OV,features = markers)+coord_flip()
FeaturePlot(OV, features = markers)



#GSE165897的情况epithelial cancer cells (WFDC2, PAX8, and EPCAM), stromal cells (COL1A2, FGFR1, and DCN), and immune cells (CD79A, FCER1G, and PTPRC)

##上皮细胞
markers <- c("WFDC2","PAX8","EPCAM")
FeaturePlot(OV, features = markers,reduction = 'umap')
DimPlot(OV, reduction = 'umap',label = T,label.size = 2)
VlnPlot(OV, features = markers)
##基质细胞
markers <- c("COL1A2","FGFR1","DCN")
FeaturePlot(OV, features = markers)
VlnPlot(OV, features = markers)
##免疫细胞
markers <- c("CD79A", "FCER1G", "PTPRC")
FeaturePlot(OV, features = markers)
VlnPlot(OV, features = markers)
##成纤维细胞
FeaturePlot(OV, features = c('FAP','ACTA2'))
VlnPlot(OV, features = c('FAP'))#4,7
##
FeaturePlot(OV, features = c('ACTA2'))
VlnPlot(OV, features = c('ACTA2'))#4,5,7

##
FeaturePlot(OV, features = c('PECAM1'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('PECAM1'))
FeaturePlot(OV, features = c('DCN'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('DCN'))
FeaturePlot(OV, features = c('COL1A2'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('COL1A2'))

##上皮细胞
FeaturePlot(OV, features = c('MUC16'))
VlnPlot(OV, features = c('MUC16'))##3,5,6,7,9
FeaturePlot(OV, features = c('KRT19'))
VlnPlot(OV, features = c('KRT19'))##
FeaturePlot(OV, features = c('EPCAM'))
VlnPlot(OV, features = c('EPCAM'))
##T细胞
FeaturePlot(OV, features = c('CD3D'))
VlnPlot(OV, features = c('CD3D'))#0,8
FeaturePlot(OV, features = c('CD8A'))
VlnPlot(OV, features = c('CD8A'))#0,CT8
FeaturePlot(OV, features = c('CD4'))
VlnPlot(OV, features = c('CD4'))## 3

##Neu
FeaturePlot(OV, features = c('MNDA'))
VlnPlot(OV, features = c('MNDA'))
FeaturePlot(OV, features = c('FCGR3B'))
VlnPlot(OV, features = c('FCGR3B')) ## 10
FeaturePlot(OV, features = c('CXCL8'))
VlnPlot(OV, features = c('CXCL8'))
##B细胞，没分出来
FeaturePlot(OV, features = c('MS4A1'))
VlnPlot(OV, features = c('MS4A1'))#4
FeaturePlot(OV, features = c('CD40'))
VlnPlot(OV, features = c('CD40'))
FeaturePlot(OV, features = c('IGHM'))
VlnPlot(OV, features = c('IGHM'))
FeaturePlot(OV, features = c('CD79A'))
VlnPlot(OV, features = c('CD79A'))#4,10

##基质细胞
FeaturePlot(OV, features = c('PECAM1'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('PECAM1'))
FeaturePlot(OV, features = c('DCN'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('DCN'))
FeaturePlot(OV, features = c('COL1A2'), alpha = c(0.1, 1))
VlnPlot(OV, features = c('COL1A2'))
##综上 
#scedata <- subset(scedata, idents = c("21"), invert = TRUE)#去掉低质量细胞群
new.cluster.ids <- c("0"="T cell", 
                     "1"="Neutrophil", 
                     "2"="Fibroblast", 
                     "3"="epithelial", 
                     "4"="B cell", 
                     "5"="epithelial", 
                     "6"="epithelial", 
                     "7"="epithelial", 
                     "8"="T cell", 
                     "9"="epithelial", 
                     "10"="B cell"
)
new.cluster.ids <- c("0"="immune cell", 
                     "1"="immune cell", 
                     "2"="stromal cell", 
                     "3"="epithelial tumor cell", 
                     "4"="epithelial tumor cell", 
                     "5"="immune cell", 
                     "6"="epithelial tumor cell", 
                     "7"="epithelial tumor cell", 
                     "8"="immune cell", 
                     "9"="stromal cell", 
                     "10"="immune cell",
                     '11'= "stromal cell",
                       '12'="immune cell",
                     '13'="epithelial tumor cell",
                     '14'="immune cell",
                     '15'="epithelial tumor cell",
                     '16'="epithelial tumor cell",
                     '17'="epithelial tumor cell",
                     '18'="epithelial tumor cell",
                     '19'="epithelial tumor cell",
                     '20'="epithelial tumor cell",
                     '21'="epithelial tumor cell",
                     '22'="epithelial tumor cell",
                     '23'="stromal cell",
                     '24'="stromal cell",
                     '25'="epithelial tumor cell",
                     '26'="epithelial tumor cell"
                     
)
markers <- c("FAP","ACTA2",
             "CD8A","CD3D",
             "FCGR3B",'CXCL8',
             "MUC16","EPCAM",
             'MS4A1','CD79A',
             'DCN','COL1A2')
##把成纤维细胞的marker单独提取出来
FeaturePlot(OV, features = c('FAP'))
FeaturePlot(OV, features = c('ACTA2'))
##命个名
Idents(OV) <- 'seurat_clusters'
OV <- RenameIdents(OV , new.cluster.ids)                        
OV$celltype <- OV@active.ident
##画个图
DimPlot(OV, group.by = "celltype",reduction = 'umap',label = T)
DimPlot(OV, group.by = "celltype",reduction = 'tsne',label = T)
DimPlot(OV, facet.highlight = TRUE, ncol = 3,
        cells.highlight = CellsByIdentities(OV))
DimPlot(OV, label = TRUE, label.size = 3,alpha = c(0.6, 1)) 

markers <- c('PRRX1','POSTN', 'SFRP2', 'CTSK', 'SULF1')
FeaturePlot(OV, features = markers)
FeaturePlot(OV, features = c('PRRX1'))
FeaturePlot(OV, features = c('POSTN'))
FeaturePlot(OV, features = c('SFRP2'))
FeaturePlot(OV, features = c('CTSK'))
FeaturePlot(OV, features = c('SULF1'))
saveRDS(OV,file='OV.rds')
