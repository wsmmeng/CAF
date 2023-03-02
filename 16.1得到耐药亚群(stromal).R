library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

OV<- readRDS('OV.rds')

##GSE63885分布
DefaultAssay(OV) <- "RNA"
load('degdrug.Rdata')
##因为对比矩阵是Resistent对sensitive,所以上调的是铂耐药
UP <- subset(degdrug,g=='UP')
gene_symbol <- rownames(UP)
gene_symbol[c(4,31)] <- c("WASIR1",'NTM')
resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(OV,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[17] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))


a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('LightSteelBlue',"CornflowerBlue",'DeepPink','Red'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

##降维图
p2 <- DimPlot(OV, reduction = 'umap', group.by = "seurat_clusters", label=T)
p1+p2

##GSE30161分布
DefaultAssay(OV) <- "RNA"
load('degdrug30161.Rdata')
##因为对比矩阵是CR对PDPR,所以下调的是铂耐药
DOWN <- subset(degdrug30161,g=='DOWN')
gene_symbol <- rownames(DOWN)
gene_symbol[c(7,9,30,33,34,36)] <- c("LINC01000",'H19','IGF2','TAX1BP3','DUXAP10','NNMT')

resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(OV,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[17] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('LightSteelBlue',"CornflowerBlue",'DeepPink','Red'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


DefaultAssay(OV) <- "RNA"
load('degdrug23554.Rdata')
##因为对比矩阵是CR对IR,所以下调的是铂耐药
DOWN <- subset(degdrug23554,g=='DOWN')
gene_symbol <- rownames(DOWN)
gene_symbol[2] <- c("NNMT")

resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(OV,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[17] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('LightSteelBlue',"CornflowerBlue",'DeepPink','Red'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

##降维图
p1 <- DimPlot(OVfsub, reduction = 'tsne', group.by = "seurat_clusters", label=T)
p2 <- DimPlot(OV, reduction = 'umap', group.by = "seurat_clusters", label=T)
p1+p2
FeaturePlot(OV,)