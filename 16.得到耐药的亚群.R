rm(list=ls())
##参考《单细胞亚群 Marker 基因热图重绘及均值展示》
##参考《如何改造你的图片，让你的单细胞测序分析图向CNS看齐？》
#取细胞子集
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

stromal<- readRDS('stromal.rds')

##GSE63885分布
DefaultAssay(stromal) <- "RNA"
load('degdrug.Rdata')
##因为对比矩阵是Resistent对sensitive,所以上调的是铂耐药
UP <- subset(degdrug,g=='UP')
gene_symbol <- rownames(UP)
gene_symbol[c(4,31)] <- c("WASIR1",'NTM')
resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(stromal,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[20] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('LightSteelBlue',"CornflowerBlue",'MediumOrchid','DeepPink'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

##降维图
p2 <- DimPlot(stromal, reduction = 'umap', group.by = "seurat_clusters", label=T)
p1+p2

##GSE30161分布
DefaultAssay(stromal) <- "RNA"
load('degdrug30161.Rdata')
##因为对比矩阵是CR对PDPR,所以下调的是铂耐药
DOWN <- subset(degdrug30161,g=='DOWN')
gene_symbol <- rownames(DOWN)
gene_symbol[c(7,9,30,33,34,36)] <- c("LINC01000",'H19','IGF2','TAX1BP3','DUXAP10','NNMT')

resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(stromal,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[20] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


DefaultAssay(stromal) <- "RNA"
load('degdrug23554.Rdata')
##因为对比矩阵是CR对IR,所以下调的是铂耐药
DOWN <- subset(degdrug23554,g=='DOWN')
gene_symbol <- rownames(DOWN)
gene_symbol[2] <- c("NNMT")

resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(stromal,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[20] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "celltype")

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

##降维图
p1 <- DimPlot(OVfsub, reduction = 'tsne', group.by = "seurat_clusters", label=T)
p2 <- DimPlot(stromal, reduction = 'umap', group.by = "seurat_clusters", label=T)
p1+p2
FeaturePlot(stromal,)