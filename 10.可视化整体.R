rm(list=ls())
library(monocle)
library(Seurat)
OV <- readRDS('OV.rds')

##Dimplot可视化
##先可视化cluster##已经画了不需要重新画

library(paletteer) #提供了 R 编程环境中提供的数百种其他调色板的组合集合，详情可以查看此包，总有你满意的方案  
Idents(OV) <- 'celltype'
table(OV$celltype)
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)] #有几群细胞需要标记就选几种颜色  
DimPlot(OV, label = T, reduction = 'umap',
        #cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)#标注有点挤，repel=T可以让排列更加合理
library(paletteer) #提供了 R 编程环境中提供的数百种其他调色板的组合集合，详情可以查看此包，总有你满意的方案  
Idents(OV) <- 'celltype'
table(OV$celltype)
#pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5)] #有几群细胞需要标记就选几种颜色  
DimPlot(OV, label = T,   
        #cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)#标注有点挤，repel=T可以让排列更加合理  

Idents(OV) <- 'sample'
#pal <- paletteer_d("ggsci::nrc_npg")[c(4:39)] #有几群细胞需要标记就选几种颜色  
DimPlot(OV, label = F,   
        #cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)#+#标注有点挤，repel=T可以让排列更加合理  
  #theme(legend.position="top")

##画个热图##热图和点图都画了，不需要改
rm(list=ls())
##参考《单细胞亚群 Marker 基因热图重绘及均值展示》
#取细胞子集
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

OV<- readRDS('OV.rds')
Idents(OV) <- 'celltype'
markers <- FindAllMarkers(OV,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
Idents(OV) <- 'seurat_clusters'
markers <- FindAllMarkers(OV,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)

# get top 10 genes
library(dplyr)
top10markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
# get color
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
col <- pal_npg()(5)
DoHeatmap(OV,features = top10markers$gene,label = F,
          group.colors = col) +
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')


###########################################

# get cells mean gene expression
mean_gene_exp <- AverageExpression(OV,
                                   features = top10markers$gene,
                                   group.by = 'celltype',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# add colnames
colnames(mean_gene_exp) <- gsub('RNA.','',colnames(mean_gene_exp))

# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))

# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))

# top annotation
anno_col <- pal_npg()(3)
anno_col <- brewer.pal(3, "Paired")
names(anno_col) <- colnames(mean_gene_exp)
column_ha = HeatmapAnnotation(celltype= colnames(htdf),
                              col = list(celltype = anno_col))
# plot
Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_title = "celltypemarkers",
        column_title = "Celltype",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        show_column_names = F,
        column_names_side = 'top',
        column_names_rot = 0,
        top_annotation = column_ha,
        # column_split = paste('clsuter ',0:8,sep = ''),
        col = col_fun)

##画个气泡图
markers <- c("FAP","ACTA2",'DCN','COL1A2',
             "CD8A","CD3D",
             "FCGR3B",'CXCL8',
             "MUC16","EPCAM",
             'MS4A1','CD79A'
)
DotPlot(OV,features = markers)+coord_flip()+
  scale_color_gradientn(values = seq(0,1,0.2),colours=c('#330066','#336699','red'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

##画个marker图
#修饰图片
##上皮细胞
markers1 <- c("WFDC2","PAX8","EPCAM")
##基质细胞
markers2 <- c("COL1A2","FGFR1","DCN")
##免疫细胞
markers3 <- c("FCER1G", "PTPRC")
marker <- c(markers1,markers2,markers3)
mycolor <- c('lightblue', 'white','red')#设置颜色 
FeaturePlot(OV, features = marker,cols = mycolor, pt.size = 1.5,max.cutoff = 2.5)#+  
#theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框
FeaturePlot(OV, features = c("ACTA2"),cols = mycolor, pt.size = 1.5,max.cutoff = 1)  

FeaturePlot(OV, features = c("FAP"),cols = mycolor, pt.size = 1.5,max.cutoff = 1)
FeaturePlot(OV, features = c("PDGFRA"),cols = mycolor, pt.size = 1.5,max.cutoff = 1)
FeaturePlot(OV, features = c("PDGFRB"),cols = mycolor, pt.size = 1.5,max.cutoff = 1.8)
FeaturePlot(OV, features = c("PDPN"),cols = mycolor, pt.size = 1.5,max.cutoff = 1.8)
FeaturePlot(OV, features = c("THY1"),cols = mycolor, pt.size = 1.5,max.cutoff = 2.5)
FeaturePlot(OV, features = c("COL1A1"),cols = mycolor, pt.size = 1.5,max.cutoff = 2.5)

DotPlot(OV,features = marker)+coord_flip()+
  scale_color_gradientn(values = seq(0,1,0.2),colours=c('#330066','#336699','red'))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

Idents(hpvp) <- 'seurat_clusters'
DotPlot(hpvp,features = markers)+coord_flip()