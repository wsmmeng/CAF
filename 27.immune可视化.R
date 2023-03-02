rm(list=ls())
##参考《单细胞亚群 Marker 基因热图重绘及均值展示》
##参考《如何改造你的图片，让你的单细胞测序分析图向CNS看齐？》
#取细胞子集
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
library(clusterProfiler)

immune<- readRDS('immune.rds')
Idents(immune) <- 'celltype'
celltypemarkersimmune <- FindAllMarkers(immune,
                          only.pos = TRUE,
                          min.pct = 0,
                          logfc.threshold = 0)
save(markers,file='celltypemarkersimmune.Rdata')
{Idents(immune) <- 'seurat_clusters'
  markers <- FindAllMarkers(immune,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
}##前面已经运行过

save(celltypemarkersimmune,file='celltypemarkersimmune.Rdata')

# get top 10 genes
ibrary(future)
options(future.globals.maxSize=4000000000)
plan('multiprocess',workers=16)
top10markers <- celltypemarkersimmune %>%
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
table(immune@meta.data$celltype)
col <- pal_npg()(7)
DoHeatmap(immune,features = top10markers$gene,
          group.colors = col,label = F) +
  scale_colour_npg() +
  scale_fill_gradient2(low = '#0099CC',mid = 'white',high = '#CC0033',
                       name = 'Z-score')




library(paletteer) #提供了 R 编程环境中提供的数百种其他调色板的组合集合，详情可以查看此包，总有你满意的方案  
Idents(immune) <- 'celltype'
table(immune@meta.data$celltype)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5,2)] #有几群细胞需要标记就选几种颜色  
pal <- pal_npg()(5)
DimPlot(immune, label = T,   
        cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)+#标注有点挤，repel=T可以让排列更加合理  
  theme(legend.position="top")
Idents(immune) <- 'sample'
pal <- paletteer_d("ggsci::nrc_npg")[c(4:35)] #有几群细胞需要标记就选几种颜色

DimPlot(immune, label = F,   
        cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)+#标注有点挤，repel=T可以让排列更加合理  
  theme(legend.position="top")
DimPlot(immune, label = F,   
        cols= pal, #设置颜色  
        pt.size = 1.5,#设置点的大小  
        repel = T)+#标注有点挤，repel=T可以让排列更加合理  
  NoLegend() 
#默认作图  
FeaturePlot(immune, features = c("POSTN","SFRP2","PRRX1","CTSK","SULF1"))  
#修饰图片  
mycolor <- c('gray', 'white','red')#设置颜色  
FeaturePlot(immune, features = c("POSTN"),cols = mycolor, pt.size = 1,max.cutoff = 2.5)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框
FeaturePlot(immune, features = c("SFRP2"),cols = mycolor, pt.size = 1,max.cutoff = 2.5)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框
FeaturePlot(immune, features = c("PRRX1"),cols = mycolor, pt.size = 1,max.cutoff = 1.8)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框
FeaturePlot(immune, features = c("CTSK"),cols = mycolor, pt.size = 1,max.cutoff = 1.8)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框
FeaturePlot(immune, features = c("SULF1"),cols = mycolor, pt.size = 1,max.cutoff = 1.8)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框

###画个点图
celltypemarker <- c('NKG7',"GNLY",##先鉴定天然免疫细胞NK ##cluster5
                    "CXCL8",'MNDA',
                    'CD14','CD163','CSF1R',##leukocyte##cluster2,3,7##只能分为白细胞了
                    'CD8A','CD8B','IL7R',"GZMB","NKG7",##再鉴定适应性免疫细胞#cluster0,1,6
                    "FOXP3","CTLA4","IL2RA"##Tregcluster4
)
Idents(immune) <- 'celltype'
DotPlot(immune, features = markers)+coord_flip()+theme_bw()+#去除背景，旋转图片  
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#文字90度呈现  
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#颜色渐变设置  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

###
ViolinPlot <- function(object, groupBy, MarkerSelected) {
  # (1)获取绘图数据1
  plot_data = FetchData(object = object, 
                        vars = c(MarkerSelected$gene, groupBy), 
                        slot = 'data') %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')
  
  # (2)获取绘图数据2
  ident_plot = MarkerSelected %>% 
    dplyr::select(cluster, gene)
  
  # (3)绘图
  figure_1 = ggplot(data = plot_data, mapping = aes(x = Expr,
                                                    y = fct_rev(factor(x = Feat, 
                                                                       levels = MarkerSelected$gene)), 
                                                    fill = group, 
                                                    label = group)) +
    geom_violin(scale = 'width', adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
    facet_grid(cols = vars(group), scales = 'free') +
    cowplot::theme_cowplot(font_family = 'Arial') +
    scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3')) +
    xlab('Expression Level') + 
    ylab('') +
    theme(legend.position = 'none', 
          panel.spacing = unit(x = 0, units = 'lines'),
          axis.line = element_blank(), #去除x和y轴坐标线(不包括axis tick)；
          panel.background = element_rect(fill = NA, color = 'black'),
          strip.background = element_blank(), #去除分页题头背景；
          strip.text = element_text(color = 'black', size = 10, family = 'Arial', face = 'bold'),
          axis.text.x = element_text(color = 'black', family = 'Arial', size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(color = 'black', family = 'Arial', size = 15),
          axis.ticks.x = element_line(color = 'black', lineend = 'round'),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(x = 0.1, units = 'cm'))
  
  figure_2 = ggplot(data = ident_plot, aes(x = 1,
                                           y = fct_rev(factor(x = gene, levels = MarkerSelected$gene)),
                                           fill = cluster)) +
    geom_tile() +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3')) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    guides(fill = guide_legend(direction = 'vertical',
                               label.position = 'right',
                               title.theme = element_blank(),
                               keyheight = 0.5,
                               nrow = 2)) +
    xlab('Feature') +
    theme(legend.text = element_text(family = 'Arial', color = 'black', size = 11),
          legend.position = 'bottom',
          legend.justification = 'left',
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,05,0,0),
          panel.spacing = unit(0, 'lines'),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(x = c(0,0,0,0), units = 'cm'),
          axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, color = 'black', family = 'Arial'),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  figure_2 + figure_1 + patchwork::plot_layout(nrow = 1, widths = c(0.03, 0.97))
}
VlnPlot(immune,features = celltypemarker,pt.size = 0,ncol = 4)

celltypemarker <- c('NKG7',"GNLY",##先鉴定天然免疫细胞NK ##cluster5
                    "CXCL8",'MNDA',
                    'CD14','CD163','CSF1R',##leukocyte##cluster2,3,7##只能分为白细胞了
                    'CD8A','CD8B','IL7R',"GZMB","NKG7",##再鉴定适应性免疫细胞#cluster0,1,6
                    "FOXP3","CTLA4","IL2RA"##Tregcluster4
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
