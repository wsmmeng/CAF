rm(list=ls())
##gene
hypoxia <- read.csv('hypoxia.csv')
angio <- read.csv('ANGIOGENESIS.csv')

#AddModuleScore，参考《新贵分析！单细胞联合空转分析，R语言手把手教学，你学废了吗？》
#选择基因列表，这里因为我们并没有感兴趣的基因，所以就常规选择了
#如果有感兴趣的基因列表（比如免疫、转移、炎症等等）可以提取感兴趣的基因列表然后映射到空间转录组上
#进行AddModuleScore打分只需要两个文件：Genelist以及Visium_Seurat对象
#转换成列表
rm(list=ls())
OV955 <- readRDS('OV955.rds')
load('~/合并/degdrug.Rdata')
##GSE63885分布
##因为对比矩阵是Resistent对sensitive,所以上调的是铂耐药
UP <- subset(degdrug,g=='UP')
gene_symbol <- rownames(UP)
gene_symbol[c(4,31)] <- c("WASIR1",'NTM')
resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(OV955,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[8] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
#SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(0.2, 1))
SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(0, 1))
library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

##降维图
p2 <- DimPlot(OVmfsub, reduction = 'umap', group.by = "seurat_clusters", label=T)
p1+p2
gene_list=list(gene=marker_10$gene)
#这个打分的含义就是这个亚群在每个细胞上的一个分数，可以进行平行比较
#至于具体有什么生物学意义就仁者见仁，智者见智了


##GSE30161分布
load('~/合并/degdrug30161.Rdata')
##因为对比矩阵是CR对PDPR,所以下调的是铂耐药
DOWN <- subset(degdrug30161,g=='DOWN')
gene_symbol <- rownames(DOWN)
gene_symbol[c(7,9,30,33,34,36)] <- c("LINC01000",'H19','IGF2','TAX1BP3','DUXAP10','NNMT')

resistance_features <- list(gene_symbol)

resistancescore <- AddModuleScore(OV955,
                                  features = resistance_features,
                                  ctrl = 100,
                                  name = "resistance_Features")

colnames(resistancescore@meta.data)
colnames(resistancescore@meta.data)[8] <- 'resistance_Score'

VlnPlot(resistancescore,features = 'resistance_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(resistancescore, features = 'resistance_Score',alpha = c(0, 1))

library(ggplot2)
mydata<- FetchData(resistancescore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


##血管生成
rm(list=ls())
OV955 <- readRDS('OV955.rds')

gene_symbol <- angio$HALLMARK_ANGIOGENESIS

angio_features <- list(gene_symbol)

angioscore <- AddModuleScore(OV955,
                                  features = angio_features,
                                  ctrl = 100,
                                  name = "angio_Features")

colnames(angioscore@meta.data)
colnames(angioscore@meta.data)[8] <- 'angioscore_Score'

VlnPlot(angioscore,features = 'angioscore_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(angioscore, features = 'angioscore_Score',alpha = c(1, 1))

library(ggplot2)
mydata<- FetchData(angioscore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))


##缺氧生成
rm(list=ls())
OV955 <- readRDS('OV955.rds')

gene_symbol <- hypoxia$HALLMARK_HYPOXIA
  
hypoxia_features <- list(gene_symbol)

hypoxiascore <- AddModuleScore(OV955,
                             features = hypoxia_features,
                             ctrl = 100,
                             name = "hypoxia_Features")

colnames(hypoxiascore@meta.data)
colnames(hypoxiascore@meta.data)[8] <- 'hypoxiascore_Score'

VlnPlot(hypoxiascore,features = 'hypoxiascore_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
SpatialFeaturePlot(hypoxiascore, features = 'hypoxiascore_Score',alpha = c(1, 1))

library(ggplot2)
mydata<- FetchData(angioscore,vars = c("UMAP_1","UMAP_2","resistance_Score"))
a <- ggplot(mydata,aes(x = UMAP_1,y =UMAP_2,colour = resistance_Score))+
  geom_point(size = 1)+scale_color_gradientn(values = seq(0,1,0.2),
                                             colours = c('#333366',"#6666FF",'#CC3333','#FFCC33'))
a+ theme_bw() + theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                      panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))

gene <- c('EPCAM',"KRT19",'MUC16')
SpatialFeaturePlot(OV955, features = gene, alpha = c(0.1, 1))
SpatialFeaturePlot(OV955, features = "Ttr", alpha = c(0.1, 1))
SpatialFeaturePlot(OV955, features = "Ttr", alpha = c(0.1, 1))
SpatialFeaturePlot(OV955, features = "Ttr", alpha = c(0.1, 1))

