rm(list=ls())
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
##提取肿瘤细胞进行分析
CAF <- subset(stromal,celltype=='demosplastic CAFs'|celltype=='inflammatory CAFs'|celltype=='myofibroblast CAFs')

cell<- names(sort(CAF@assays$RNA@counts['SFRP2',],
                  decreasing =T))
SFRP2expr <- c(rep('high',3668),rep('low',3668)) 
newmeta<-data.frame(cell,SFRP2expr)

##名字对齐
rownames(newmeta) <- newmeta[,1]
cellname <- rownames(CAF@meta.data)
newmeta <- newmeta[cellname,]

##加回去
CAF <- AddMetaData(CAF,metadata = newmeta[,2],col.name = 'SFRP2expr')
Idents(CAF) <- 'SFRP2expr'
p1 <- VlnPlot(CAF,features = 'SFRP2')
p1

#比较B_cell和T_cells的差异表达基因
Idents(CAF) <- 'SFRP2expr'
dge.cellSFRP2p <- FindMarkers(CAF, ident.1 = 'high', ident.2 = 'low',logfc.threshold = 0.05,group.by = 'SFRP2expr')
save(dge.cellSFRP2p,file='dge.cellSFRP2p.Rdata')

##
markers <- rownames(dge.cellSFRP2p)

geneID <- bitr(markers, fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)#进行id转换，可在symbol,ENTREZID,SYMBOL中转换
geneID <- geneID$ENTREZID

geneIDall <- bitr(rownames(stromal), fromType = "SYMBOL",
                  toType = c( "ENTREZID"),
                  OrgDb = org.Hs.eg.db)#进行id转换，可在symbol,ENTREZID,SYMBOL中转换
geneIDall <- geneIDall$ENTREZID
ego_ALL <- enrichGO(gene          = markers,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)

ego_CC <- enrichGO(gene          = markers,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = markers ,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = markers,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

ego_BPresultDC <- ego_BP@result
p_BP <- barplot(ego_BP,showCategory = 5) + ggtitle("barplot for BiologDCal process")
p_CC <- barplot(ego_CC,showCategory = 5) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 5) + ggtitle("barplot for Molecular function")
ego_CCresult <- ego_CC@result
ego_MFresult <- ego_MF@result
ego_BPresult <- ego_BP@result
plotc <- p_BP/p_CC/p_MF
plotc

R.utils::setOption( "clusterProfiler.download.method",'auto' )
kk <- enrichKEGG(gene         = geneID,
                 organism     = 'hsa',
                 universe = geneIDall ,
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.1)##必须得挂VPN

head(kk)[,1:6]
KKSFRP2result <- kk@result
pathway <- KKSFRP2result[c("hsa04141","hsa04510","hsa03050","hsa04218"),]
library(DOSE)
pathway$GeneRatio <- parse_ratio(pathway$GeneRatio)
pathway$GeneRatio <- round(pathway$GeneRatio,digits = 2)
head(pathway)
ggplot(pathway,aes(GeneRatio,Description))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low="green",high = "red")+
  theme(axis.text.y=element_text(vjust=1,size=60,face = "bold"))+#
  labs(color=expression(pvalue),size="Count", ##expression函数定义函数样式 []添加下标，^添加上标
       x="GeneRatio", ##自定义标轴
       y="",
       title="KEGG Pathway Enrichment")+ ##自定义坐标轴
  theme_bw() #去掉背景
