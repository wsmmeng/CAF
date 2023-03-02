##
rm(list=ls())
stromal <- readRDS('stromal.rds')
load('celltypemarkers.Rdata')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

table(celltypemarkers$cluster)
ICmarker <- subset(celltypemarkers,cluster=='inflammatory CAFs')
ICgene <- ICmarker$gene
gene <- ICmarker$gene
geneID <- bitr(gene, fromType = "SYMBOL",
               toType = c( "ENTREZID"),
               OrgDb = org.Hs.eg.db)#进行id转换，可在symbol,ENTREZID,SYMBOL中转换
geneID <- geneID$ENTREZID

geneIDall <- bitr(rownames(stromal), fromType = "SYMBOL",
                  toType = c( "ENTREZID"),
                  OrgDb = org.Hs.eg.db)#进行id转换，可在symbol,ENTREZID,SYMBOL中转换
geneIDall <- geneIDall$ENTREZID
ego_ALL <- enrichGO(gene          = gene,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)

ego_CC <- enrichGO(gene          = gene,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = gene ,
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = gene,
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
p_BP <- barplot(ego_BP,showCategory = 5) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 5) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 5) + ggtitle("barplot for Molecular function")
ego_CCresult <- ego_CC@result
ego_MFresult <- ego_MF@result
ego_BPresult <- ego_BP@result
plotc <- p_BP/p_CC/p_MF
plotc
ego_BPresultIC <- ego_BP@result
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kk <- enrichKEGG(gene         = geneID,
                 organism     = 'hsa',
                 universe = geneIDall ,
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.1)##必须得挂VPN

head(kk)[,1:6]
KKICresult <- kk@result
pathway <- KKICresult[c("hsa04657","hsa04668","hsa04210","hsa04064","hsa04218"),]
library(DOSE)
pathway$GeneRatio <- parse_ratio(pathway$GeneRatio)
pathway$GeneRatio <- round(pathway$GeneRatio,digits = 2)
head(pathway)
ggplot(pathway,aes(GeneRatio,Description))+
  geom_point(aes(size=Count,color=pvalue))+
  scale_color_gradient(low="green",high = "red")+
  theme(axis.text.y=element_text(vjust=1,size=40,face = "bold"))+#
  labs(color=expression(pvalue),size="Count", ##expression函数定义函数样式 []添加下标，^添加上标
       x="GeneRatio", ##自定义标轴
       y="",
       title="KEGG Pathway enrichment")+ ##自定义坐标轴
  theme_bw() #去掉背景

save(ego_BPresultIC,file = 'ego_BHresultIC.Rdata')
###用p.adjust画功能分析的热图,这个需要整合
p.adjust <- ego_BHresult[,c(1,2,3,4,5)]
p.adjustXX <- p.adjust$p.adjust
p.adjustXX <- -log10(p.adjustXX)
p <- data.frame(p.adjustXX1,p.adjustXX2)
rownames(p) <- rownames(p.adjust)

