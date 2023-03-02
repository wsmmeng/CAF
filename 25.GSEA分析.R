##GSEA分析
 ##得到基因的表达排序
load('degdrug.Rdata')
geneList= degdrug$logFC
names(geneList)= toupper(rownames(degdrug))
geneList=sort(geneList,decreasing = T)
head(geneList)
library(GSEABase) # BiocManager::install('GSEABase')

##得到geneset
load('celltypemarkers.Rdata')
DCmarkers <- rownames(subset(celltypemarkers,cluster=='demosplastic CAFs'))
geneset <- data.frame(c(rep('DC',399)),DCmarkers)
colnames(geneset)[1] <- 'term'
length(unique(geneset$term))
geneList = degdrug$logFC
names(geneList)=toupper(rownames(degdrug))
geneList=sort(geneList,decreasing = T)
head(geneList)
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'gsea_results_df.csv')
library(enrichplot)
gseaplot2(egmt,geneSetID = 'DC',pvalue_table=T)
ridgeplot(egmt) 

##GSE30161
load('degdrug30161.Rdata')
geneList= degdrug30161$logFC
names(geneList)= toupper(rownames(degdrug30161))
geneList=sort(geneList,decreasing = T)
head(geneList)
library(GSEABase) # BiocManager::install('GSEABase')

##得到geneset
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'gsea_results_df.csv')
library(enrichplot)
gseaplot2(egmt,geneSetID = 'DC',pvalue_table=T)
ridgeplot(egmt) 


##这个不是耐药，先不分析了
load('degdrug23554.Rdata')
geneset <- c(rep('term',),marker)  
length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'gsea_results_df.csv')
library(enrichplot)
gseaplot2(egmt,geneSetID = 'resistance',pvalue_table=T)
ridgeplot(egmt) 
