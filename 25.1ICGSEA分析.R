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
table(celltypemarkers$cluster)
ICmarkers <- rownames(subset(celltypemarkers,cluster=='inflammatory CAFs'))
geneset <- data.frame(c(rep('IC',381)),ICmarkers)
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
gseaplot2(egmt,geneSetID = 'IC',pvalue_table=T)
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
gseaplot2(egmt,geneSetID = 'IC',pvalue_table=T)
ridgeplot(egmt)