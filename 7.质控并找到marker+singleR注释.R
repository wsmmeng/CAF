setwd('~/合并/')
rm(list=ls())
library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(mindr))install.packages("tidyverse")
library(ggplot2)
OV <- readRDS('OV.rds')
##计算质控指标
#计算细胞中线粒体基因比例
OV[["percent.mt"]] <- PercentageFeatureSet(OV, pattern = "^MT-")
#无
#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(OV@assays$RNA)) 
HB.genes <- rownames(OV@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
OV[["percent.HB"]]<-PercentageFeatureSet(OV, features=HB.genes) 
#无
#计算ERCC比例
OV[["percent.ERCC"]] <- PercentageFeatureSet(OV, pattern = "^ERCC-")
#无
#head(OV@meta.data)
#观察count,feature,mt,hb,ercc的值
Idents(OV) <- 'orig.ident'
col.num <- length(levels(OV@active.ident))
violin <- VlnPlot(OV,
                  features = c("nFeature_RNA", "nCount_RNA","percent.mt","percent.HB","percent.ERCC"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 3) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin <- VlnPlot(OV,
                  features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #不需要显示点，可以设置pt.size = 0
                  ncol = 3
                  ) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 12, height = 6) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 12, height = 6)  
plot1 <- FeatureScatter(OV, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OV, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(OV, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot1+plot2

name <- data.frame(rownames(OV))

##设置质控标准
print(c("请输入允许基因数和核糖体比例，示例如下：", "minGene=500", "maxGene=4000", "pctMT=20"))
minGene=200

##ncount基因小于25000？
##数据质控
OV <- subset(OV, subset = nFeature_RNA > 200 & nCount_RNA > 600 & percent.mt < 20)#还可以加
dim(OV)
col.num <- length(levels(OV@active.ident))
violin <- VlnPlot(OV,
                 features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                 cols =rainbow(col.num), 
                 pt.size = 0.1, 
                 ncol = 3) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
ggsave("质控标准/vlnplot_after_qc.pdf", plot = violin, width = 12, height = 6) 

##创建一个标准的seurat对象
#一
library(future)
options(future.globals.maxSize=4000000000)
plan('multiprocess',workers=16)
plan()
OV <- NormalizeData(OV, normalization.method = "LogNormalize", scale.factor = 10000)
#二
plan('multiprocess',workers=16)
OV <- FindVariableFeatures(OV,selection.method = "vst", nfeatures = 2000)
#三
scale.genes <-  rownames(OV)
plan('multiprocess',workers=16)

OV <- ScaleData(OV, features = scale.genes)
#看需不需要消除细胞周期的影响
saveRDS(OV,'OV.rds')
#查看我们选择的高变基因中有哪些细胞周期相关基因：
CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(OV))
#细胞周期评分
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(OV))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(OV))
plan('multiprocess',workers=16)
OV <- CellCycleScoring(object=OV,g2m.features=g2m_genes,s.features=s_genes)
#查看细胞周期基因对细胞聚类的影响
OVa <- RunPCA(OV, features = c(s_genes, g2m_genes))
p <- DimPlot(OVa, reduction = "pca", group.by = "Phase")
p
#感觉不太好
##如果需要消除细胞周期的影响
plan('multiprocess',workers=1)
OV$CC.difference <- OV$S.Score-OV$G2M.Score
OV <- ScaleData(OV, vars.to.regress = 'CC.difference', features = rownames(OV))
OV <- SCTransform(OV,vars.to.regress = c('S.Score','G2M.Score'))
ggsave("cluster/cellcycle_pca.png", p, width = 8, height = 6)
#四
OV <- readRDS('OV.rds')
plan('multiprocess',workers=10)
OV <- RunPCA(OV,features = VariableFeatures(OV))
##作为后续选择维度的依据
plot1 <- DimPlot(OV, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(OV, ndims=50, reduction="pca") 
plot1+plot2
saveRDS(OV,'OV.rds')
##用harmony去除批次效应
##在meta.data加上一列
#赋值条件变量
##看没去批次效应之前的情况
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = OV, reduction = "pca", pt.size = .1,group.by = "group")
p2 <- VlnPlot(object = OV, features = "PC_1", group.by = "group",pt.size = .1)
library(cowplot)
plot_grid(p1,p2)

##harmony去批次效应
library(harmony)
options(repr.plot.height = 2.5, repr.plot.width = 6)
OV <- OV %>%
  RunHarmony("group", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(OV, 'harmony')
harmony_embeddings[1:5, 1:5]
##看去批次效应之后的情况
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = OV, reduction = "harmony", pt.size = .1, group.by = "group")
p2 <- VlnPlot(object = OV, features = "harmony_1", group.by = "group",pt.size = .1)
plot_grid(p1,p2)

saveRDS(OV,'OV.rds')

##降维
##感觉用harmony不好,还是用pca
plan()
OV <- FindNeighbors(OV,reduction = "harmony",dims = 1:50) 
OV <- FindClusters(OV,reduction = "harmony",resolution = 0.5)
table(OV@meta.data$seurat_clusters)
metadata <- OV@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

OV <- RunUMAP(OV,reduction = "pca", dims = 1:50)
OV <- RunTSNE(OV,reduction = "pca", dims = 1:50)

plot2 = DimPlot(OV, reduction = "umap")
plot1 = DimPlot(OV, reduction = "tsne")
plot1+plot2
saveRDS(OV,'OV.rds')

rm(list=ls())
OV <- readRDS('OV.rds')
row
##注释
diff.wilcox = FindAllMarkers(OV)##默认使用wilcox方法挑选差异基因，大概4-5min
all.markers = diff.wilcox %>% select(gene, everything()) %>%
  subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
#An adjusted P value < 0.05and | log 2 [fold change (FC)] | > 0.5 
#were considered the 2 cutoff criteria for identifying marker genes.
save(diff.wilcox,file = "allmarkerovariancluster.Rdata")
load("allmarkerovariancluster.Rdata")

#SingleR鉴定细胞类型
if(!require(SingleR))BiocManager::install(SingleR)
if(!require(matrixStats))BiocManager::install('matrixStats')

if(!require(celldex))BiocManager::install('celldex')#有一些版本冲突，需要重新安装一些包。
memory.limit(100000)

load('E:/R语言代码(Jimmy代码)/原版代码(单细胞jimmy)/原版代码/rawdata/HumanPrimaryCellAtlasData.Rdata')

refdata <- ImmGenData()
#refdata <- MonacoImmuneData()
refdata <- HumanPrimaryCellAtlasData
testdata <- GetAssayData(OV, slot="data")
clusters <- OV@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltypeSR = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
OV@meta.data$celltypeSR = "NA"
#OV@meta.data= OV@meta.data[,-c(13,14)]#就这么调整metadata
for(i in 1:nrow(celltypeSR)){
  OV@meta.data[which(OV@meta.data$seurat_clusters == celltypeSR$ClusterID[i]),'celltypeSR'] <- celltypeSR$celltypeSR[i]}
table(OV@meta.data$celltypeSR)
#结果展示
p1 = DimPlot(OV, group.by="seurat_clusters", label=T, label.size=5, reduction='tsne')
p2 = DimPlot(OV, group.by="celltype", label=T, label.size=5, reduction='umap')


saveRDS(OV, file="OV.rds")
