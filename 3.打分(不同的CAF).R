setwd('~/合并/spatial/123467/A_955_OvarianTumor/')
load('~/合并/celltypemarkers.Rdata')


OV955 <- readRDS('OV955.rds')

##GSE63885分布
##因为对比矩阵是Resistent对sensitive,所以上调的是铂耐药
des <- subset(celltypemarkers,cluster=='demosplastic CAFs')
gene_symbol <- rownames(des)
des_features <- list(gene_symbol)
desscore <- AddModuleScore(OV955,
                                  features = des_features,
                                  ctrl = 100,
                                  name = "des_Features")
colnames(desscore@meta.data)
colnames(desscore@meta.data)[8] <- 'des_Score'

VlnPlot(desscore,features = 'des_Score',pt.size = 0, adjust = 2,group.by = "seurat_clusters")
#SpatialFeaturePlot(desscore, features = 'des_Score',alpha = c(0.2, 1))

SpatialFeaturePlot(desscore, features = 'des_Score',alpha = c(0, 1))
SpatialFeaturePlot(desscore, features = "PRRX1",alpha = c(0, 1))
SpatialFeaturePlot(desscore, features = "SFRP2",alpha = c(0, 1))
SpatialFeaturePlot(desscore, features = "POSTN",alpha = c(0, 1))

inf <- subset(celltypemarkers,cluster=='inflammatory CAFs')
gene_symbol <- rownames(inf)
inf_features <- list(gene_symbol)
infscore <- AddModuleScore(OV955,
                           features = inf_features,
                           ctrl = 100,
                           name = "inf_Features")
colnames(infscore@meta.data)
colnames(infscore@meta.data)[8] <- 'inf_Score'
#SpatialFeaturePlot(desscore, features = 'des_Score',alpha = c(0.2, 1))
SpatialFeaturePlot(infscore, features = 'inf_Score',alpha = c(0, 1))


myf <- subset(celltypemarkers,cluster=='myofibroblast CAFs')
gene_symbol <- rownames(myf)
myf_features <- list(gene_symbol)
myfscore <- AddModuleScore(OV955,
                           features = myf_features,
                           ctrl = 100,
                           name = "myf_Features")
colnames(myfscore@meta.data)
colnames(myfscore@meta.data)[8] <- 'myf_Score'
#SpatialFeaturePlot(desscore, features = 'des_Score',alpha = c(0.2, 1))
SpatialFeaturePlot(myfscore, features = 'myf_Score',alpha = c(0, 1))





