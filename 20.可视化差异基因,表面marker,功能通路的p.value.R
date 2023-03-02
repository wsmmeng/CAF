rm(list = ls (all = TRUE))
library(monocle)
library(Seurat)
OVfall <- readRDS('OVfall.rds')
# 获取平均表达量
Idents(OVfall) <- "celltype"   # 这一步可以指定要计算哪一个分组的平均表达量，可以选择细胞类型（"CellType"）cluster（"seurat_clusters"）或者是样本类型（"orig.ident")，要注意这里的变量名称不一定正确，要根据数据中的具体变量来指定
AverageExp <- AverageExpression(OVfall)
expr <- AverageExp$RNA

