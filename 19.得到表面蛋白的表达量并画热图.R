rm(list=ls())
library(dplyr)
stromal <- readRDS('stromal.rds')
load('celltypemarkers.Rdata')
library(dplyr)
top30markers <- celltypemarkers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
top50 = celltypemarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
surfacemarker <- c('CTSK','SULF2',##dCAF
                   '','',)

# 获取平均表达量
Idents(stromal) <- "celltype"   # 这一步可以指定要计算哪一个分组的平均表达量，可以选择细胞类型（"CellType"）cluster（"seurat_clusters"）或者是样本类型（"orig.ident")，要注意这里的变量名称不一定正确，要根据数据中的具体变量来指定
AverageExp <- AverageExpression(stromal)
expr <- AverageExp$RNA

Collagens <- c('COL1A1','COL1A2','COL3A1','COL4A1','COL4A2','COL5A1','COL5A2','COL6A1','COL7A1',
               'COL8A1','COL10A1','COL11A1','COL12A1','COL13A1','COL14A1','COL15A1','COL16A1','COL18A1')
ECM <- c('BGN','DCN','LUM','TAGLN','ELN','FN1')
Inflammatory <- c('CFD','CFI','C3','C7','CXCL14','CXCL12','IL33','CXCL3','CXCL2','CXCL1','CCL2','IL6','IL7')
Contractile <- c('ACTA2','MYL6','MYH9','MYH11','TPM1','SORBS2')

Collagens <- c('COL4A1','COL4A2','COL5A1','COL5A2','COL6A1','COL7A1',
               'COL8A1','COL10A1','COL11A1','COL12A1','COL13A1','COL14A1','COL15A1','COL16A1','COL18A1')

Inflammatory <- c('CFD','CFI','C3','C7','CXCL14','CXCL12','IL33','CXCL3','CXCL2','CXCL1','CCL2','IL6','IL7')
Contractile <- c('ACTA2','MYL6','MYH9','MYH11','TPM1','SORBS2')

all <- c(Collagens,Inflammatory,Contractile)
which(all %in% rownames(expr)=='FALSE')
all[c(29,43)]

allexpr <- expr[all,]
allexpr <- allexpr[,c(1,3,5)]
colnames(allexpr)[1] <- 'dCAF'
colnames(allexpr)[2] <- 'iCAF'
colnames(allexpr)[3] <- 'mCAF'
library(ComplexHeatmap)
##先进行标准化
data <- allexpr
exp <- apply(data, 2, scale)##1行2列
rownames(exp) <- rownames(data)

Heatmap(exp,
        cluster_columns = F,
        cluster_rows = F,
        border = 'black',
        rect_gp = gpar(col = "white", lwd = 2),
        column_names_rot = 45)
