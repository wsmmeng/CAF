setwd('~/合并')
rm(list=ls())
library(Seurat)
library(future)
plan()
plan('multiprocess',workers=10)
install.packages('/home/shpc_101561/monocle_2.26.0.tar.gz',
                 repos = NULL, type='source',
                 lib="/home/xiyou/rpackage/")

.libPaths(c("/home/shpc_101561/rpackage", .libPaths()))

##参考《单细胞转录组基础分析六：伪时间分析》
##关于ordercell的解决方法，参考《解决monocle2的ordercell报错的两种方法》《https://www.jianshu.com/p/d6644a8cb404》
##
library(monocle)
library(Seurat)
OV <- readRDS('OV.rds')
table(OV@meta.data$celltype2)
##提取细胞子集
Cells.sub <- subset(OV@meta.data, celltype2=="demosplastic CAFs"|celltype2=="inflammatory CAFs"|celltype2=='myofibroblast CAFs')
CAF <- subset(OV, cells=row.names(Cells.sub))

##数据导入，创建对象
dir.create("monocle")
data <- as(as.matrix(CAF@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CAF@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
##做一些中间处理
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)


##我一般用第三种方法
##使用clusters差异表达基因
##寻找差异表达的基因
diffCAFr = FindAllMarkers(OVrfsub)
diff.genes <- subset(diffCAFr,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
##使用seurat选择的高变基因
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
p3
p1|p2|p3


#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
trace('project2MST',edit = T,where = asNamespace('monocle'))
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("pseudotime/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/State.jpg", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("pseudotime/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("pseudotime/Cluster.jpg", plot = plot2, width = 6, height = 5)
##celltype轨迹分布图
plot3 <- plot_cell_trajectory(mycds, color_by = "celltype2")
ggsave("pseudotime/celltype.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Celltype.jpg", plot = plot3, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Pseudotime.jpg", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("pseudotime/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("pseudotime/Combination.jpg", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime/pseudotime.csv")


#cluster差异基因
diffCAFr = FindAllMarkers(OVrfsub)
library(dplyr)
top20 = diffCAFr %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#sig_diff.genes <- subset(diffCAFr,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene
sig_diff.genes <- top20$gene
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=3,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap1.jpg", plot = p1, width = 5, height = 8)
#高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap2.jpg", plot = p2, width = 5, height = 10)