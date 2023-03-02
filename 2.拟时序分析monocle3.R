rm(list=ls())
library(Seurat)

library(monocle3)

library(tidyverse)

library(patchwork)

rm(list=ls())

dir.create("Monocle3")

setwd('~/CAF整合')



#//创建CDS对象并预处理数据
OVfall <- readRDS("OVfall.rds")
data <- GetAssayData(OVfall, assay ='RNA', slot = 'counts')
cell_metadata <- OVfall@meta.data
gene_annotation <-data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <-rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata =cell_metadata, 
                         gene_metadata =gene_annotation)

#1. 标准化和PCA降维

#（RNA-seq是使用PCA，如果是处理ATAC-seq的数据用Latent Semantic Indexing)
#//preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)

#2. 可视化
#umap降维
cds <- reduce_dimension(cds, reduction_method="UMAP",preprocess_method = "PCA")
colData(cds)
plot_cells(cds,color_cells_by="celltype",cell_size = 2,graph_label_size = 5,group_label_size= 3)
cds <- cluster_cells(cds) 
#轨迹学习Learn the trajectory graph（使用learn_graph()函数）
## 识别轨迹
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=FALSE, label_branch_points=FALSE)
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster=FALSE,label_leaves=TRUE, label_branch_points=TRUE,graph_label_size=3)
##可视化基因
markers <- c('SFRP2','POSTN','PRRX1')
p <- plot_cells(cds,genes=markers, show_trajectory_graph=FALSE,
                
                label_cell_groups=FALSE,  label_leaves=FALSE)
p
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,  label_branch_points = FALSE)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T,label_branch_points = T)
saveRDS(cds, file = "cds.rds")


#寻找拟时轨迹差异基因

#//graph_test分析最重要的结果是莫兰指数（morans_I），其值在-1至1之间，0代表此基因没有
##//空间共表达效应，1代表此基因在空间距离相近的细胞中表达值高度相似。

Track_genes <- graph_test(cds,neighbor_graph="principal_graph", cores=6)
Track_genes <-Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes,"Trajectory_genes.csv", row.names = F)

##挑选top10画图展示
Track_genes_sig <- Track_genes %>%top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()


##基因表达趋势图
plot_genes_in_pseudotime(cds[Track_genes_sig,],color_cells_by="celltype",min_expr=0.5, ncol= 2)

#FeaturePlot图
p <- plot_cells(cds,genes=Track_genes_sig, show_trajectory_graph=FALSE,
                
                label_cell_groups=FALSE,  label_leaves=FALSE)
p$facet$params$ncol <- 5
p

ggsave("Genes_Featureplot.pdf",plot = p, width = 20, height = 8)


#寻找共表达基因模块
Track_genes <-read.csv("Trajectory_genes.csv")
genelist <- pull(Track_genes,gene_short_name) %>% as.character()
gene_module <-find_gene_modules(cds[genelist,], resolution=1e-1, cores = 6)
write.csv(gene_module,"Genes_Module.csv", row.names = F)
cell_group <-tibble::tibble(cell=row.names(colData(cds)),
                            cell_group=colData(cds)$seurat_clusters)

agg_mat <-aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <-stringr::str_c("Module ", row.names(agg_mat))
p <- pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")
ggsave("Genes_Module.pdf", plot =p, width = 8, height = 8)


#提取拟时分析结果返回seurat对象
pseudotime <- pseudotime(cds,reduction_method = 'UMAP')
pseudotime <-pseudotime[rownames(OVfall@meta.data)]
OVfall$pseudotime <- pseudotime
p = FeaturePlot(OVfall, reduction ="umap", features = "pseudotime")
p
ggsave("Pseudotime_Seurat.pdf",plot = p, width = 8, height = 6)
saveRDS(OVfall, file ="sco_pseudotime.rds")
