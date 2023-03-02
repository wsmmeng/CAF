library(CellChat)
library(Seurat)

##参考《CellChat学习笔记【一】——通讯网络构建》
##准备输入文件
data <- GetAssayData(object = OV, slot = 'data')
meta <- OV@meta.data
##可以unique一下去掉无用的因子，不然后面会报错
meta$celltype2 = droplevels(meta$celltype2, exclude = setdiff(levels(meta$celltype),unique(meta$celltype)))
table(meta$celltype2)

cellchat <- createCellChat(object = data,
                           meta = meta,
                           group.by = 'celltype2')
##选择数据库
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB

##选择第一个secreted signalling数据库，也可以不选，我这里就不选
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# 将数据库传递给cellchat对象
#cellchat@DB <- CellChatDB.use 
library(future)
plan(strategy = 'multiprocess', workers = 14)
#subset
cellchat <- subsetData(cellchat)#为什么这里要做subset，因为为了节省运算时间和空间，通过这一步，我们后面的分析将只关注于与细胞通讯有关的基因（从数据框中提取出来的），所以在这里针对于表达矩阵的基因取了子集。

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#根据表达值推测细胞互作的概率


cellchat <- computeCommunProb(cellchat)

##过滤细胞数目比较少的细胞类型
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
#推断信号通路水平的细胞通讯网络（结果储存在@netP下面，有一个概率值和对应的pval）
cellchat <- computeCommunProbPathway(cellchat)

df.netp <- subsetCommunication(cellchat, slot.name = "netP")
save(df.net,df.netp,file='net.Rdata')
save(cellchat,file='cellchat.Rdata')

####构建细胞通讯网络
load('cellchat.Rdata')
cellchat <- aggregateNet(cellchat)
##气泡图(全部受体和配体)
levels(cellchat@idents)
cellchat@netP$pathways
###气泡图
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
##DC
p = netVisual_bubble(cellchat, sources.use = c(2),
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12), remove.isolate = FALSE)
p
##IC
p2 = netVisual_bubble(cellchat, sources.use = c(4),
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12), remove.isolate = FALSE)
p2
##MC
p3 = netVisual_bubble(cellchat, sources.use = c(6),
                      targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12), remove.isolate = FALSE)
p3

##和弦图
netVisual_chord_gene(cellchat,
                     sources.use = 2,
                     targets.use = c(1:12),
                     lab.cex = 0.5,legend.pos.y = 30)
