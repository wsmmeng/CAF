rm(list=ls())

###提取对象中的count数并合并
OVPPP <- readRDS('OVPPP.rds')
OVI <- readRDS('OVI.rds')
OVP <- readRDS('OVP.rds')
OVH <- readRDS('OVH.rds')

##注意样本的细胞数目
OVPPP <- as.data.frame(OVPPP[["RNA"]]@counts)#
OVI <- as.data.frame(OVI[["RNA"]]@counts)#
OVP <- as.data.frame(OVP[["RNA"]]@counts)#
OVH <- as.data.frame(OVH[["RNA"]]@counts)#

OVPPPgene <- as.data.frame(rownames(OVPPP))
OVIgene  <- as.data.frame(rownames(OVI))
OVPgene <- as.data.frame(rownames(OVP))
OVHgene <- as.data.frame(rownames(OVH))

{OVPPPgene <- rownames(OVPPP)
  OVIgene <- rownames(OVI)
  OVPgene <- rownames(OVP)
  OVHgene <- rownames(OVH)
 
  row <- intersect(OVPPPgene,OVIgene) %>% intersect(OVPgene) %>% intersect(OVHgene)
  OVPPP<- OVPPP[row,]
  OVI <- OVI[row,]
  OVP <- OVP[row,]
  OVH <- OVH[row,]
}
colnames(OVPPP)[1:5]
colnames(OVI)[1:5]
colnames(OVP)[1:5]
colnames(OVH)[1:5]
OV <- CreateSeuratObject(counts = cbind(OVPPP,OVI,OVP,OVH), project = "OVall", min.cells = 5)
OV@meta.data$group <- c(rep("OVPPP",21761), rep("OVI", 30025),rep("OVP", 9885),rep("OVH", 15202))#赋值条件变量
saveRDS(OV,'OV.rds')

