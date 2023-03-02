immune <- readRDS()
stromal <- readRDS("stromal.rds")
OV <- readRDS('OV.rds')

metaOV <- OV@meta.data
metaimmune <- immune@meta.data
metastromal <- stromal@meta.data

celltypeOV <- data.frame(metaOV[,c('celltype')])
rownames(celltypeOV)<- rownames(metaOV)
colnames(celltypeOV) <- 'celltype'
table(celltypeOV$celltype)

celltypeep <- subset(celltypeOV,celltype=='epithelial tumor cell')

celltypestromal <- data.frame(metastromal[,c('celltype')])
rownames(celltypestromal)<- rownames(metastromal)
colnames(celltypestromal)<- 'celltype'

celltypeimmune <- data.frame(metaimmune[,c('celltype')])
rownames(celltypeimmune)<- rownames(metaimmune)
colnames(celltypeimmune)<- 'celltype'

identical(rownames(celltypeall),rownames(celltypeOV))
celltypeall <- rbind(celltypeep,celltypestromal,celltypeimmune)
celltypeall$XX <- c(rep('xx',67201))##随便加一行,不然重排序的时候会出错
celltypeall <- celltypeall[rownames(metaOV),]
celltype2 <- celltypeall$celltype
celltype2 <- data.frame(celltype2)
rownames(celltype2) <- rownames(celltypeall)

OV <- AddMetaData(OV, celltype2)
table(OV@meta.data$celltype2)
metaOV <- OV@meta.data##
OV@meta.data$celltype2 = droplevels(OV@meta.data$celltype2, exclude = setdiff(levels(OV@meta.data$celltype2),unique(OV@meta.data$celltype2)))
saveRDS(OV,'OV.rds')
