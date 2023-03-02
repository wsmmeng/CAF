setwd('~/合并/')
CAF <- subset(stromal,celltype=='demosplastic CAFs'|celltype=='inflammatory CAFs'|celltype=='myofibroblast CAFs')
##对三个基因可视化
mycolor <- c('gray', 'white','red')#设置颜色 
#theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框FeaturePlot(OV, features = c("ACTA2"),cols = mycolor, pt.size = 1.5,max.cutoff = 1)  
FeaturePlot(OV, features = c("PRRX1"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
FeaturePlot(stromal, features = c("PRRX1"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
VlnPlot(stromal, features = c("PRRX1"))

FeaturePlot(OV, features = c("SFRP2"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
FeaturePlot(stromal, features = c("SFRP2"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
VlnPlot(stromal, features = c("SFRP2"))

FeaturePlot(OV, features = c("POSTN"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
FeaturePlot(stromal, features = c("POSTN"),cols = mycolor, pt.size = 0.5,max.cutoff = 1)
VlnPlot(stromal, features = c("POSTN"))


