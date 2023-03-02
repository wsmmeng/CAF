##表面marker要画小提琴图,用CAF
##提取细胞子集
##参考《美化你的单细胞各个亚群特异性高表达基因小提琴图》生信技能树
Cells.sub <- subset(stromal@meta.data, celltype=="demosplastic CAFs"|celltype=="inflammatory CAFs"|celltype=='myofibroblast CAFs')
CAF <- subset(stromal, cells=row.names(Cells.sub))
##这个函数最上面的标题有问题，需要自己加
ViolinPlot <- function(object, groupBy, MarkerSelected) {
  # (1)获取绘图数据1
  plot_data = FetchData(object = object, 
                        vars = c(MarkerSelected$gene, groupBy), 
                        slot = 'data') %>% 
    dplyr::rename(group = as.name(groupBy)) %>% 
    tidyr::pivot_longer(cols = -group, names_to = 'Feat', values_to = 'Expr')
  
  # (2)获取绘图数据2
  ident_plot = MarkerSelected %>% 
    dplyr::select(cluster, gene)
  
  # (3)绘图
  figure_1 = ggplot(data = plot_data, mapping = aes(x = Expr,
                                                    y = fct_rev(factor(x = Feat, 
                                                                       levels = MarkerSelected$gene)), 
                                                    fill = group, 
                                                    label = group)) +
    geom_violin(scale = 'width', adjust = 1, trim = TRUE) +
    scale_x_continuous(expand = c(0, 0), labels = function(x)
      c(rep(x = '', times = length(x) - 2), x[length(x) - 1], '')) +
    facet_grid(cols = vars(group), scales = 'free') +
    cowplot::theme_cowplot(font_family = 'Arial') +
    scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3')) +
    xlab('Expression Level') + 
    ylab('') +
    theme(legend.position = 'none', 
          panel.spacing = unit(x = 0, units = 'lines'),
          axis.line = element_blank(), #去除x和y轴坐标线(不包括axis tick)；
          panel.background = element_rect(fill = NA, color = 'black'),
          strip.background = element_blank(), #去除分页题头背景；
          strip.text = element_text(color = 'black', size = 20),
          axis.text.x = element_text(color = 'black', size = 11),
          axis.text.y = element_blank(),
          axis.title.x = element_text(color = 'black', size = 15),
          axis.ticks.x = element_line(color = 'black'),
          axis.ticks.y = element_blank(),
          axis.ticks.length = unit(x = 0.1, units = 'cm'))
  
  figure_2 = ggplot(data = ident_plot, aes(x = 1,
                                           y = fct_rev(factor(x = gene, levels = MarkerSelected$gene)),
                                           fill = cluster)) +
    geom_tile() +
    theme_bw(base_size = 12) +
    scale_fill_manual(values = paletteer::paletteer_d('ggsci::category20c_d3')) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    guides(fill = guide_legend(direction = 'vertical',
                               label.position = 'right',
                               title.theme = element_blank(),
                               keyheight = 0.5,
                               nrow = 2)) +
    xlab('Feature') +
    theme(legend.text = element_text(color = 'black', size = 11),
          legend.position = 'bottom',
          legend.justification = 'left',
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,05,0,0),
          panel.spacing = unit(0, 'lines'),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          plot.margin = unit(x = c(0,0,0,0), units = 'cm'),
          axis.title.y = element_blank(),
          axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, color = 'black'),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  figure_2 + figure_1 + patchwork::plot_layout(nrow = 1, widths = c(0.03, 0.97))
}
load('celltypemarkers.Rdata')

marker <- c('CTSK','SULF1','RGS5','APOE')
library(patchwork)
marker_selected_1 <- celltypemarkers[marker,]
ViolinPlot(object = CAF, groupBy = 'celltype', MarkerSelected = marker_selected_1)
