require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)

#make list of dfs
cpm.SL.figs4 = list("Hyb2HP"=cpm.SL[,c("Gene", "l2.B73xMo17.HP", "l10.ttest.HP")],
                     "Hyb2MP"=cpm.SL[,c("Gene", "l2.B73xMo17.MP", "l10.ttest.MP")],
                     "Par2Par"=cpm.SL[,c("Gene", "l2.B73.Mo17", "l10.ttest.P")])
#name plot panels
cpm.SL.figs4$Hyb2HP$plot = "A"
cpm.SL.figs4$Hyb2MP$plot = "B"
cpm.SL.figs4$Par2Par$plot = "C"
#rename column to match
cpm.SL.figs4 = lapply(cpm.SL.figs4, setNames, c("Gene","l2.ratio", "l10.ttest","plot"))
#bind
cpm.SL.figs4 = do.call("rbind",cpm.SL.figs4)
plot=unique(cpm.SL.figs4[,"plot"])
plotList=vector(mode = "list", length = length(plot))
for(i in 1:length(plot)){
  figS4=cpm.SL.figs4[which(cpm.SL.figs4[,"plot"]==plot[i]),]
  ###Set plot axis limits
  figS4$l2.ratio[figS4$l2.ratio <= -1] = -1
  figS4$l2.ratio[figS4$l2.ratio >= 1] = 1
  ###Set GOIs
  figS4 = figS4
  figS4$Category = "Other"
  figS4[(figS4$Gene %in% GOIs$plastid),"Category"] = "Plastid-localized"
  figS4[(figS4$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
  figS4[(figS4$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
  ###factor
  figS4$Category = factor(figS4$Category, levels = c("Other", "Plastid-localized", "PhANGs", "NE plastid ribosome"))
  figS4 = figS4[(order(figS4$Category)),]
  #colorblind-friendly palette
  colors.p <- c("Other"="#D8D8D8", "Plastid-localized"="black",
                "NE plastid ribosome"="#E5A3FF",
                "PhANGs"="#81E401")
  #scatter plot with marginal density curves
  figS4.scat = ggplot(figS4, mapping=aes(x=l2.ratio, y=l10.ttest, color=Category, shape=Category))+
    geom_point(aes(color=Category), size=1.5)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_shape_manual(values=c(1,19,19,19), guide="none")+
    theme_classic()+
    xlim(c(-1,1))+
    geom_vline(xintercept=0, linetype="dashed", color = "gray40", size=0.5)+
    geom_hline(yintercept= 1.30103, linetype="dashed", color = "gray40", size=0.5)+
    ylab("-log10(p-value)")+
    ggtitle(plot[i])+
    theme(aspect.ratio=1)+
    coord_fixed(expand = F)+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
    theme(axis.text =element_text(size=10), legend.text = element_text(size=10), axis.title = element_text(size=10),
          legend.position="bottom", legend.title = element_blank(), plot.title = element_text(size=10, vjust = -3, hjust = -0.2), legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(-3, 'cm'))+
    guides(col=guide_legend(ncol = 2))
  if(i==1){
    figS4.scat = figS4.scat+ xlab("log2 BxM/HP")
  }
  else if(i==2){
    figS4.scat = figS4.scat+ xlab("log2 BxM/MP")
  }
  else if(i==3){
    figS4.scat = figS4.scat+ xlab("log2 B73/Mo17")
  }
  
  legend = get_legend(figS4.scat)
  figS4.scat = figS4.scat+theme(legend.position = "none")
  ###Add density curves
  densityplotX=axis_canvas(figS4.scat, axis="x")+
    geom_density(aes(figS4[,"l2.ratio"], color=figS4[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    theme_void()+
    theme(legend.title = element_blank())+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))
  #set up marginal density plot Y axis
  densityplotY=axis_canvas(figS4.scat, axis="y", coord_flip = T)+
    geom_density(aes(figS4[,"l10.ttest"], color=figS4[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_linetype_manual(values=c("dotted", "solid"))+
    theme_void()+
    theme(legend.title = element_blank())+
    coord_flip(expand = F)+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))
  #put two plots together
  q = insert_xaxis_grob(figS4.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
  q = insert_yaxis_grob(q, densityplotY, grid::unit((0.275), "null"), position = "right")
  if(i==1){
    plotList=list()
    plotList[[1]]=q
  }
  else{
    plotList = c(plotList, list(q))
  }
}
plotList = c(plotList, list(legend))

#setup layout for panels
des = "A
       B
       C
       D"
r=patchwork::wrap_plots(A=plotList[[1]], B=plotList[[2]], C=plotList[[3]], D=plotList[[4]], heights = c(1,1,1,0.5), nrow = 4, design = des)

ggsave(plot=r, filename =paste0(PathToPlots, "FigureS4.pdf"), dpi=300, width=3.4252, height=8)