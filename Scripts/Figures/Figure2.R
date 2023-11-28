require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)

#select hyb/HP cols in order of plant height heterosis levels
tmt.6H.f2.hyb2HP = tmt.6H[,c("Gene","l2.A682xB73.HP","l2.B84xMo17.HP","l2.Mo17xB73.HP","l2.B73xMo17.HP","l2.B84xB73.HP","l2.A682xMo17.HP")]
#select par/par cols in same order
tmt.6H.f2.par2par = tmt.6H[,c("Gene","l2.A682.B73","l2.B84.Mo17","l2.Mo17.B73","l2.B73.Mo17","l2.B84.B73","l2.A682.Mo17")]
#rename cols to match figure panels
colnames(tmt.6H.f2.hyb2HP) = c("Gene", LETTERS[seq(from=1, to=6)])
colnames(tmt.6H.f2.par2par) = c("Gene", LETTERS[seq(from=1, to=6)])
#combine into single df
tmt.6H.f2.hyb2HP = melt(tmt.6H.f2.hyb2HP)
tmt.6H.f2.par2par = melt(tmt.6H.f2.par2par)
colnames(tmt.6H.f2.hyb2HP) = c("Gene","plot","l2.HYB.HP")
colnames(tmt.6H.f2.par2par) = c("Gene","plot","l2.parent.parent")
tmt.6H.f2 = tmt.6H.f2.hyb2HP %>% inner_join(tmt.6H.f2.par2par, by=c('Gene','plot'))
#colorblind-friendly palette
colors.p <- c("Other"="#D8D8D8", "Plastid-localized"="black",
              "Cytosolic ribosome"="#B50000",
              "PE plastid ribosome"="#6A00EA","NE plastid ribosome"="#E5A3FF",
              "PhAPGs"="#008275", "PhANGs"="#81E401")

plot=unique(tmt.6H.f2[,"plot"])
plotList=vector(mode = "list", length = length(plot))
plotLab=unique(tmt.6H.f2[,"plot"])
for(i in 1:length(plot)){
  f2=tmt.6H.f2[which(tmt.6H.f2[,"plot"]==plot[i]),]
  ###Set limits
  f2$l2.HYB.HP[f2$l2.HYB.HP <= -0.75] = -0.75
  f2$l2.HYB.HP[f2$l2.HYB.HP >= 0.75] = 0.75
  
  f2$l2.parent.parent[f2$l2.parent.parent <= -1] = -1
  f2$l2.parent.parent[f2$l2.parent.parent >= 1] = 1
  ###Set GOIs
  f2 = f2
  f2$Category = "Other"
  f2[(f2$Gene %in% GOIs$plastid),"Category"] = "Plastid-localized"
  f2[(f2$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
  f2[(f2$Gene %in% GOIs$PhAPGs),"Category"] = "PhAPGs"
  f2[(f2$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
  f2[(f2$Gene %in% GOIs$PE.ribo),"Category"] = "PE plastid ribosome"
  ###factor
  f2$Category = factor(f2$Category, levels = c("Other", "Plastid-localized", "PhANGs", "PhAPGs", "NE plastid ribosome", "PE plastid ribosome"))
  f2 = f2[(order(f2$Category)),]
  #scatter plot with marginal density curves
  f2.scat = ggplot(f2, mapping=aes(x=l2.HYB.HP, y=l2.parent.parent, color=Category, shape=Category))+
    geom_point(aes(color=Category), size=1)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_shape_manual(values=c(1,19,19,19,19,19), guide="none")+
    theme_classic()+
    xlim(c(-0.75,0.75))+
    ylim(c(-1,1))+
    geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
    geom_hline(yintercept=0, linetype="dashed", color = "gray60", size=0.5)+
    ggtitle(plotLab[i])+
    theme(aspect.ratio=1)+
    coord_fixed(expand = F)+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
    theme(axis.text=element_text(size=8), legend.text=element_text(size=8),axis.title = element_text(size = 8),
          legend.position="bottom", legend.title = element_blank(), plot.title = element_text(size=10, vjust = -3.5, hjust = -0.485), legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(-3, 'cm'))+
    guides(col=guide_legend(ncol = 2))
  legend = get_legend(f2.scat)
  f2.scat = f2.scat+theme(legend.position = "none")
  if(i==1){
    f2.scat = f2.scat+ xlab("log2 A682xB73/HP")+ ylab("log2 A682/B73")
  }
  else if(i==2){
    f2.scat = f2.scat+ xlab("log2 B84xMo17/HP")+ ylab("log2 B84/Mo17")
  }
  else if(i==3){
    f2.scat = f2.scat+ xlab("log2 Mo17xB73/HP")+ ylab("log2 Mo17/B73")
  }
  else if(i==4){
    f2.scat = f2.scat+ xlab("log2 B73xMo17/HP")+ ylab("log2 B73/Mo17")
  }
  else if(i==5){
    f2.scat = f2.scat+ xlab("log2 B84xB73/HP")+ ylab("log2 B84/B73")
  }
  else if(i==6){
    f2.scat = f2.scat+ xlab("log2 A682xMo17/HP")+ ylab("log2 A682/Mo17")
  }
  ###Add density curves
  densityplotX=axis_canvas(f2.scat, axis="x")+
    geom_density(aes(f2[,"l2.HYB.HP"], color=f2[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
    theme_void()+
    theme(legend.title = element_blank())+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))
  #set up marginal density plot Y axis
  densityplotY=axis_canvas(f2.scat, axis="y", coord_flip = T)+
    geom_density(aes(f2[,"l2.parent.parent"], color=f2[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_linetype_manual(values=c("dotted", "solid"))+
    geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
    theme_void()+
    theme(legend.title = element_blank())+
    coord_flip()+
    theme(plot.margin=grid::unit(c(-2,-2,0,-2), "mm"))
  #put two plots together
  q = insert_xaxis_grob(f2.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
  q = insert_yaxis_grob(q, densityplotY, grid::unit(0.275, "null"), position = "right")
  if(i==1){
    plotList=list()
    plotList[[1]]=q
  }
  else{
    plotList = c(plotList, list(q))
  }
}
plotList = c(plotList, list(legend))
#setup layout of figure panels
des = "AB
       CD
       EG
       HH"
r=patchwork::wrap_plots(A=plotList[[1]], B=plotList[[2]], C=plotList[[3]], D=plotList[[4]], E=plotList[[5]], G=plotList[[6]], H=plotList[[7]], nrow = 4, heights = c(1,1,1,0.5), design = des)

ggsave(plot=r, filename =paste0(PathToPlots, "Figure2.pdf"), dpi=300, width=3.4252, height=5.5)  