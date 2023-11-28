require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)
require(BSDA)

#combine seedling leaf (SL) and leaf blade (LB) protein data
fs2 = merge.data.frame(tmt.LB[,c("Accession", "l2.B73xMo17.MP", "l2.B73.Mo17")], tmt.SL[,c("Accession", "l2.B73xMo17.MP", "l2.B73.Mo17")], by.x = "Accession", by.y = "Accession", all=F)
#Set plot axis limits
fs2[,2:5][fs2[,2:5] <= -1] = -1
fs2[,2:5][fs2[,2:5] >= 1] = 1
#Set GOIs
fs2$Gene = substr(fs2$Accession, start = 1, stop = 14)
fs2$Category = "Other"
fs2[(fs2$Gene %in% GOIs$plastid),"Category"] = "Plastid-localized"
fs2[(fs2$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
fs2[(fs2$Gene %in% GOIs$PE.ribo),"Category"] = "PE plastid ribosome"
fs2[(fs2$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
fs2[(fs2$Gene %in% GOIs$PhAPGs),"Category"] = "PhAPGs"
#factor
fs2$Category = factor(fs2$Category, levels = c("Other", "Plastid-localized", "NE plastid ribosome","PE plastid ribosome", "PhANGs", "PhAPGs"))
fs2 = fs2[(order(fs2$Category)),]
#colorblind-friendly palette
colors.p <- c("Other"="#D8D8D8", "Plastid-localized"="black",
              "PE plastid ribosome"="#6A00EA","NE plastid ribosome"="#E5A3FF",
              "PhAPGs"="#008275", "PhANGs"="#81E401")
#Plot Fig S2 A
fs2a.scat = ggplot(fs2, mapping=aes(x=l2.B73xMo17.MP.x, y=l2.B73xMo17.MP.y, color=Category, shape=Category))+
  geom_point(aes(color=Category), size=1)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_shape_manual(values=c(1,19,19,19,19,19,19),guide="none")+
  theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  theme_classic()+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  geom_hline(yintercept=0, linetype="dashed", color = "gray50")+
  theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title = element_text(size=10),
        legend.position="bottom", legend.title = element_blank(),  plot.title = element_text(size=10, vjust = -3, hjust = -0.4), plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  ggtitle("A")+
  xlab("log2 BxM/MP LB")+
  ylab("log2 BxM/MP SL")+
  theme(aspect.ratio=1)+
  coord_fixed(expand = F)
fs2a.scat = fs2a.scat+theme(legend.position = "none")
##Add density curves
###set up marginal density plot X axis
densityplotX=axis_canvas(fs2a.scat, axis="x")+
  geom_density(aes(fs2[,"l2.B73xMo17.MP.x"], color=fs2[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())
###set up marginal density plot Y axis
densityplotY=axis_canvas(fs2a.scat, axis="y", coord_flip = T)+
  geom_density(aes(fs2[,"l2.B73xMo17.MP.y"], color=fs2[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())+
  coord_flip()
###put two plots together
A = insert_xaxis_grob(fs2a.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
A = insert_yaxis_grob(A, densityplotY, grid::unit(0.275, "null"), position = "right")
#Plot Fig S2 B
fs2b.scat = ggplot(fs2, mapping=aes(x=l2.B73.Mo17.x, y=l2.B73.Mo17.y, color=Category, shape=Category))+
  geom_point(aes(color=Category), size=1)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_shape_manual(values=c(1,19,19,19,19,19,19),guide="none")+
  theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  theme_classic()+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  geom_hline(yintercept=0, linetype="dashed", color = "gray50")+
  theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title = element_text(size=10),
        legend.position="bottom", legend.title = element_blank(),  plot.title = element_text(size=10, vjust = -3, hjust = -0.4), plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  guides(col=guide_legend(ncol=3))+
  ggtitle("B")+
  xlab("log2 B73/Mo17 LB")+
  ylab("log2 B73/Mo17 SL")+
  theme(aspect.ratio=1)+
  coord_fixed(expand = F)
legend = get_legend(fs2b.scat)
fs2b.scat = fs2b.scat+theme(legend.position = "none")
##Add density curves
###set up marginal density plot X axis
densityplotX=axis_canvas(fs2b.scat, axis="x")+
  geom_density(aes(fs2[,"l2.B73.Mo17.x"], color=fs2[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())
###set up marginal density plot Y axis
densityplotY=axis_canvas(fs2b.scat, axis="y", coord_flip = T)+
  geom_density(aes(fs2[,"l2.B73.Mo17.y"], color=fs2[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())+
  coord_flip()
###put two plots together
B = insert_xaxis_grob(fs2b.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
B = insert_yaxis_grob(B, densityplotY, grid::unit(0.275, "null"), position = "right")
#define panel layout
fs2.design = "AB
              CC"
fs2ab=patchwork::wrap_plots(A=A, B=B, C=legend,  nrow = 2, heights = c(3,1.2), design = fs2.design)

ggsave(plot=fs2ab, filename =paste0(PathToPlots, "FigureS2.pdf"), dpi=300, width=4.488189, height=3)