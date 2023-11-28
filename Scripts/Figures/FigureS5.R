require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)

#Combine LB and SL data
fs5 = merge.data.frame(cpm.LB[,c("Gene", "l2.B73xMo17.MP", "l2.B73.Mo17")], cpm.SL[,c("Gene", "l2.B73xMo17.MP", "l2.B73.Mo17")], by.x = "Gene", by.y = "Gene", all=F)
###Set GOIs
fs5$Category = "Other"
fs5[(fs5$Gene %in% GOIs$plastid),"Category"] = "Plastid-localized"
fs5[(fs5$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
fs5[(fs5$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
###factor
fs5$Category = factor(fs5$Category, levels = c("Other", "Plastid-localized", "NE plastid ribosome","PhANGs"))
fs5 = fs5[(order(fs5$Category)),]
#Panel A
##Set plot axis limits
fs5a = fs5
fs5a[fs5a <= -2] = -2
fs5a[fs5a >= 2] = 2
colors.p <- c("Other"="#D8D8D8", "Plastid-localized"="black",
              "NE plastid ribosome"="#E5A3FF",
              "PhANGs"="#81E401")
##Scatterplot
fs5a.scat = ggplot(fs5a, mapping=aes(x=l2.B73xMo17.MP.x, y=l2.B73xMo17.MP.y, color=Category, shape=Category))+
  geom_point(aes(color=Category), size=1)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_shape_manual(values=c(1,19,19,19,19))+
  theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  theme_classic()+
  xlim(c(-1,1))+
  ylim(c(-1,1))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  geom_hline(yintercept=0, linetype="dashed", color = "gray50")+
  theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title = element_text(size=10),
        legend.position="none", legend.title = element_blank(),  plot.title = element_text(size=10, vjust = -3, hjust = -0.4), plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
  ggtitle("A")+
  xlab("log2 BxM/MP LB")+
  ylab("log2 BxM/MP SL")+
  theme(aspect.ratio=1)+
  coord_fixed(expand = F)
fs5a.scat = fs5a.scat+theme(legend.position = "none")
##Add density curves
###set up marginal density plot X axis
densityplotX=axis_canvas(fs5a.scat, axis="x")+
  geom_density(aes(fs5[,"l2.B73xMo17.MP.x"], color=fs5[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())
###set up marginal density plot Y axis
densityplotY=axis_canvas(fs5a.scat, axis="y", coord_flip = T)+
  geom_density(aes(fs5[,"l2.B73xMo17.MP.y"], color=fs5[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())+
  coord_flip()
###put two plots together
A = insert_xaxis_grob(fs5a.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
A = insert_yaxis_grob(A, densityplotY, grid::unit(0.275, "null"), position = "right")

#Panel B
##Set limits
fs5b = fs5
fs5b[fs5b <= -1] = -1
fs5b[fs5b >= 1] = 1
##Set GOIs
fs5b$Category = "Other"
fs5b[(fs5b$Gene %in% GOIs$plastid),"Category"] = "Plastid-localized"
fs5b[(fs5b$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
fs5b[(fs5b$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
##factor
fs5b$Category = factor(fs5b$Category, levels = c("Other", "Plastid-localized", "NE plastid ribosome", "PhANGs"))
fs5b = fs5b[(order(fs5b$Category)),]
##Scatterplot
fs5b.scat = ggplot(fs5b, mapping=aes(x=l2.B73.Mo17.x, y=l2.B73.Mo17.y, color=Category, shape=Category))+
  geom_point(aes(color=Category), size=1)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_shape_manual(values=c(1,19,19,19,19))+
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
legend = get_legend(fs5b.scat)
fs5b.scat = fs5b.scat+theme(legend.position = "none")
###Add density curves
####set up marginal density plot X axis
densityplotX=axis_canvas(fs5b.scat, axis="x")+
  geom_density(aes(fs5b[,"l2.B73.Mo17.x"], color=fs5b[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())
####set up marginal density plot Y axis
densityplotY=axis_canvas(fs5b.scat, axis="y", coord_flip = T)+
  geom_density(aes(fs5b[,"l2.B73.Mo17.y"], color=fs5b[,"Category"]), size=0.7)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray50")+
  theme_void()+
  theme(legend.title = element_blank())+
  coord_flip()
####put two plots together
B = insert_xaxis_grob(fs5b.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
B = insert_yaxis_grob(B, densityplotY, grid::unit(0.275, "null"), position = "right")

fs5.des = "AB
           CC"
fs5 =patchwork::wrap_plots(A=A, B=B, C=legend,  nrow = 2, heights = c(3,1.2), design = fs5.des)

ggsave(plot=fs5, filename =paste0(PathToPlots, "FigureS5.pdf"), dpi=300, width=4.488189, height=3)