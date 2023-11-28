require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)
require(BSDA)

#define proteins significantly above mid-parent (AMP) and below mid-parent (BMP) in A682xB73 hybrid
tmt.6H.A682B.Sig.AMP = tmt.6H[tmt.6H$l2.A682xB73.MP > 0 & tmt.6H$ttest.A682xB73.MP < 0.05,]
tmt.6H.A682B.Sig.BMP = tmt.6H[tmt.6H$l2.A682xB73.MP < 0 & tmt.6H$ttest.A682xB73.MP < 0.05,]
#make column on acs data indicating whether the protein is AMP or BMP in A682xB73
f4a = tmt.acs[,c("Gene", "l2.ACS.B73","l10.pvalue")]
f4a$A682B = "MP"
f4a$A682B[f4a$Gene %in% tmt.6H.A682B.Sig.AMP$Gene] = "AMP"
f4a$A682B[f4a$Gene %in% tmt.6H.A682B.Sig.BMP$Gene] = "BMP"
f4a$A682B = factor(f4a$A682B, levels = c("MP","AMP", "BMP"))
f4a = f4a[order(f4a$A682B),]
#set plotting axis limits
f4a$l2.ACS.B73[f4a$l2.ACS.B73 < -2] = -2
f4a$l2.ACS.B73[f4a$l2.ACS.B73 > 2] = 2
f4a$l10.pvalue[f4a$l10.pvalue > 7.5] = 7.5
#merge acs with A682xB73 data
f4b = merge.data.frame(tmt.acs[,c("Accession", "l2.ACS.B73")], tmt.6H[,c("Accession", "l2.A682xB73.MP")], by.x = "Accession", by.y = "Accession", all=F)
#convert accession to gene
f4b$Gene = substr(f4b$Accession, start = 1, stop = 14)
#Set plotting axis limits
f4b$l2.A682xB73.MP[f4b$l2.A682xB73.MP <= -.75] = -.75
f4b$l2.A682xB73.MP[f4b$l2.A682xB73.MP >= .75] = .75
f4b$l2.ACS.B73[f4b$l2.ACS.B73 <= -.75] = -.75
f4b$l2.ACS.B73[f4b$l2.ACS.B73 >= .75] = .75
#Set GOIs
f4b = f4b
f4b$Category = "Other"
f4b[(f4b$Gene %in% GOIs$CytoRibo),"Category"] = "Cytosolic ribosome"
f4b[(f4b$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
f4b[(f4b$Gene %in% GOIs$PhAPGs),"Category"] = "PhAPGs"
f4b[(f4b$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
f4b[(f4b$Gene %in% GOIs$PE.ribo),"Category"] = "PE plastid ribosome"
f4b[(f4b$Gene %in% GOIs$EtBio),"Category"] = "Ethylene biosynthesis"
#factor
f4b$Category = factor(f4b$Category, levels = c("Other", "Cytosolic ribosome", "PhANGs", "PhAPGs", "NE plastid ribosome", "PE plastid ribosome", "Ethylene biosynthesis"))
f4b = f4b[(order(f4b$Category)),]

#Plot
#Plot A panel
f4a.vol = ggplot(f4a, mapping = aes(x=l2.ACS.B73, y=l10.pvalue, color=A682B))+
  geom_point(size=1)+
  scale_color_manual(values = c("MP"="grey", "AMP"="red", "BMP"="blue"), name="")+
  theme_bw()+
  xlab("log2 acs/B73")+
  ylab("-log10 p-value")+
  ggtitle("A")+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = -log10(0.05))+
  theme(aspect.ratio=1, legend.position = "bottom", axis.title = element_text(size=10), axis.text=element_text(size=10), legend.text = element_text(size=10), plot.title = element_text(hjust = -.05, vjust = -2))+
  coord_cartesian(expand = F)
f4a.vol.l = get_legend(f4a.vol)
f4a.vol = f4a.vol+theme(legend.position = "none")
#Plot B panel
#colorblind-friendly palette
colors.p <- c("Other"="#D8D8D8",
              "Cytosolic ribosome"="#B50000",
              "PE plastid ribosome"="#6A00EA","NE plastid ribosome"="#E5A3FF",
              "PhAPGs"="#008275", "PhANGs"="#81E401",
              "Ethylene biosynthesis"="#EF8F6D")
f4b.scat = ggplot(f4b, mapping=aes(x=l2.A682xB73.MP, y=l2.ACS.B73, color=Category, shape=Category))+
  geom_point(aes(color=Category), size=1)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_shape_manual(values=c(1,19,19,19,19,19,19))+
  theme_classic()+
  xlim(c(-.75,.75))+
  ylim(c(-.75,.75))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray60", size=0.5)+
  ylab("log2 acs/B73")+
  xlab("log2 A682xB73/MP")+
  ggtitle("B")+
  theme(aspect.ratio=1)+
  coord_fixed(expand = F)+
  theme(axis.text=element_text(size=10), legend.text=element_text(size=10),axis.title = element_text(size = 10),
        legend.position="bottom", legend.title = element_blank(), plot.title = element_text(hjust = -.55, vjust = -2))+
  guides(col=guide_legend(ncol = 2))
f4b.scat.l = get_legend(f4b.scat)
f4b.scat = f4b.scat+theme(legend.position = "none")

###Add density curves
densityplotX = axis_canvas(f4b.scat, axis="x")+
  geom_density(aes(f4b[,"l2.A682xB73.MP"], color=f4b[,"Category"]), size=0.5)+
  scale_color_manual(values=colors.p, name="Category")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
  theme_void()+
  theme(legend.title = element_blank())
#set up marginal density plot Y axis
densityplotY = axis_canvas(f4b.scat, axis="y", coord_flip = T)+
  geom_density(aes(f4b[,"l2.ACS.B73"], color=f4b[,"Category"]), size=0.5)+
  scale_color_manual(values=colors.p, name="Category")+
  scale_linetype_manual(values=c("dotted", "solid"))+
  geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
  theme_void()+
  theme(legend.title = element_blank())+
  coord_flip()
#put two plots together
q = insert_xaxis_grob(f4b.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
q = insert_yaxis_grob(q, densityplotY, grid::unit(0.275, "null"), position = "right")
#define panel layout
des = "A
       B
       C
       D"
r= patchwork::wrap_plots(A=f4a.vol, B=f4a.vol.l, C=q, D=f4b.scat.l,nrow = 4, heights = c(2.5,0.3,3,1.5), ncol = 1, design = des)

ggsave(plot=r, filename =paste0(PathToPlots, "Figure4.pdf"), dpi=300, width = 3.4252, height = 7)

#Calculate statistics
GOIs.acs = GOIs[c("CytoRibo","PhANGs","PhAPGs","NE.PtRibo","PE.ribo","EtBio")]
GOIs.acs.names = c("Cytosolic ribosome", "PhANGs", "PhAPGs", "NE plastid ribosome",  "PE plastid ribosome", "Ethylene biosynthesis")
tmt.acs.Ztest <- list()
for(i in 1:6){
  tmt.acs.x <- sd(tmt.acs$l2.ACS.B73[tmt.acs$Gene %in% GOIs.acs[[i]]])
  tmt.acs.y <- sd(tmt.acs$l2.ACS.B73[!tmt.acs$Gene %in% GOIs.acs[[i]]])
  tmt.acs.Ztest[[i]] <- as.data.frame(unlist(z.test(x = tmt.acs$l2.ACS.B73[tmt.acs$Gene %in% GOIs.acs[[i]]], y = tmt.acs$l2.ACS.B73[!tmt.acs$Gene %in% GOIs.acs[[i]]], sigma.x = tmt.acs.x, sigma.y = tmt.acs.y)))
  colnames(tmt.acs.Ztest[[i]]) <- GOIs.acs.names[i]
  tmt.acs.Ztest[[i]]$stat <- rownames(tmt.acs.Ztest[[i]])
}

tmt.acs.Ztest <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, by.x = "stat", by.y = "stat"),tmt.acs.Ztest)
tmt.acs.Ztest <- melt(tmt.acs.Ztest, id.vars = "stat")
tmt.acs.Ztest <- tmt.acs.Ztest[,c("variable", "stat", "value")]
tmt.acs.Ztest <- tmt.acs.Ztest[order(tmt.acs.Ztest$variable, tmt.acs.Ztest$stat),]
tmt.acs.Ztest$figure <- "Figure 4"

