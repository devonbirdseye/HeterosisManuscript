require(ggplot2)
require(ggridges)
require(dplyr)
require(reshape2)
require(BSDA)

#combine RIL and 6 hybrids Protein
tmt.RIL.6H.Hyb2MP = merge.data.frame(tmt.RIL[c("Accession","l2.B73xMo17.MP","l2.M0014.14B.MP","l2.M0014.14M.MP","l2.M0016.16B.MP","l2.M0016.16M.MP","l2.M0021.21B.MP","l2.M0021.21M.MP","l2.M0317.317B.MP","l2.M0317.317M.MP")], tmt.6H[,c("Accession","l2.B73xMo17.MP","l2.Mo17xB73.MP","l2.B84xB73.MP","l2.B84xMo17.MP","l2.A682xB73.MP","l2.A682xB73.MP","l2.A682xMo17.MP")], by.x = "Accession", by.y = "Accession", all = F)
#convert Accession column to gene column
tmt.RIL.6H.Hyb2MP$Gene = substr(tmt.RIL.6H.Hyb2MP$Accession, start = 1, stop = 14)
#run correlation function
tmt.RIL.6H.Hyb2MP.cor = CorFun(tmt.RIL.6H.Hyb2MP, 2:16, c(PlantHeights.ril.MP, PlantHeights.6h.MP))
#Combine data for plotting
f3 = tmt.RIL.6H.Hyb2MP.cor[,c("Gene", "cor")]
#Set GOIs
f3$Category = "Other"
f3[(f3$Gene %in% GOIs$CytoRibo),"Category"] = "Cytosolic ribosome"
f3[(f3$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
f3[(f3$Gene %in% GOIs$PE.ribo),"Category"] = "PE plastid ribosome"
f3[(f3$Gene %in% GOIs$TPR),"Category"] = "TPR"
f3[(f3$Gene %in% GOIs$PTAC),"Category"] = "PTAC"
f3[(f3$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
f3[(f3$Gene %in% GOIs$PhAPGs),"Category"] = "PhAPGs"
f3[(f3$Gene %in% GOIs$ProteinBioSyn),"Category"] = "Protein biosynthesis"
f3[(f3$Gene %in% GOIs$CarbonFixation),"Category"] = "Carbon fixation"
f3[(f3$Gene %in% GOIs$Oxidoreductase),"Category"] = "Oxidoreductase"
f3[(f3$Gene %in% GOIs$Protease),"Category"] = "Protease"
f3[(f3$Gene %in% GOIs$BiosynSecMet),"Category"] = "Biosynthesis of secondary metabolites"
f3[(f3$Gene %in% GOIs$AlphaLinAcidMet),"Category"] = "Alpha-linoleic acid metabolism"

f3$Category = factor(f3$Category, levels = c("Other", "PE plastid ribosome", "NE plastid ribosome","Cytosolic ribosome", "TPR", "Protein biosynthesis", "PTAC", "PhANGs", "PhAPGs", "Carbon fixation","Oxidoreductase","Protease","Biosynthesis of secondary metabolites","Alpha-linoleic acid metabolism"))

#Plot
ggplot(f3, mapping=aes(y=Category, x=cor, fill = stat(x)))+
  geom_density_ridges_gradient(alpha=0.7) +
  scale_fill_viridis_c(option = "B")+
  theme_ridges() +
  theme(axis.text =element_text(size=8), axis.title = element_text(size=8), legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5, size=8))+
  coord_cartesian(expand = F)+
  xlab("Pearson correlation")+
  xlim(-1,1)
ggsave(filename =paste0(PathToPlots, "Figure3.pdf"), dpi=300, width=3.4252, height=4)

#Calculate statistics
GOIs.cor = GOIs[c("AlphaLinAcidMet","BiosynSecMet","Protease","Oxidoreductase","CarbonFixation","PhAPGs","PhANGs","PTAC","ProteinBioSyn","TPR","CytoRibo","NE.PtRibo","PE.ribo")]
GOIs.cor.names = c("Alpha-linoleic acid metabolism", "Biosynthesis of secondary metabolites", "Protease","Oxidoreductase","Carbon fixation","PhAPGs", "PhANGs","PTAC", "Protein biosynthesis", "TPR","Cytosolic ribosome", "NE plastid ribosome",  "PE plastid ribosome")
tmt.RIL.6H.Hyb2MP.cor.Ztest <- list()
for(i in 1:13){
  tmt.RIL.6H.Hyb2MP.cor.x <- sd(tmt.RIL.6H.Hyb2MP.cor$cor[tmt.RIL.6H.Hyb2MP.cor$Gene %in% GOIs.cor[[i]]])
  tmt.RIL.6H.Hyb2MP.cor.y <- sd(tmt.RIL.6H.Hyb2MP.cor$cor[!tmt.RIL.6H.Hyb2MP.cor$Gene %in% GOIs.cor[[i]]])
  tmt.RIL.6H.Hyb2MP.cor.Ztest[[i]] <- as.data.frame(unlist(z.test(x = tmt.RIL.6H.Hyb2MP.cor$cor[tmt.RIL.6H.Hyb2MP.cor$Gene %in% GOIs.cor[[i]]], y = tmt.RIL.6H.Hyb2MP.cor$cor[!tmt.RIL.6H.Hyb2MP.cor$Gene %in% GOIs.cor[[i]]], sigma.x = tmt.RIL.6H.Hyb2MP.cor.x, sigma.y = tmt.RIL.6H.Hyb2MP.cor.y)))
  colnames(tmt.RIL.6H.Hyb2MP.cor.Ztest[[i]]) <- GOIs.cor.names[i]
  tmt.RIL.6H.Hyb2MP.cor.Ztest[[i]]$stat <- rownames(tmt.RIL.6H.Hyb2MP.cor.Ztest[[i]])
}
tmt.RIL.6H.Hyb2MP.cor.Ztest <- Reduce(function(x, y, ...) merge(x, y, all = TRUE, by.x = "stat", by.y = "stat"),tmt.RIL.6H.Hyb2MP.cor.Ztest)
tmt.RIL.6H.Hyb2MP.cor.Ztest <- melt(tmt.RIL.6H.Hyb2MP.cor.Ztest, id.vars = "stat")
tmt.RIL.6H.Hyb2MP.cor.Ztest <- tmt.RIL.6H.Hyb2MP.cor.Ztest[,c("variable", "stat", "value")]
tmt.RIL.6H.Hyb2MP.cor.Ztest <- tmt.RIL.6H.Hyb2MP.cor.Ztest[order(tmt.RIL.6H.Hyb2MP.cor.Ztest$variable, tmt.RIL.6H.Hyb2MP.cor.Ztest$stat),]
tmt.RIL.6H.Hyb2MP.cor.Ztest$figure <- "Figure 3"