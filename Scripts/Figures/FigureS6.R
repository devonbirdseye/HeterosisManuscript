require(ggplot2)
require(dplyr)

#combine RIL and 6 hybrids Protein
fs6.df = merge.data.frame(select(cpm.RIL, c(Gene, contains(".MP"))), select(cpm.6H, c(Gene, contains(".MP"))), by.x = "Gene", by.y = "Gene", all = F)
#run correlation function
fs6.df = CorFun(fs6.df, 2:16, c(PlantHeights.ril.MP, PlantHeights.6h.MP))
#Combine data for plotting
fs6.df$data = "Hybrid/MP"
fs6.df = fs6.df[,c("Gene", "cor", "data")]
#Set GOIs
fs6.df$Category = "Other"
fs6.df[(fs6.df$Gene %in% GOIs$CytoRibo),"Category"] = "Cytosolic ribosome"
fs6.df[(fs6.df$Gene %in% GOIs$NE.PtRibo),"Category"] = "NE plastid ribosome"
fs6.df[(fs6.df$Gene %in% GOIs$TPR),"Category"] = "TPR"
fs6.df[(fs6.df$Gene %in% GOIs$PTAC),"Category"] = "PTAC"
fs6.df[(fs6.df$Gene %in% GOIs$PhANGs),"Category"] = "PhANGs"
fs6.df[(fs6.df$Gene %in% GOIs$ProteinBioSyn),"Category"] = "Protein biosynthesis"
fs6.df[(fs6.df$Gene %in% GOIs$CarbonFixation),"Category"] = "Carbon fixation"
fs6.df[(fs6.df$Gene %in% GOIs$Oxidoreductase),"Category"] = "Oxidoreductase"
fs6.df[(fs6.df$Gene %in% GOIs$Protease),"Category"] = "Protease"
fs6.df[(fs6.df$Gene %in% GOIs$BiosynSecMet),"Category"] = "Biosynthesis of secondary metabolites"
fs6.df[(fs6.df$Gene %in% GOIs$lphaLinAcidMet),"Category"] = "Alpha-linoleic acid metabolism"
#Factor
fs6.df$Category = factor(fs6.df$Category, levels = c("Other", "NE plastid ribosome","Cytosolic ribosome", "TPR","Protein biosynthesis", "PTAC", "PhANGs", "Carbon fixation","Oxidoreductase","Protease","Biosynthesis of secondary metabolites","Alpha-linoleic acid metabolism"))
#Plot
ggplot(fs6.df, mapping=aes(y=Category, x=cor, fill = stat(x)))+
  geom_density_ridges_gradient(alpha=0.7) +
  scale_fill_viridis_c(option = "B")+
  theme_ridges() +
  theme(axis.text =element_text(size=8), axis.title = element_text(size=8), legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_text(hjust = 0.5, size=10))+
  coord_cartesian(expand = F)+
  xlab("Pearson correlation")+
  xlim(-1,1)
ggsave(filename =paste0(PathToPlots, "FigureS6.pdf"), dpi=300, width=3.4252, height=4)