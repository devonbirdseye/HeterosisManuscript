#Combine SL and LB data
cpm.SL.LB = merge.data.frame(cpm.SL, cpm.LB, by.x = "Gene", by.y = "Gene", all=F)
#Calculate r-squared between each replicate
fs3v = c(
  summary(lm(B73.1~Mo17.1, data = cpm.SL))$adj.r.squared,
  summary(lm(B73xMo17.1~B73.1, data = cpm.SL))$adj.r.squared,
  summary(lm(B73xMo17.1~Mo17.1, data = cpm.SL))$adj.r.squared,
  summary(lm(B73xMo17.1~B73xMo17.2, data = cpm.SL))$adj.r.squared,
  summary(lm(B73.1~Mo17.1, data = cpm.LB))$adj.r.squared,
  summary(lm(B73xMo17.1~B73.1, data = cpm.LB))$adj.r.squared,
  summary(lm(B73xMo17.1~Mo17.1, data = cpm.LB))$adj.r.squared,
  summary(lm(B73xMo17.1~B73xMo17.2, data = cpm.LB))$adj.r.squared,
  summary(lm(B73.1.x~B73.1.y, data = cpm.SL.LB))$adj.r.squared,
  summary(lm(Mo17.1.x~Mo17.1.y, data = cpm.SL.LB))$adj.r.squared,
  summary(lm(B73xMo17.1.x~B73xMo17.1.y, data = cpm.SL.LB))$adj.r.squared
)
fs3df = data.frame("Rsquared"=fs3v, "comparison"=c("SL B73.1 vs SL Mo17.1", "SL BxM.1 vs SL B73.1", "SL BxM.1 vs SL Mo17.1", "SL.BxM.1 vs SL BxM.2", "LB B73.1 vs LB Mo17.1", "LB BxM.1 vs LB B73.1", "LB BxM.1 vs LB Mo17.1", "LB.BxM.1 vs LB BxM.2", "SL B73.1 vs LB B73.1", "SL Mo17.1 vs LB Mo17.1", "SL BxM.1 vs LB BxM.1"))
fs3df$comparison = factor(fs3df$comparison, levels = fs3df$comparison)
#Plot
ggplot(fs3df, mapping = aes(x=comparison, y=Rsquared))+
  geom_col()+
  theme_bw()+
  ylab("R squared")+
  xlab("")+
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = 8))
ggsave(paste0(PathToPlots,"FigureS3.pdf"), width = 7.007874, height = 6)