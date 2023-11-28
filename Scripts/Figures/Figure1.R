require(ggplot2)
require(cowplot)
require(patchwork)
require(dplyr)
require(reshape2)
require(BSDA)

#make list of dfs
tmt.SL.fig1.l = list("Hyb2HP"=tmt.SL[,c("Gene", "l2.B73xMo17.HP", "l10.ttest.HP")],
                    "Hyb2MP"=tmt.SL[,c("Gene", "l2.B73xMo17.MP", "l10.ttest.MP")],
                    "Par2Par"=tmt.SL[,c("Gene", "l2.B73.Mo17", "l10.ttest.P")])
#name plot panels
tmt.SL.fig1.l$Hyb2HP$plot = "A"
tmt.SL.fig1.l$Hyb2MP$plot = "B"
tmt.SL.fig1.l$Par2Par$plot = "C"
#rename column to match
tmt.SL.fig1.l = lapply(tmt.SL.fig1.l, setNames, c("Gene","l2.ratio", "l10.ttest","plot"))
#bind
tmt.SL.fig1 = do.call("rbind",tmt.SL.fig1.l)
#subset GOIs
GOIs.fig1 = GOIs[c("PhANGs","PhAPGs","NE.PtRibo","PE.ribo")]
names(GOIs.fig1) = c("PhANGs","PhAPGs","NEplastidRibosome","PEplastidRibosome")
#colorblind-friendly palette
colors.p <- c("Other"="#D8D8D8", "Plastid-localized"="black",
              "Cytosolic ribosome"="#B50000",
              "PE plastid ribosome"="#6A00EA","NE plastid ribosome"="#E5A3FF",
              "PhAPGs"="#008275", "PhANGs"="#81E401")
#Plot
plot=unique(tmt.SL.fig1[,"plot"])
plotList=vector(mode = "list", length = length(plot))
plotLab=unique(tmt.SL.fig1[,"plot"])
for(i in 1:length(plot)){
  f1=tmt.SL.fig1[which(tmt.SL.fig1[,"plot"]==plot[i]),]
  ###Set limits
  f1$l2.ratio[f1$l2.ratio <= -1] =-1
  f1$l2.ratio[f1$l2.ratio >= 1] =1
  ###Set GOIs
  f1$Category ="Other"
  f1[(f1$Gene %in% GOIs$plastid),"Category"] ="Plastid-localized"
  f1[(f1$Gene %in% GOIs$PhANGs),"Category"] ="PhANGs"
  f1[(f1$Gene %in% GOIs$PhAPGs),"Category"] ="PhAPGs"
  f1[(f1$Gene %in% GOIs$NE.PtRibo),"Category"] ="NE plastid ribosome"
  f1[(f1$Gene %in% GOIs$PE.ribo),"Category"] ="PE plastid ribosome"
  ###factor
  f1$Category =factor(f1$Category, levels = c("Other", "Plastid-localized", "PhANGs", "PhAPGs", "NE plastid ribosome", "PE plastid ribosome"))
  f1 =f1[(order(f1$Category)),]
  #scatter plot with marginal density curves
  f1.scat =ggplot(f1, mapping=aes(x=l2.ratio, y=l10.ttest, color=Category, shape=Category))+
    geom_point(aes(color=Category), size=1.5)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_shape_manual(values=c(1,19,19,19,19,19))+
    theme_classic()+
    xlim(c(-1,1))+
    geom_vline(xintercept=0, linetype="dashed", color = "gray40", size=0.5)+
    geom_hline(yintercept= 1.30103, linetype="dashed", color = "gray40", size=0.5)+
    ylab("-log10(p-value)")+
    ggtitle(plotLab[i])+
    theme(aspect.ratio=1)+
    coord_fixed(expand = F)+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))+
    theme(axis.text =element_text(size=10), legend.text = element_text(size=10), axis.title = element_text(size=10),
          legend.position="bottom", legend.title = element_blank(), plot.title = element_text(size=10, vjust = -3, hjust = -0.2), legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(-3, 'cm'))+
    guides(col=guide_legend(ncol = 2))
  if(i==1){
    f1.scat =f1.scat+ xlab("log2 BxM/HP")
  }
  else if(i==2){
    f1.scat =f1.scat+ xlab("log2 BxM/MP")
  }
  else if(i==3){
    f1.scat =f1.scat+ xlab("log2 B73/Mo17")
  }
  
  legend =get_legend(f1.scat)
  f1.scat =f1.scat+theme(legend.position = "none")
  ###Add density curves
  densityplotX=axis_canvas(f1.scat, axis="x")+
    geom_density(aes(f1[,"l2.ratio"], color=f1[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    theme_void()+
    theme(legend.title = element_blank())+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))
  #set up marginal density plot Y axis
  densityplotY=axis_canvas(f1.scat, axis="y", coord_flip = T)+
    geom_density(aes(f1[,"l10.ttest"], color=f1[,"Category"]), size=0.5)+
    scale_color_manual(values=colors.p, name="Category")+
    scale_linetype_manual(values=c("dotted", "solid"))+
    theme_void()+
    theme(legend.title = element_blank())+
    coord_flip(expand = F)+
    theme(plot.margin=grid::unit(c(-3,-3,0,-3), "mm"))
  #put two plots together
  q =insert_xaxis_grob(f1.scat, densityplotX, grid::unit(0.275, "null"), position = "top")
  q =insert_yaxis_grob(q, densityplotY, grid::unit((0.275), "null"), position = "right")
  if(i==1){
    plotList=list()
    plotList[[1]]=q
  }
  else{
    plotList =c(plotList, list(q))
  }
}
plotList = c(plotList, list(legend))
#setup layout of figure panels
des = "A
       B
       C
       D"
r=patchwork::wrap_plots(A=plotList[[1]], B=plotList[[2]], C=plotList[[3]], D=plotList[[4]], heights = c(1,1,1,0.5), nrow = 4, design = des)

ggsave(plot=r, filename =paste0(PathToPlots, "Figure1.pdf"), dpi=300, width=3.4252, height=8)  

#Calculate statistics
ztestfun = function(df, name){
  df.l = list()
  for(i in 1:length(GOIs.fig1)){
    #subset to GOI
    df.x = df$l2.ratio[df$Gene %in% GOIs.fig1[[i]]]
    #subset background
    df.y = df$l2.ratio[!df$Gene %in% GOIs.fig1[[i]]]
    #Calculate sigma.x
    sigx = sd(df.x)
    #Calculate sigma.y
    sigy = sd(df.y)
    #Calculate Z test
    z = z.test(x = df.x, y = df.y, sigma.x = sigx, sigma.y = sigy)
    #convert to dataframe
    z = z%>%
      unlist()%>%
      as.data.frame()
    colnames(z) = "value"
    z$stat = rownames(z)
    z$variable = paste0(name, "_", names(GOIs.fig1[i]))
    df.l[[i]] = z[1:9,]
  }
  names(df.l) = names(GOIs.fig1)
  df.l
}

#calculate z test for each figure panel
tmt.SL.Hyb2HP.ztest = ztestfun(tmt.SL.fig1.l$Hyb2HP, "Hyb/HP")
tmt.SL.Hyb2MP.ztest = ztestfun(tmt.SL.fig1.l$Hyb2MP, "Hyb/MP")
tmt.SL.Par2Par.ztest = ztestfun(tmt.SL.fig1.l$Par2Par, "B73/Mo17")

#combine statistics into single dataframe
tmt.SL.Hyb2HP.ztest = do.call("rbind", tmt.SL.Hyb2HP.ztest)
tmt.SL.Hyb2MP.ztest = do.call("rbind", tmt.SL.Hyb2MP.ztest)
tmt.SL.Par2Par.ztest = do.call("rbind", tmt.SL.Par2Par.ztest)

tmt.SL.Ztest <- do.call("rbind", list(tmt.SL.Hyb2HP.ztest, tmt.SL.Hyb2MP.ztest, tmt.SL.Par2Par.ztest))
tmt.SL.Ztest <- tmt.SL.Ztest[order(tmt.SL.Ztest$variable, tmt.SL.Ztest$stat),c("variable","stat","value")]
tmt.SL.Ztest$figure <- "Figure 1"
rownames(tmt.SL.Ztest) <- NULL