#Process Seedling Leaf (SL) RNA data (CPM)
cpm.SL <- as.data.frame(cpm.OG.data[,c("gid","BR193","BR194","BR195","BR196","BR197","BR198","BR199", "BR200","BR201")])
colnames(cpm.SL) <- c("Gene", "B73.1", "B73.2", "B73.3", "Mo17.1", "Mo17.2", "Mo17.3", "B73xMo17.1", "B73xMo17.2", "B73xMo17.3")
##Remove genes not detected in one or more rep
cpm.SL[cpm.SL==0] <- NA
cpm.SL <- na.omit(cpm.SL)
##Calculate averages of bio reps
cpm.SL$B73 <- rowMeans(cpm.SL[,c("B73.1", "B73.2", "B73.3")])
cpm.SL$Mo17 <- rowMeans(cpm.SL[,c("Mo17.1", "Mo17.2", "Mo17.3")])
cpm.SL$B73xMo17 <- rowMeans(cpm.SL[,c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3")])
##calculate ratios
cpm.SL$l2.B73xMo17.MP <- log2(cpm.SL$B73xMo17/rowMeans(cpm.SL[,c("B73", "Mo17")]))
cpm.SL$l2.B73.Mo17 <-log2(cpm.SL$B73/cpm.SL$Mo17)
##Calculate hyb/MP significance
cpm.SL$MP1 <- rowMeans(cpm.SL[,c("B73.1","Mo17.1")])
cpm.SL$MP2 <- rowMeans(cpm.SL[,c("B73.2","Mo17.2")])
cpm.SL$MP3 <- rowMeans(cpm.SL[,c("B73.3","Mo17.3")])
cpm.SL$ttest.MP <- ttestFun(cpm.SL, c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3"), c("MP1", "MP2", "MP3"))
cpm.SL$l10.ttest.MP <- -(log10(cpm.SL$ttest.MP))
##Define high-parent(HP)
cpm.SL$HP1 <- pmax(cpm.SL$B73.1, cpm.SL$Mo17.1)
cpm.SL$HP2 <- pmax(cpm.SL$B73.2, cpm.SL$Mo17.2)
cpm.SL$HP3 <- pmax(cpm.SL$B73.3, cpm.SL$Mo17.3)
##Define low-parent(LP)
cpm.SL$LP1 <- pmin(cpm.SL$B73.1, cpm.SL$Mo17.1)
cpm.SL$LP2 <- pmin(cpm.SL$B73.2, cpm.SL$Mo17.2)
cpm.SL$LP3 <- pmin(cpm.SL$B73.3, cpm.SL$Mo17.3)
##Calculate hyb/HP
cpm.SL$HP <- rowMeans(cpm.SL[,c("HP1", "HP2", "HP3")])
cpm.SL$l2.B73xMo17.HP <- log2(cpm.SL$B73xMo17/cpm.SL$HP)
##Calculate hyb/LP
cpm.SL$LP <- rowMeans(cpm.SL[,c("LP1", "LP2", "LP3")])
cpm.SL$l2.B73xMo17.LP <- log2(cpm.SL$B73xMo17/cpm.SL$LP)
##Calculate hyb/HP significance
cpm.SL$ttest.HP <- ttestFun(cpm.SL, c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3"), c("HP1", "HP2", "HP3"))
cpm.SL$l10.ttest.HP <- -(log10(cpm.SL$ttest.HP))
##Calculate hyb/LP significance
cpm.SL$ttest.LP <- ttestFun(cpm.SL, c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3"), c("LP1", "LP2", "LP3"))
cpm.SL$l10.ttest.LP <- -(log10(cpm.SL$ttest.LP))
##Calculate parent/parent significance
cpm.SL$ttest.P <- ttestFun(cpm.SL, c("B73.1", "B73.2", "B73.3"), c("Mo17.1", "Mo17.2", "Mo17.3"))
cpm.SL$l10.ttest.P <- -(log10(cpm.SL$ttest.P))
##subset to genes detected as proteins
cpm.SL.protein <- cpm.SL[cpm.SL$Gene%in%tmt.SL$Gene,]

#Process Leaf Blade (LB)
cpm.LB <- as.data.frame(cpm.OG.data[,c("gid","BR001", "BR002", "BR004", "BR003", "BR005", "BR007", "BR006", "BR008", "BR009")])
colnames(cpm.LB) <- c("Gene", "B73.1", "B73.2", "B73.3", "Mo17.1", "Mo17.2", "Mo17.3", "B73xMo17.1", "B73xMo17.2", "B73xMo17.3")
##Remove genes not detected in one or more rep
cpm.LB[cpm.LB==0] <- NA
cpm.LB <- na.omit(cpm.LB)
##Calculate averages of bio reps
cpm.LB$B73 <- rowMeans(cpm.LB[,c("B73.1", "B73.2", "B73.3")])
cpm.LB$Mo17 <- rowMeans(cpm.LB[,c("Mo17.1", "Mo17.2", "Mo17.3")])
cpm.LB$B73xMo17 <- rowMeans(cpm.LB[,c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3")])
##calculate ratios
cpm.LB$l2.B73xMo17.MP <- log2(cpm.LB$B73xMo17/rowMeans(cpm.LB[,c("B73", "Mo17")]))
cpm.LB$l2.B73.Mo17 <-log2(cpm.LB$B73/cpm.LB$Mo17)
##Calculate hyb/MP significance
cpm.LB$MP1 <- rowMeans(cpm.LB[,c("B73.1","Mo17.1")])
cpm.LB$MP2 <- rowMeans(cpm.LB[,c("B73.2","Mo17.2")])
cpm.LB$MP3 <- rowMeans(cpm.LB[,c("B73.3","Mo17.3")])
cpm.LB$ttest.MP <- ttestFun(cpm.LB, c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3"), c("MP1", "MP2", "MP3"))
cpm.LB$l10.ttest.MP <- -(log10(cpm.LB$ttest.MP))
##Define high-parent(HP)
cpm.LB$HP1 <- pmax(cpm.LB$B73.1, cpm.LB$Mo17.1)
cpm.LB$HP2 <- pmax(cpm.LB$B73.2, cpm.LB$Mo17.2)
cpm.LB$HP3 <- pmax(cpm.LB$B73.3, cpm.LB$Mo17.3)
##Calculate hyb/HP
cpm.LB$HP <- rowMeans(cpm.LB[,c("HP1", "HP2", "HP3")])
cpm.LB$l2.B73xMo17.HP <- log2(cpm.LB$B73xMo17/cpm.LB$HP)
##Calculate hyb/HP significance
cpm.LB$ttest.HP <- ttestFun(cpm.LB, c("B73xMo17.1", "B73xMo17.2", "B73xMo17.3"), c("HP1", "HP2", "HP3"))
cpm.LB$l10.ttest.HP <- -(log10(cpm.LB$ttest.HP))
##Calculate parent/parent significance
cpm.LB$ttest.P <- ttestFun(cpm.LB, c("B73.1", "B73.2", "B73.3"), c("Mo17.1", "Mo17.2", "Mo17.3"))
cpm.LB$l10.ttest.P <- -(log10(cpm.LB$ttest.P))
##subset to genes detected as protein
cpm.LB.protein <- cpm.LB[cpm.LB$Gene%in%tmt.LB$Gene,]

#Process 6 hybrid (6H) data
##Remove genes not detected in one or more rep
cpm.6H[cpm.6H==0] <- NA
cpm.6H <- na.omit(cpm.6H)
##Calculate averages
cpm.6H$B73 <- rowMeans(cpm.6H[,c("B73.2", "B73.3", "B73.4")])
cpm.6H$B84 <- rowMeans(cpm.6H[,c("B84.2", "B84.3", "B84.4")])
cpm.6H$Mo17 <- rowMeans(cpm.6H[,c("Mo17.2", "Mo17.3", "Mo17.4")])
cpm.6H$A682 <- rowMeans(cpm.6H[,c("A682.2", "A682.3", "A682.4")])
cpm.6H$B73xMo17 <- rowMeans(cpm.6H[,c("B73xMo17.2", "B73xMo17.3", "B73xMo17.4")])
cpm.6H$Mo17xB73 <- rowMeans(cpm.6H[,c("Mo17xB73.2", "Mo17xB73.3", "Mo17xB73.4")])
cpm.6H$B84xB73 <- rowMeans(cpm.6H[,c("B84xB73.2", "B84xB73.3", "B84xB73.4")])
cpm.6H$B84xMo17 <- rowMeans(cpm.6H[,c("B84xMo17.2", "B84xMo17.3", "B84xMo17.4")])
cpm.6H$A682xB73 <- rowMeans(cpm.6H[,c("A682xB73.2", "A682xB73.3", "A682xB73.4")])
cpm.6H$A682xMo17 <- rowMeans(cpm.6H[,c("A682xMo17.2", "A682xMo17.3", "A682xMo17.4")])
##calculate ratio of hybrid to midparent
cpm.6H$l2.B73xMo17.MP <- log2(cpm.6H$B73xMo17/rowMeans(cpm.6H[,c("B73","Mo17")]))
cpm.6H$l2.Mo17xB73.MP <- log2(cpm.6H$Mo17xB73/rowMeans(cpm.6H[,c("B73","Mo17")]))
cpm.6H$l2.B84xB73.MP <- log2(cpm.6H$B84xB73/rowMeans(cpm.6H[,c("B84","B73")]))
cpm.6H$l2.B84xMo17.MP <- log2(cpm.6H$B84xMo17/rowMeans(cpm.6H[,c("B84","Mo17")]))
cpm.6H$l2.A682xB73.MP <- log2(cpm.6H$A682xB73/rowMeans(cpm.6H[,c("A682","B73")]))
cpm.6H$l2.A682xMo17.MP <- log2(cpm.6H$A682xMo17/rowMeans(cpm.6H[,c("A682","Mo17")]))
##Calculate ratio of parent to parent
cpm.6H$l2.B73.Mo17 <- log2(cpm.6H$B73/cpm.6H$Mo17)
cpm.6H$l2.Mo17.B73 <- log2(cpm.6H$Mo17/cpm.6H$B73)
cpm.6H$l2.B84.B73 <- log2(cpm.6H$B84/cpm.6H$B73)
cpm.6H$l2.B84.Mo17 <- log2(cpm.6H$B84/cpm.6H$Mo17)
cpm.6H$l2.A682.B73 <- log2(cpm.6H$A682/cpm.6H$B73)
cpm.6H$l2.A682.Mo17 <- log2(cpm.6H$A682/cpm.6H$Mo17)
##Subset to genes with detected protein
cpm.6H.protein <- cpm.6H[cpm.6H$Gene%in%tmt.6H$Gene,]

#Process RIL data
##Remove genes not detected in one or mor rep
cpm.RIL[cpm.RIL==0] <- NA
cpm.RIL <- na.omit(cpm.RIL)
##Calculate averages
cpm.RIL$B73 <- rowMeans(cpm.RIL[,c("B73.1", "B73.2", "B73.3", "B73.4")])
cpm.RIL$Mo17 <- rowMeans(cpm.RIL[,c("Mo17.1", "Mo17.2", "Mo17.3", "Mo17.4")])
cpm.RIL$B73xMo17 <- rowMeans(cpm.RIL[,c("BxM.1", "BxM.2", "BxM.3", "BxM.4")])
cpm.RIL$M0014 <- rowMeans(cpm.RIL[,c("M0014.1", "M0014.2", "M0014.3", "M0014.4")])
cpm.RIL$M0014xB73 <- rowMeans(cpm.RIL[,c("M0014xB73.1", "M0014xB73.2", "M0014xB73.3", "M0014xB73.4")])
cpm.RIL$M0014xMo17 <- rowMeans(cpm.RIL[,c("M0014xMo17.1", "M0014xMo17.2", "M0014xMo17.3", "M0014xMo17.4")])
cpm.RIL$M0016 <- rowMeans(cpm.RIL[,c("M0016.1", "M0016.2", "M0016.3", "M0016.4")])
cpm.RIL$M0016xB73 <- rowMeans(cpm.RIL[,c("M0016xB73.1", "M0016xB73.2", "M0016xB73.3")])
cpm.RIL$M0016xMo17 <- rowMeans(cpm.RIL[,c("M0016xMo17.1", "M0016xMo17.2", "M0016xMo17.3", "M0016xMo17.4")])
cpm.RIL$M0021 <- rowMeans(cpm.RIL[,c("M0021.1", "M0021.2", "M0021.3", "M0021.4")])
cpm.RIL$M0021xB73 <- rowMeans(cpm.RIL[,c("M0021xB73.1", "M0021xB73.2", "M0021xB73.3", "M0021xB73.4")])
cpm.RIL$M0021xMo17 <- rowMeans(cpm.RIL[,c("M0021xMo17.1", "M0021xMo17.2", "M0021xMo17.3", "M0021xMo17.4")])
cpm.RIL$M0317 <- rowMeans(cpm.RIL[,c("M0317.1", "M0317.2", "M0317.3", "M0317.4")])
cpm.RIL$M0317xB73 <- rowMeans(cpm.RIL[,c("M0317xB73.1", "M0317xB73.2", "M0317xB73.3", "M0317xB73.4")])
cpm.RIL$M0317xMo17 <- rowMeans(cpm.RIL[,c("M0317xMo17.1", "M0317xMo17.2", "M0317xMo17.3", "M0317xMo17.4")])
##calculate ratio of hybrid to midparent
cpm.RIL$l2.B73xMo17.MP <- log2(cpm.RIL$B73xMo17/rowMeans(cpm.RIL[,c("B73","Mo17")]))
cpm.RIL$l2.M0014xB73.MP <- log2(cpm.RIL$M0014xB73/rowMeans(cpm.RIL[,c("M0014", "B73")]))
cpm.RIL$l2.M0014xMo17.MP <- log2(cpm.RIL$M0014xMo17/rowMeans(cpm.RIL[,c("M0014", "Mo17")]))
cpm.RIL$l2.M0016xB73.MP <- log2(cpm.RIL$M0016xB73/rowMeans(cpm.RIL[,c("M0016", "B73")]))
cpm.RIL$l2.M0016xMo17.MP <- log2(cpm.RIL$M0016xMo17/rowMeans(cpm.RIL[,c("M0016", "Mo17")]))
cpm.RIL$l2.M0021xB73.MP <- log2(cpm.RIL$M0021xB73/rowMeans(cpm.RIL[,c("M0021", "B73")]))
cpm.RIL$l2.M0021xMo17.MP <- log2(cpm.RIL$M0021xMo17/rowMeans(cpm.RIL[,c("M0021", "Mo17")]))
cpm.RIL$l2.M0317xB73.MP <- log2(cpm.RIL$M0317xB73/rowMeans(cpm.RIL[,c("M0317", "B73")]))
cpm.RIL$l2.M0317xMo17.MP <- log2(cpm.RIL$M0317xMo17/rowMeans(cpm.RIL[,c("M0317", "Mo17")]))
##Calculate ratio of parent to parent
cpm.RIL$l2.B73.Mo17 <- log2(cpm.RIL$B73/cpm.RIL$Mo17)
cpm.RIL$l2.M0014.B73 <- log2(cpm.RIL$M0014/cpm.RIL$B73)
cpm.RIL$l2.M0014.Mo17 <- log2(cpm.RIL$M0014/cpm.RIL$Mo17)
cpm.RIL$l2.M0016.B73 <- log2(cpm.RIL$M0016/cpm.RIL$B73)
cpm.RIL$l2.M0016.Mo17 <- log2(cpm.RIL$M0016/cpm.RIL$Mo17)
cpm.RIL$l2.M0021.B73 <- log2(cpm.RIL$M0021/cpm.RIL$B73)
cpm.RIL$l2.M0021.Mo17 <- log2(cpm.RIL$M0021/cpm.RIL$Mo17)
cpm.RIL$l2.M0317.B73 <- log2(cpm.RIL$M0317/cpm.RIL$B73)
cpm.RIL$l2.M0317.Mo17 <- log2(cpm.RIL$M0317/cpm.RIL$Mo17)
##Subset to genes with detected protein
cpm.RIL.protein <- cpm.RIL[cpm.RIL$Gene%in%tmt.RIL$Gene,]