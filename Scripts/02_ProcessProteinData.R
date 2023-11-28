#Subset Seedling Leaf (SL) Data
tmt.SL <- select(tmt,c(Group:Gene, SeedlingLeaf.B73.r1:SeedlingLeaf.BM.r3))
##Remove proteins not detected in one or more samples
tmt.SL <- na.omit(tmt.SL)
##Calculate averages of bio reps
tmt.SL$B73 <- rowMeans(select(tmt.SL,contains("B73")))
tmt.SL$Mo17 <- rowMeans(select(tmt.SL,contains("Mo17")))
tmt.SL$B73xMo17 <- rowMeans(select(tmt.SL,contains("BM")))
tmt.SL$MP <- rowMeans(tmt.SL[,c("B73", "Mo17")])
##calculate ratios
tmt.SL$l2.B73xMo17.MP <- log2(tmt.SL$B73xMo17/tmt.SL$MP)
tmt.SL$l2.B73.Mo17 <-log2(tmt.SL$B73/tmt.SL$Mo17)
##Calculate hyb/MP significance
tmt.SL$MP1 <- rowMeans(tmt.SL[,c("SeedlingLeaf.B73.r1","SeedlingLeaf.Mo17.r1")])
tmt.SL$MP2 <- rowMeans(tmt.SL[,c("SeedlingLeaf.B73.r2","SeedlingLeaf.Mo17.r2")])
tmt.SL$MP3 <- rowMeans(tmt.SL[,c("SeedlingLeaf.B73.r3","SeedlingLeaf.Mo17.r3")])
tmt.SL$ttest.MP <- ttestFun(tmt.SL, c("SeedlingLeaf.BM.r1", "SeedlingLeaf.BM.r2", "SeedlingLeaf.BM.r3"), c("MP1", "MP2", "MP3"))
tmt.SL$l10.ttest.MP <- -(log10(tmt.SL$ttest.MP))
##Define high-parent(HP)
tmt.SL$HP1 <- pmax(tmt.SL$SeedlingLeaf.B73.r1, tmt.SL$SeedlingLeaf.Mo17.r1)
tmt.SL$HP2 <- pmax(tmt.SL$SeedlingLeaf.B73.r2, tmt.SL$SeedlingLeaf.Mo17.r2)
tmt.SL$HP3 <- pmax(tmt.SL$SeedlingLeaf.B73.r3, tmt.SL$SeedlingLeaf.Mo17.r3)
##Calculate hyb/HP
tmt.SL$HP <- rowMeans(tmt.SL[,c("HP1", "HP2", "HP3")])
tmt.SL$l2.B73xMo17.HP <- log2(tmt.SL$B73xMo17/tmt.SL$HP)
##Calculate hyb/HP significance
tmt.SL$ttest.HP <- ttestFun(tmt.SL, c("SeedlingLeaf.BM.r1", "SeedlingLeaf.BM.r2", "SeedlingLeaf.BM.r3"), c("HP1", "HP2", "HP3"))
tmt.SL$l10.ttest.HP <- -(log10(tmt.SL$ttest.HP))
##Calculate Par/Par significance
tmt.SL$ttest.P <- ttestFun(tmt.SL, c("SeedlingLeaf.B73.r1", "SeedlingLeaf.B73.r2", "SeedlingLeaf.B73.r3"), c("SeedlingLeaf.Mo17.r1", "SeedlingLeaf.Mo17.r2", "SeedlingLeaf.Mo17.r3"))
tmt.SL$l10.ttest.P <- -(log10(tmt.SL$ttest.P))



#Subset Leaf Blade (LB) data
tmt.LB <- select(tmt,c(Group:Gene, MatureLeaf.B73.r1:MatureLeaf.BM.r3))
##Remove proteins not detected in one or more samples
tmt.LB <- na.omit(tmt.LB)
##Calculate averages of bio reps
tmt.LB$B73 <- rowMeans(select(tmt.LB,contains("B73")))
tmt.LB$Mo17 <- rowMeans(select(tmt.LB,contains("Mo17")))
tmt.LB$B73xMo17 <- rowMeans(select(tmt.LB,contains("BM")))
tmt.LB$MP <- rowMeans(tmt.LB[,c("B73", "Mo17")])
##calculate ratios
tmt.LB$l2.B73xMo17.MP <- log2(tmt.LB$B73xMo17/tmt.LB$MP)
tmt.LB$l2.B73.Mo17 <-log2(tmt.LB$B73/tmt.LB$Mo17)
##Calculate hyb/MP significance
tmt.LB$MP1 <- rowMeans(tmt.LB[,c("MatureLeaf.B73.r1","MatureLeaf.Mo17.r1")])
tmt.LB$MP2 <- rowMeans(tmt.LB[,c("MatureLeaf.B73.r2","MatureLeaf.Mo17.r2")])
tmt.LB$MP3 <- rowMeans(tmt.LB[,c("MatureLeaf.B73.r3","MatureLeaf.Mo17.r3")])
tmt.LB$ttest.MP <- ttestFun(tmt.LB, c("MatureLeaf.BM.r1", "MatureLeaf.BM.r2", "MatureLeaf.BM.r3"), c("MP1", "MP2", "MP3"))
tmt.LB$l10.ttest.MP <- -(log10(tmt.LB$ttest.MP))
##Define high-parent(HP)
tmt.LB$HP1 <- pmax(tmt.LB$MatureLeaf.B73.r1, tmt.LB$MatureLeaf.Mo17.r1)
tmt.LB$HP2 <- pmax(tmt.LB$MatureLeaf.B73.r2, tmt.LB$MatureLeaf.Mo17.r2)
tmt.LB$HP3 <- pmax(tmt.LB$MatureLeaf.B73.r3, tmt.LB$MatureLeaf.Mo17.r3)
##Calculate hyb/HP
tmt.LB$HP <- rowMeans(tmt.LB[,c("HP1", "HP2", "HP3")])
tmt.LB$l2.B73xMo17.HP <- log2(tmt.LB$B73xMo17/tmt.LB$HP)
##Calculate hyb/HP significance
tmt.LB$ttest.HP <- ttestFun(tmt.LB, c("MatureLeaf.BM.r1", "MatureLeaf.BM.r2", "MatureLeaf.BM.r3"), c("HP1", "HP2", "HP3"))
tmt.LB$l10.ttest.HP <- -(log10(tmt.LB$ttest.HP))
##Calculate Par/Par significance
tmt.LB$ttest.P <- ttestFun(tmt.LB, c("MatureLeaf.B73.r1", "MatureLeaf.B73.r2", "MatureLeaf.B73.r3"), c("MatureLeaf.Mo17.r1", "MatureLeaf.Mo17.r2", "MatureLeaf.Mo17.r3"))
tmt.LB$l10.ttest.P <- -(log10(tmt.LB$ttest.P))


#Subset 6-hybrid experiment (6H)
tmt.6H <- select(tmt,c(Group:Gene, contains("X6Hybrid")))
##Remove proteins not detected in one or more samples
tmt.6H <- na.omit(tmt.6H)
##calculate averages
tmt.6H$B73 <- rowMeans(select(tmt.6H,contains(".B73.")))
tmt.6H$Mo17 <- rowMeans(select(tmt.6H,contains(".Mo17.")))
tmt.6H$B73xMo17 <- rowMeans(select(tmt.6H,contains(".B73xMo17.")))
tmt.6H$Mo17xB73 <- rowMeans(select(tmt.6H,contains(".Mo17xB73.")))
tmt.6H$B84 <- rowMeans(select(tmt.6H,contains(".B84.")))
tmt.6H$B84xB73 <- rowMeans(select(tmt.6H,contains(".B84xB73.")))
tmt.6H$B84xMo17 <- rowMeans(select(tmt.6H,contains(".B84xMo17.")))
tmt.6H$A682 <- rowMeans(select(tmt.6H,contains(".A682.")))
tmt.6H$A682xB73 <- rowMeans(select(tmt.6H,contains(".A682xB73.")))
tmt.6H$A682xMo17 <- rowMeans(select(tmt.6H,contains(".A682xMo17.")))
##Calculate hyb/MP ratios
tmt.6H$l2.B73xMo17.MP <- log2(tmt.6H$B73xMo17/rowMeans(tmt.6H[,c("B73", "Mo17")]))
tmt.6H$l2.Mo17xB73.MP <- log2(tmt.6H$Mo17xB73/rowMeans(tmt.6H[,c("B73", "Mo17")]))
tmt.6H$l2.B84xB73.MP <- log2(tmt.6H$B84xB73/rowMeans(tmt.6H[,c("B84", "B73")]))
tmt.6H$l2.B84xMo17.MP <- log2(tmt.6H$B84xMo17/rowMeans(tmt.6H[,c("B84", "Mo17")]))
tmt.6H$l2.A682xB73.MP <- log2(tmt.6H$A682xB73/rowMeans(tmt.6H[,c("A682", "B73")]))
tmt.6H$l2.A682xMo17.MP <- log2(tmt.6H$A682xMo17/rowMeans(tmt.6H[,c("A682", "Mo17")]))
##Calculate parent/parent ratios
tmt.6H$l2.B73.Mo17 <- log2(tmt.6H$B73/tmt.6H$Mo17)
tmt.6H$l2.Mo17.B73 <- log2(tmt.6H$Mo17/tmt.6H$B73)
tmt.6H$l2.B84.B73 <- log2(tmt.6H$B84/tmt.6H$B73)
tmt.6H$l2.B84.Mo17 <- log2(tmt.6H$B84/tmt.6H$Mo17)
tmt.6H$l2.A682.B73 <- log2(tmt.6H$A682/tmt.6H$B73)
tmt.6H$l2.A682.Mo17 <- log2(tmt.6H$A682/tmt.6H$Mo17)
##Calculate MP reps
tmt.6H$MP1.B73Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r1","X6Hybrid.Mo17.r1")])
tmt.6H$MP2.B73Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r2","X6Hybrid.Mo17.r2")])
tmt.6H$MP3.B73Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r3","X6Hybrid.Mo17.r3")])
tmt.6H$MP1.B84B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r1","X6Hybrid.B84.r1")])
tmt.6H$MP2.B84B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r2","X6Hybrid.B84.r2")])
tmt.6H$MP3.B84B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r3","X6Hybrid.B84.r3")])
tmt.6H$MP1.B84Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r1","X6Hybrid.B84.r1")])
tmt.6H$MP2.B84Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r2","X6Hybrid.B84.r2")])
tmt.6H$MP3.B84Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r3","X6Hybrid.B84.r3")])
tmt.6H$MP1.A682B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r1","X6Hybrid.A682.r1")])
tmt.6H$MP2.A682B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r2","X6Hybrid.A682.r2")])
tmt.6H$MP3.A682B73 <- rowMeans(tmt.6H[,c("X6Hybrid.B73.r3","X6Hybrid.A682.r3")])
tmt.6H$MP1.A682Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r1","X6Hybrid.A682.r1")])
tmt.6H$MP2.A682Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r2","X6Hybrid.A682.r2")])
tmt.6H$MP3.A682Mo17 <- rowMeans(tmt.6H[,c("X6Hybrid.Mo17.r3","X6Hybrid.A682.r3")])
##Calculate pval
tmt.6H$ttest.B73xMo17.MP <- (ttestFun(tmt.6H, c("X6Hybrid.B73xMo17.r1", "X6Hybrid.B73xMo17.r2", "X6Hybrid.B73xMo17.r3"), c("MP1.B73Mo17", "MP2.B73Mo17", "MP3.B73Mo17")))
tmt.6H$ttest.Mo17xB73.MP <- (ttestFun(tmt.6H, c("X6Hybrid.Mo17xB73.r1", "X6Hybrid.Mo17xB73.r2", "X6Hybrid.Mo17xB73.r3"), c("MP1.B73Mo17", "MP2.B73Mo17", "MP3.B73Mo17")))
tmt.6H$ttest.B84xB73.MP <- (ttestFun(tmt.6H, c("X6Hybrid.B84xB73.r1", "X6Hybrid.B84xB73.r2", "X6Hybrid.B84xB73.r3"), c("MP1.B84B73", "MP2.B84B73", "MP3.B84B73")))
tmt.6H$ttest.B84xMo17.MP <- (ttestFun(tmt.6H, c("X6Hybrid.B84xMo17.r1", "X6Hybrid.B84xMo17.r2", "X6Hybrid.B84xMo17.r3"), c("MP1.B84Mo17", "MP2.B84Mo17", "MP3.B84Mo17")))
tmt.6H$ttest.A682xB73.MP <- (ttestFun(tmt.6H, c("X6Hybrid.A682xB73.r1", "X6Hybrid.A682xB73.r2", "X6Hybrid.A682xB73.r3"), c("MP1.A682B73", "MP2.A682B73", "MP3.A682B73")))
tmt.6H$ttest.A682xMo17.MP <- (ttestFun(tmt.6H, c("X6Hybrid.A682xMo17.r1", "X6Hybrid.A682xMo17.r2", "X6Hybrid.A682xMo17.r3"), c("MP1.A682Mo17", "MP2.A682Mo17", "MP3.A682Mo17")))
##Define high-parent(HP)
tmt.6H$HP.B73Mo17 <- pmax(tmt.6H$B73, tmt.6H$Mo17)
tmt.6H$HP.B84B73 <- pmax(tmt.6H$B73, tmt.6H$B84)
tmt.6H$HP.B84Mo17 <- pmax(tmt.6H$Mo17, tmt.6H$B84)
tmt.6H$HP.A682B73 <- pmax(tmt.6H$B73, tmt.6H$A682)
tmt.6H$HP.A682Mo17 <- pmax(tmt.6H$Mo17, tmt.6H$A682)
##Calculate Hyb/HP
tmt.6H$l2.B73xMo17.HP <- log2(tmt.6H$B73xMo17/tmt.6H$HP.B73Mo17)
tmt.6H$l2.Mo17xB73.HP <- log2(tmt.6H$Mo17xB73/tmt.6H$HP.B73Mo17)
tmt.6H$l2.B84xB73.HP <- log2(tmt.6H$B84xB73/tmt.6H$HP.B84B73)
tmt.6H$l2.B84xMo17.HP <- log2(tmt.6H$B84xMo17/tmt.6H$HP.B84Mo17)
tmt.6H$l2.A682xB73.HP <- log2(tmt.6H$A682xB73/tmt.6H$HP.A682B73)
tmt.6H$l2.A682xMo17.HP <- log2(tmt.6H$A682xMo17/tmt.6H$HP.A682Mo17)

#Subset RIL experiment (RIL)
tmt.RIL <- select(tmt,c(Group:Gene, B73.1:X317M.4))
##Remove proteins not detected in one or more samples
tmt.RIL <- na.omit(tmt.RIL)
##Calculate averages
tmt.RIL$BM.B73 <- rowMeans(tmt.RIL[,c("B73.1", "B73.3", "B73.4")])
tmt.RIL$BM.Mo17 <- rowMeans(tmt.RIL[,c("Mo17.1", "Mo17.3", "Mo17.4")])
tmt.RIL$BM.BM <- rowMeans(tmt.RIL[,c("BM.1", "BM.2", "BM.3", "BM.4")])
tmt.RIL$M0014.B73 <- rowMeans(tmt.RIL[,c("B73.1.1", "B73.2", "B73.3.1", "B73.4.1")])
tmt.RIL$M0014.Mo17<- rowMeans(tmt.RIL[,c("Mo17.1.1", "Mo17.2", "Mo17.3.1", "Mo17.4.1")])
tmt.RIL$M0014.14 <- rowMeans(tmt.RIL[,c("X14.1", "X14.2", "X14.3", "X14.4")])
tmt.RIL$M0014.14B <- rowMeans(tmt.RIL[,c("X14B.1", "X14B.2", "X14B.3", "X14B.4")])
tmt.RIL$M0014.14M <- rowMeans(tmt.RIL[,c("X14M.1", "X14M.2", "X14M.3", "X14M.4")])
tmt.RIL$M0016.B73 <- rowMeans(tmt.RIL[,c("B73.1.2", "B73.2.1", "B73.3.2", "B73.4.2")])
tmt.RIL$M0016.Mo17<- rowMeans(tmt.RIL[,c("Mo17.1.2", "Mo17.2.1", "Mo17.3.2", "Mo17.4.2")])
tmt.RIL$M0016.16 <- rowMeans(tmt.RIL[,c("X16.1", "X16.2", "X16.3", "X16.4")])
tmt.RIL$M0016.16B <- rowMeans(tmt.RIL[,c("X16B.1", "X16B.2", "X16B.3", "X16B.4")])
tmt.RIL$M0016.16M <- rowMeans(tmt.RIL[,c("X16M.1", "X16M.2", "X16M.3", "X16M.4")])
tmt.RIL$M0021.B73 <- rowMeans(tmt.RIL[,c("B73.1.3", "B73.2.2", "B73.3.3", "B73.4.3")])
tmt.RIL$M0021.Mo17<- rowMeans(tmt.RIL[,c("Mo17.1.3", "Mo17.2.2", "Mo17.3.3", "Mo17.4.3")])
tmt.RIL$M0021.21 <- rowMeans(tmt.RIL[,c("X21.1", "X21.2", "X21.3", "X21.4")])
tmt.RIL$M0021.21B <- rowMeans(tmt.RIL[,c("X21B.1", "X21B.2", "X21B.3", "X21B.4")])
tmt.RIL$M0021.21M <- rowMeans(tmt.RIL[,c("X21M.1", "X21M.2", "X21M.3", "X21M.4")])
tmt.RIL$M0317.B73 <- rowMeans(tmt.RIL[,c("B73.1.4", "B73.2.3", "B73.3.4", "B73.4.4")])
tmt.RIL$M0317.Mo17<- rowMeans(tmt.RIL[,c("Mo17.1.4", "Mo17.2.3", "Mo17.3.4", "Mo17.4.4")])
tmt.RIL$M0317.317 <- rowMeans(tmt.RIL[,c("X317.1", "X317.2", "X317.3", "X317.4")])
tmt.RIL$M0317.317B <- rowMeans(tmt.RIL[,c("X317B.1", "X317B.2", "X317B.3", "X317B.4")])
tmt.RIL$M0317.317M <- rowMeans(tmt.RIL[,c("X317M.1", "X317M.2", "X317M.3", "X317M.4")])
##Calculate hyb/mp ratios
tmt.RIL$l2.B73xMo17.MP <- log2(tmt.RIL$BM.BM/rowMeans(tmt.RIL[,c("BM.B73", "BM.Mo17")]))
tmt.RIL$l2.M0014.14B.MP <- log2(tmt.RIL$M0014.14B/rowMeans(tmt.RIL[,c("M0014.14", "M0014.B73")]))
tmt.RIL$l2.M0014.14M.MP <- log2(tmt.RIL$M0014.14M/rowMeans(tmt.RIL[,c("M0014.14", "M0014.Mo17")]))
tmt.RIL$l2.M0016.16B.MP <- log2(tmt.RIL$M0016.16B/rowMeans(tmt.RIL[,c("M0016.16", "M0016.B73")]))
tmt.RIL$l2.M0016.16M.MP <- log2(tmt.RIL$M0016.16M/rowMeans(tmt.RIL[,c("M0016.16", "M0016.Mo17")]))
tmt.RIL$l2.M0021.21B.MP <- log2(tmt.RIL$M0021.21B/rowMeans(tmt.RIL[,c("M0021.21", "M0021.B73")]))
tmt.RIL$l2.M0021.21M.MP <- log2(tmt.RIL$M0021.21M/rowMeans(tmt.RIL[,c("M0021.21", "M0021.Mo17")]))
tmt.RIL$l2.M0317.317B.MP <- log2(tmt.RIL$M0317.317B/rowMeans(tmt.RIL[,c("M0317.317", "M0317.B73")]))
tmt.RIL$l2.M0317.317M.MP <- log2(tmt.RIL$M0317.317M/rowMeans(tmt.RIL[,c("M0317.317", "M0317.Mo17")]))
##Calculate parent/parent ratios
tmt.RIL$l2.B73.Mo17 <- log2(tmt.RIL$BM.B73/tmt.RIL$BM.Mo17)
tmt.RIL$l2.M0014.14.B73 <- log2(tmt.RIL$M0014.14/tmt.RIL$M0014.B73)
tmt.RIL$l2.M0014.14.Mo17  <- log2(tmt.RIL$M0014.14/tmt.RIL$M0014.Mo17)
tmt.RIL$l2.M0016.16.B73 <- log2(tmt.RIL$M0016.16/tmt.RIL$M0016.B73)
tmt.RIL$l2.M0016.16.Mo17  <- log2(tmt.RIL$M0016.16/tmt.RIL$M0016.Mo17)
tmt.RIL$l2.M0021.21.B73 <- log2(tmt.RIL$M0021.21/tmt.RIL$M0021.B73)
tmt.RIL$l2.M0021.21.Mo17  <- log2(tmt.RIL$M0021.21/tmt.RIL$M0021.Mo17)
tmt.RIL$l2.M0317.317.B73 <- log2(tmt.RIL$M0317.317/tmt.RIL$M0317.B73)
tmt.RIL$l2.M0317.317.Mo17  <- log2(tmt.RIL$M0317.317/tmt.RIL$M0317.Mo17)

#Subset acs mutant experiment (acs)
tmt.acs <- select(tmt,c(Group:Gene, contains("acs")))
#@Remove proteins not detected in one or more samples
tmt.acs <- na.omit(tmt.acs)
#@Calculate averages of bio reps - acs
tmt.acs$B73avg <- rowMeans(select(tmt.acs,contains("B73.r")))
tmt.acs$ACSavg <- rowMeans(select(tmt.acs,contains("acs.r")))
#@calculate log2 ratio
tmt.acs$l2.ACS.B73 <- log2(tmt.acs$ACSavg/tmt.acs$B73avg)
#@calculate p-value
tmt.acs$ttest <- ttestFun(tmt.acs, 14:18, 19:23)
tmt.acs$l10.pvalue <- -log10(tmt.acs$ttest)