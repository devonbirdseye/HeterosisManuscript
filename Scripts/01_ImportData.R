#Import protein data
tmt <- read.csv(paste0(PathToData,"SupplementalTable9.csv"), stringsAsFactors = F)
##re-define protein column to shortened name
tmt <- subset(tmt, select=-Protein)
colnames(tmt)[colnames(tmt) == 'Gene'] <- 'Protein'
##make new column for gene accession without transcript id
tmt$Gene <- substr(tmt$Accession, start = 1, stop = 14)
##Remove polymorphic peptides
tmt.withPolymorphic <- tmt
tmt <- tmt[-(grep("only", tmt$Accession)),]
##Convert to numeric
tmt[1:170] <- sapply(tmt[1:170],as.numeric)

#Import RNA data
##Original experiment
cpm.OG.data <- read.table(paste0(PathToData, "cpm.tsv"), sep = '\t', header = TRUE)
##6-hybrid experiment
cpm.6H <- read.csv(paste0(PathToData, "SupplementalTable7.csv"), stringsAsFactors = F)
##RIL experiment
cpm.RIL <- read.csv(paste0(PathToData, "SupplementalTable8.csv"), stringsAsFactors = F)

#Import Plant Height data
PlantHeights <- read.csv(paste0(PathToData, "SupplementalTable4.csv"), stringsAsFactors = F)
#Make vector of 6H hyb/MP ratios
PlantHeights.6h.MP <- PlantHeights$Hyb.MP[PlantHeights$Experiment=="6H"]
#Make vector of RIL hyb/MP ratios
PlantHeights.ril.MP <- PlantHeights$Hyb.MP[PlantHeights$Experiment=="RIL"]

#Import Gene Of Interest (GOI) data
GOIs <- readRDS(paste0(PathToData, "GOIs"))
#Import mgdb accession conversion table
##Downloaded from MaizeGDB 2020/05/18, 07:11pm
mgdb <- read.csv(paste0(PathToData,"mgdb.csv"), stringsAsFactors = F, header = F)
colnames(mgdb) <- mgdb[5,]
mgdb <- mgdb[-c(1:5),]
#Import Entrez conversion table
Entrez <- readRDS(paste0(PathToData,"Entrez"))