require(plyr)
require(pathview)

#Retrieve entrez ID
##rename accession column to gene in conversion table
Entrez2 = Entrez
colnames(Entrez2)[1] = "Gene"
fs1 = join(tmt.SL[,c("Gene", "l2.B73xMo17.MP")], Entrez2)
#remove proteins with no entrez or duplicate entrez
fs1.nodups = fs1[!is.na(fs1$Entrez),]
fs1.nodups = fs1.nodups[!fs1.nodups$Entrez=="",]
fs1.nodups = fs1.nodups[!duplicated(fs1.nodups$Entrez),]
#make new dataframe with just l2(BM/MP) values and Entrez as rownames
fs1.nodups.trim = data.frame("l2.B73xMo17.MP"=fs1.nodups$l2.B73xMo17.MP)
rownames(fs1.nodups.trim) = fs1.nodups$Entrez
#Make map
pathview(gene.data = fs1.nodups.trim, pathway.id = "zma00710", species = "zma",  bins = 20, limit = 1, out.suffix="tmt.SL.BxM")