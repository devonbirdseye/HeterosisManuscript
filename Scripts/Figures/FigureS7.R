require(pathview)
require(plyr)

#Plot fig S7A
#merge with mdb to get entrez 
cpm.SL.entrez = join(cpm.SL[,c("Gene", "l2.B73xMo17.MP")], Entrez2)
#remove proteins with no entrez or duplicate entrez
cpm.SL.entrez.nodups = cpm.SL.entrez[!is.na(cpm.SL.entrez$Entrez),]
cpm.SL.entrez.nodups = cpm.SL.entrez.nodups[!cpm.SL.entrez.nodups$Entrez=="",]
cpm.SL.entrez.nodups = cpm.SL.entrez.nodups[!duplicated(cpm.SL.entrez.nodups$Entrez),]
#make new dataframe with just l2(BM/MP) values and Entrez as rownames
cpm.SL.entrez.nodups.trim = data.frame("l2.B73xMo17.MP"=cpm.SL.entrez.nodups$l2.B73xMo17.MP)
rownames(cpm.SL.entrez.nodups.trim) = cpm.SL.entrez.nodups$Entrez
#Map
cpm.SL.entrez.nodups.trim.pv = pathview(gene.data = cpm.SL.entrez.nodups.trim, pathway.id = "zma04712", species = "zma",  bins = 20, limit = 1, out.suffix="cpm.SL.BxM")

#Plot fig S7B
#merge with mdb to get entrez 
cpm.RIL.6H.Hyb2MP.cor.entrez = merge.data.frame(cpm.RIL.6H.Hyb2MP.cor, mgdb[,c(1,12)], by.x = "Gene", by.y = "v4_gene_model", all = F)
#remove proteins with no entrez or duplicate entrez
cpm.RIL.6H.Hyb2MP.cor.entrez.nodups = cpm.RIL.6H.Hyb2MP.cor.entrez[!is.na(cpm.RIL.6H.Hyb2MP.cor.entrez$Entrez),]
cpm.RIL.6H.Hyb2MP.cor.entrez.nodups = cpm.RIL.6H.Hyb2MP.cor.entrez.nodups[!cpm.RIL.6H.Hyb2MP.cor.entrez.nodups$Entrez=="",]
cpm.RIL.6H.Hyb2MP.cor.entrez.nodups = cpm.RIL.6H.Hyb2MP.cor.entrez.nodups[!duplicated(cpm.RIL.6H.Hyb2MP.cor.entrez.nodups$Entrez),]
#make new dataframe with just l2(BM/MP) values and Entrez as rownames
cpm.RIL.6H.Hyb2MP.cor.entrez.nodups.trim = data.frame("cor"=cpm.RIL.6H.Hyb2MP.cor.entrez.nodups$cor)
rownames(cpm.RIL.6H.Hyb2MP.cor.entrez.nodups.trim) = cpm.RIL.6H.Hyb2MP.cor.entrez.nodups$Entrez
#map
cpm.RIL.6H.Hyb2MP.cor.entrez.nodups.trim.pv = pathview(gene.data = cpm.RIL.6H.Hyb2MP.cor.entrez.nodups.trim, pathway.id = "zma04712", species = "zma",  bins = 20, limit = 1, out.suffix="cpm.cor", node.sum = "max.abs")