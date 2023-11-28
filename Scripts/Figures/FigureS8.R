require(pathview)
require(plyr)

#Panel A
#merge with mgdb to get entrez 
tmt.A682B.entrez = merge.data.frame(tmt.6H, mgdb[,c(1,12)], by.x = "Gene", by.y = "v4_gene_model", all = F)
#remove proteins with no entrez or duplicate entrez
tmt.A682B.entrez.nodups = tmt.A682B.entrez[!is.na(tmt.A682B.entrez$Entrez),]
tmt.A682B.entrez.nodups = tmt.A682B.entrez.nodups[!tmt.A682B.entrez.nodups$Entrez=="",]
tmt.A682B.entrez.nodups = tmt.A682B.entrez.nodups[!duplicated(tmt.A682B.entrez.nodups$Entrez),]
#make new dataframe with just l2(BM/MP) values and Entrez as rownames
tmt.A682B.entrez.nodups.trim = data.frame("l2.A682xB73.MP"=tmt.A682B.entrez.nodups$l2.A682xB73.MP)
rownames(tmt.A682B.entrez.nodups.trim) = tmt.A682B.entrez.nodups$Entrez
#map
tmt.A682B.entrez.nodups.trim.pv = pathview(gene.data = tmt.A682B.entrez.nodups.trim, pathway.id = "zma00592", species = "zma",  bins = 20, limit = 0.5, out.suffix="tmt.A682B")

#Panel B
#merge with mdb to get entrez 
tmt.acs.entrez = merge.data.frame(tmt.acs, mgdb[,c(1,12)], by.x = "Gene", by.y = "v4_gene_model", all = F)
#remove proteins with no entrez or duplicate entrez
tmt.acs.entrez.nodups = tmt.acs.entrez[!is.na(tmt.acs.entrez$Entrez),]
tmt.acs.entrez.nodups = tmt.acs.entrez.nodups[!tmt.acs.entrez.nodups$Entrez=="",]
tmt.acs.entrez.nodups = tmt.acs.entrez.nodups[!duplicated(tmt.acs.entrez.nodups$Entrez),]
#make new dataframe with just l2(acs/B73) values and Entrez as rownames
tmt.acs.entrez.nodups.trim = data.frame("l2.ACS.B73"=tmt.acs.entrez.nodups$l2.ACS.B73)
rownames(tmt.acs.entrez.nodups.trim) = tmt.acs.entrez.nodups$Entrez
#map
tmt.acs.entrez.nodups.trim.pv = pathview(gene.data = tmt.acs.entrez.nodups.trim, pathway.id = "zma00592", species = "zma",  bins = 20, limit = 0.5, out.suffix="tmt.acs")
