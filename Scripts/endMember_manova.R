


load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Data+Phyto.RData")

AllMacs <- c("Brasenia schreberi", "Chara", "Najas flexilis", "Nuphar variegata", "Nymphaea odorata", "Potamogeton amplifolius", "Potamogeton nodosus", "Potamogeton pusillus")
FloatMacs <- c("Brasenia schreberi", "Nuphar variegata", "Nymphaea odorata", "Potamogeton nodosus")
SubMacs <- c("Chara", "Najas flexilis", "Potamogeton amplifolius", "Potamogeton pusillus")
AllTerr <- c("Alder", "Sedge", "Tamarack", "Tree")
LocalTerr <- c("Alder", "Sedge", "Tamarack")


SourceOpts <- list("Floating Macrophytes"=FloatMacs, "Submersed Macrophytes"=SubMacs, "All Terrestrial"=AllTerr, "All Phytoplankton"=c("EpiPhyto", "MetaPhyto"),  "Periphyton"="Periphyton")

EndMembers[,"Taxon"] <- as.character(EndMembers[,"Taxon"])
EndMembers <- EndMembers[EndMembers[,"Taxon"]%in%unlist(SourceOpts),]
EndMembers[,"Resource"] <- as.character(NA)
for(i in 1:length(SourceOpts)){
	t.cat <- names(SourceOpts)[i]
	t.tax <- SourceOpts[[i]]
	EndMembers[EndMembers[,"Taxon"]%in%t.tax,"Resource"] <- t.cat
}
EndMembers[EndMembers[,"Taxon"]%in%,"Taxon"]


iso.vals <- as.matrix(EndMembers[,c("d13C","d15N","dD")])
res.cat <- as.matrix(EndMembers[,"Resource"])
yr.cat <- as.matrix(as.character(EndMembers[,"Year"]))
# summary(manova(iso.vals ~ res.cat*yr.cat), test="Wilks")
summary(manova(iso.vals ~ res.cat), test="Wilks")

summary(manova(iso.vals ~ res.cat), test="Wilks")$stats[1,6]


cat.pairs0 <- expand.grid(names(SourceOpts), names(SourceOpts))
cat.pairs <- cat.pairs0[apply(cat.pairs0, 1, function(x)!any(duplicated(x))),]
row.names(cat.pairs) <- NULL

p.vals <- c()

for(i in 1:nrow(cat.pairs)){
	t.cat.name <- paste(cat.pairs[i,1], cat.pairs[i,2], sep=" & ")
	t.cat.index <- EndMembers[,"Resource"]%in%c(as.character(cat.pairs[i,1]), as.character(cat.pairs[i,2]))
	
	t.iso.vals <- as.matrix(EndMembers[t.cat.index,c("d13C","d15N","dD")])
	t.res.cat <- as.matrix(EndMembers[t.cat.index,"Resource"])
		
	print(t.cat.name)
	print(summary(manova(t.iso.vals ~ t.res.cat), test="Wilks"))
	
	p.vals[i] <- summary(manova(t.iso.vals ~ t.res.cat), test="Wilks")$stats[1,6]
	print(p.vals[i])
	
}



