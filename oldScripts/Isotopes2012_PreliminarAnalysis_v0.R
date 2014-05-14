rm(list=ls())
graphics.off()

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data/IsotopeData2012")
DataAll <- read.csv("WardIsotopes_2010&2012_4Dec2012.csv")

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")


ZoopD <- subset(DataAll, Type=="Zooplankton" & Taxon!="Chaoborus")
ExploreZoop <- lm(dD~Taxon*as.factor(Year)*Habitat*Week, data=ZoopD)
step(ExploreZoop)




dev.new(width=8, height=7)
par(mfrow=c(2,2), mar=c(2.5,4,1,1))
# =============
# = Calanoids =
# =============
CalanoidD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Calanoid")
summary(lm(dD~as.factor(Year)*Week, data=CalanoidD))

# dev.new()
boxplot(dD~Week*Year*as.character(Habitat), data=CalanoidD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4))
abline(v=8.5)
text(c(3, 11), -240, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
mtext("Calanoid dD (‰)", side=2, line=2.5)
legend("topright", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")

# =============
# = Chaoborus =
# =============
ChaobD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Chaoborus")
summary(lm(dD~as.factor(Year)*Week, data=CalanoidD))

# dev.new()
boxplot(dD~Week*Year, data=ChaobD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),2))
# abline(v=8.5)
# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
mtext("Chaoborus dD (‰)", side=2, line=2.5)
# legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")
text(3, -205, c("Depth Integrated"), font=4, cex=0.9)

# =======
# = POM =
# =======
PomD <- subset(DataAll, Type=="POM" & Habitat!="Hypo" & Habitat!="Littoral")
summary(lm(dD~as.factor(Year)*Week*Habitat, data=PomD))

# dev.new()
boxplot(dD~Week*Year*as.character(Habitat), data=PomD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4))
abline(v=8.5)
text(c(3, 14), -225, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
mtext("POM dD (‰)", side=2, line=2.5)
# legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")


# =============
# = Mesocyclops =
# =============
MesoD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Mesocyclops")
summary(lm(dD~as.factor(Year)*Week, data=MesoD))

# dev.new()
boxplot(dD~Week*Year, data=MesoD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),2))
# abline(v=8.5)
# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
mtext("Mesocyclops dD (‰)", side=2, line=2.5)
# legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")
text(3, -189, "Metalimnion Only", font=4, cex=0.9)



# =============
# = Snail =
# =============
SnailD <- subset(DataAll, Type=="Snail")
summary(lm(dD~as.factor(Year)*Week, data=SnailD))

dev.new()
SnailMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(SnailD[,"Week"])), sort(unique(SnailD[,"Year"])))[,1]]
SnailYears <- as.character(expand.grid(sort(unique(SnailD[,"Week"])), sort(unique(SnailD[,"Year"])))[,2])
BGcols <- c("2010"="#FA807225", "2012"="#3A5FCD25")[SnailYears]
FG2010 <- c("2010"="red", "2012"="blue")[SnailYears]
mtext("SNail dD (‰)", side=2, line=2.5)
# RepYearCol <- length(FishMonths)/2

# dev.new()
boxplot(dD~Week*Year, data=SnailD, col=BGcols, border=FG2010, names=SnailMonths)




# =============
# = Fish =
# =============
FishSpecies <- as.character(unique(subset(DataAll, Type=="Fish" & Tissue=="Whole", select=Taxon))[,])
dev.new(width=9, height=6)
par(mfrow=c(2,3), mar=c(2,3.5,1,1), ps=10, cex=1)

for(i in 1:length(FishSpecies)){
	ThisFish <- FishSpecies[i]
	FishD <- subset(DataAll, Type=="Fish" & Taxon==ThisFish & Tissue=="Whole")
	
	
	FishMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(FishD[,"Week"])), sort(unique(FishD[,"Year"])))[,1]]
	RepYearCol <- length(FishMonths)/2
	
	# dev.new()
	boxplot(dD~Week*Year, data=FishD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=FishMonths)
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	mtext(paste(ThisFish, "dD (‰)"), side=2, line=2.5)
	if(i==2){legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")}
}





# =============
# = DOM =
# =============
DomD <- subset(DataAll, Type=="DOM" & Habitat!="Hypo")
summary(lm(dD~as.factor(Year)*Week, data=DomD))

DOMmonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(DomD[,"Week"])), sort(unique(DomD[,"Year"])))[,1]]
RepYearCol <- length(DOMmonths)/2

dev.new()
boxplot(dD~Week*Year, data=DomD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=FishMonths)
# boxplot(dD~Week*Year, data=DomD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),2))
# abline(v=8.5)
# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
mtext("DOM dD (‰)", side=2, line=2.5)
# legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")
text(3, -189, "Metalimnion Only", font=4, cex=0.9)




# =============
# = Water =
# =============
WaterD <- subset(DataAll, Type=="Water" & Habitat!="Hypo")
summary(lm(dD~as.factor(Year)*Week, data=WaterD))

Watermonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(WaterD[,"Week"])), sort(unique(WaterD[,"Year"])), sort(unique(WaterD[,"Habitat"])))[,1]]
RepYearCol <- length(Watermonths)/2

dev.new()
boxplot(dD~Week*Year*as.character(Habitat), data=WaterD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4))
mtext("Water dD (‰)", side=2, line=2.5)
abline(v=8.5)
text(c(3, 14), -75, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")



# =============
# = Vegetation =
# =============
VegSpecies <- as.character(unique(subset(DataAll, is.element(Type, c("Macrophyte", "Terrestrial")), select=Taxon))[,])
VegKey <- c("Nuphar"="Nuphar", "Brasenia schreberi"="B. schreberi", "Sedge"="Sedge", "Alder"="Alder", "Potamogeton pusillus"="P. pusillus", "Chara"="Chara", "Nymphaea odorata"="N. odorata", "Najas flexilis"="N. flexilis", "Nuphar variegata"="N. variegata", "Tamarack"="Tamarack", "Potamogeton amplifolius"="P. amplif", "Pickerell Weed"="Pick Weed", "Potamogeton nodosus"="P. nodosus")
dev.new(width=7, height=6)
par(mfrow=c(4,3), mar=c(2,3,1,1), ps=9, cex=1)
for(i in 1:length(VegSpecies)){
	ThisVeg <- VegSpecies[i]
	VegD <- subset(DataAll, is.element(Type, c("Macrophyte", "Terrestrial")) & Taxon==ThisVeg)
	
	
	VegMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(VegD[,"Week"])), sort(unique(VegD[,"Year"])))[,1]]
	VegYears <- as.character(expand.grid(sort(unique(VegD[,"Week"])), sort(unique(VegD[,"Year"])))[,2]) #as.character(VegD[,"Year"])
	BGcols <- c("2010"="#FA807225", "2012"="#3A5FCD25")[VegYears]
	FG2010 <- c("2010"="red", "2012"="blue")[VegYears]
	
	# RepYearCol <- length(FishMonths)/2
	
	# dev.new()
	boxplot(dD~Week*Year, data=VegD, col=BGcols, border=FG2010, names=VegMonths)
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	mtext(VegKey[ThisVeg], side=2, line=2)
	if(i==2){legend("topright", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")}
}





solve(a=matrix(data=c(-238,1,-135,1), nrow=2), b=matrix(data=c(-175,1), nrow=2))