rm(list=ls())
graphics.off()

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data/IsotopeData2012")
DataAll <- read.csv("WardIsotopes_2010&2012_4Dec2012.csv")

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")


ZoopD <- subset(DataAll, Type=="Zooplankton" & Taxon!="Chaoborus")
ExploreZoop <- lm(dD~Taxon*as.factor(Year)*Habitat*Week, data=ZoopD)
step(ExploreZoop)



# ========================
# = Plot "invertebrates" =
# ========================
dev.new(width=8, height=7)
par(mfrow=c(2,2), mar=c(2.5,4,1,1))
# CalanoidD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Calanoid")
# ChaobD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Chaoborus")
# MesoD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Mesocyclops")
# SnailD <- subset(DataAll, Type=="Snail")
InverYlim <- NULL #range(subset(DataAll, is.element(Type, c("Zooplankton", "Snail")), select="dD"))
# =============
# = Calanoids =
# =============
CalanoidD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Calanoid")
summary(lm(dD~Week*as.factor(Year)+Habitat+as.factor(Year):Habitat, data=CalanoidD))
#Calanoid story: The seasonal trend of enrichment was much stronger in 2012 than in 2010. This caused the epilimnetic calanoids to finish the season more enriched in 2012 than in 2010, despite the Epi Cals begining the season more depleted in 2012 than in 2010.  Epilimnetic and Metalimnetic calanoids had similar dD throughout 2010.  On the other hand, in 2012 calanoids from the two layers had different starting values (they were more depleted in the metalimnion), but a similar amount of enrichment over the growing season.

boxplot(dD~Week*Year*as.character(Habitat), data=CalanoidD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4), ylim=InverYlim)
abline(v=8.5)
text(c(3, 11), -240, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
mtext("Calanoid dD (‰)", side=2, line=2.5)
legend("topright", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")

# =============
# = Chaoborus =
# =============
ChaobD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Chaoborus")
summary(lm(dD~Week*as.factor(Year), data=ChaobD))
#No difference was observed between 2010 chaoborus and 2012 chaoborus.  Chaoborus in both began both seasons at similar dD, became more enriched during the growing season, and ended the season about 30‰ more enriched in dD.
boxplot(dD~Week*Year, data=ChaobD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),2), ylim=InverYlim)
mtext("Chaoborus dD (‰)", side=2, line=2.5)
text(3, -240, c("Depth Integrated"), font=4, cex=0.9)

# =============
# = Mesocyclops =
# =============
MesoD <- subset(DataAll, Type=="Zooplankton" & Taxon=="Mesocyclops")#not enough points to detect seasonal differences.  Years appear similar.
boxplot(dD~Week*Year, data=MesoD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),2), ylim=InverYlim)
mtext("Mesocyclops dD (‰)", side=2, line=2.5)
text(3, -240, "Metalimnion Only", font=4, cex=0.9)

# =============
# = Snail =
# =============
SnailD <- subset(DataAll, Type=="Snail")
#Snails were very similar within and between years, with the exception that in June of 2012 the snails were ~25‰ more depleted than at any other sampling period.
SnailMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(SnailD[,"Week"])), sort(unique(SnailD[,"Year"])))[,1]]
SnailYears <- as.character(expand.grid(sort(unique(SnailD[,"Week"])), sort(unique(SnailD[,"Year"])))[,2])
BGcols <- c("2010"="#FA807225", "2012"="#3A5FCD25")[SnailYears]
FG2010 <- c("2010"="red", "2012"="blue")[SnailYears]
boxplot(dD~Week*Year, data=SnailD, col=BGcols, border=FG2010, names=SnailMonths, ylim=InverYlim)
mtext("SNail dD (‰)", side=2, line=2.5)




dev.new(width=8, height=7)
par(mfrow=c(2,2), mar=c(2.5,4,1,1))
# =======
# = POM =
# =======
PomD <- subset(DataAll, Type=="POM" & Habitat!="Hypo" & Habitat!="Littoral")
summary(lm(dD~relevel(as.factor(Week),4)*as.factor(Year)+Habitat+relevel(as.factor(Week),4):Habitat, data=PomD))#Shows that the week-to-week (seasonal) pattern was not different between years (is the same conclusion if you remove the week:habitat interaction terms)
summary(lm(dD~relevel(as.factor(Week),4)+as.factor(Year)*Habitat, data=PomD))
#POM Story: POM was significantly more enriched in 2012 than in 2010 (p=5.36e-07), especially in the metalimnion (p=0.01106).  Furthermore, there was a strong month-to-month variation in POM dD that was consistent among layers and years--- POM tended to become more depleted over the season, until the last month when it became slightly more enriched than the previous month, but still more depleted than the first 2 months (think reverse check mark).
boxplot(dD~Week*Year*as.character(Habitat), data=PomD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4))
abline(v=8.5)
text(c(3, 14), -225, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
mtext("POM dD (‰)", side=2, line=2.5)

# =============
# = Periphyton =
# =============
PeriD <- subset(DataAll, Type=="Periphyton")
summary(lm(dD~as.factor(Year), data=PeriD))
#Periphyton was ~10.417‰ more depleted in 2012 than in 2010 (p=0.00756).
PeriMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(PeriD[,"Week"])), sort(unique(PeriD[,"Year"])))[,1]]
RepYearCol <- length(PeriMonths)/2
boxplot(dD~Week*Year, data=PeriD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=PeriMonths)
mtext("Periphyton dD (‰)", side=2, line=2.5)

# =============
# = DOM =
# =============
DomD <- subset(DataAll, Type=="DOM" & Habitat!="Hypo")
summary(lm(dD~as.factor(Year), data=DomD))
#DOM was ~40.77‰ more depleted in 2012 than in 2010.  The grand mean and Year-to-Year variability (i.e., the intercept and the Year effect) explain 96.43% of the observed variation in DOM dD.
DomMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(DomD[,"Week"])), sort(unique(DomD[,"Year"])), sort(unique(DomD[,"Habitat"])))[,1]]
RepYearCol <- length(DomMonths)/2
DomYears <- as.character(expand.grid(sort(unique(DomD[,"Week"])), sort(unique(DomD[,"Year"])), sort(unique(DomD[,"Habitat"])))[,2])
BGcols <- c("2010"="#FA807225", "2012"="#3A5FCD25")[DomYears]
FG2010 <- c("2010"="red", "2012"="blue")[DomYears]
boxplot(dD~Week*Year*as.character(Habitat), data=DomD, col=BGcols, border=FG2010, names=DomMonths)
abline(v=6.5)
text(c(1.75, 8.25), -164, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
mtext("DOM dD (‰)", side=2, line=2.5)

# =============
# = Water =
# =============
WaterD <- subset(DataAll, Type=="Water" & Habitat!="Hypo")
summary(lm(dD~as.factor(Year)+Week+Habitat, data=WaterD))
#Water: Water deuterium was markedly more enriched in 2012 than in 2010 (8.789‰, p=2.79e-08).  In both years and in both layers, water became more enriched over time (3.1‰ per month, p=3.94e-06).  In both years, the metalimnion was slightly more depleted than the epilimnion (3.9‰ more depleted, p=0.000421)
Watermonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(WaterD[,"Week"])), sort(unique(WaterD[,"Year"])), sort(unique(WaterD[,"Habitat"])))[,1]]
RepYearCol <- length(Watermonths)/2
boxplot(dD~Week*Year*as.character(Habitat), data=WaterD, col=c(rep("#FA807225",4), rep("#3A5FCD25",4)), border=c(rep("red",4), rep("blue",4)), names=rep(c("May","Jun","Jul","Aug"),4))
mtext("Water dD (‰)", side=2, line=2.5)
abline(v=8.5)
text(c(3, 14), -75, c("Epilimnion", "Metalimnion"), font=4, cex=0.9)
legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")




# =============
# = Fish =
# =============
AllFishD <- subset(DataAll, Type=="Fish" & Tissue=="Whole" & Taxon!="YWP") #& Taxon!="PKS"
summary(lm(dD~as.factor(Year)*Taxon+Week, data=AllFishD))
#In general, Ward fish tend to become more enriched in deuterium over the growing season (~5.8‰ per month, p=0.00108; Pumpkinseeds excluded from this regression due to little informatino on 2012 seasonality).  No species-specific seasonal trend in dD was detectable in our data. Brown bullheads, dace, and fathead minnows were slightly more enriched in 2012 than in 2010 (12‰, 5‰, and 5‰, respectively; the change in dD between 2010 and 2012 is not significantly different for these three species), whereas the central mud minnow and pumpkinseed sunfish were depleted in 2012 than in 2010,  (4%  and 19‰ more depleted, respectively; significantly different from the bullhead, p=0.0279 and p=0.000081, respectively).  The central mud minnow was the most enriched in 2H. The central mud minnow had the most enriched dD of any of these 5 species of Ward fish, and was one of two species that was more depleted in deuterium in 2012 than in 2010.  The pumpkinseed was the second most enriched in dD, and had the greatest decrease in dD between 2010 and 2012.

FishSpecies <- as.character(unique(subset(DataAll, Type=="Fish" & Tissue=="Whole", select=Taxon))[,])
dev.new(width=9, height=6)
par(mfrow=c(2,3), mar=c(2,3.5,1,1), ps=10, cex=1)
for(i in 1:length(FishSpecies)){
	ThisFish <- FishSpecies[i]
	FishD <- subset(DataAll, Type=="Fish" & Taxon==ThisFish & Tissue=="Whole")
	FishMonths <- c("May","Jun","Jul","Aug")[expand.grid(sort(unique(FishD[,"Week"])), sort(unique(FishD[,"Year"])))[,1]]
	RepYearCol <- length(FishMonths)/2
	boxplot(dD~Week*Year, data=FishD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=FishMonths)
	mtext(paste(ThisFish, "dD (‰)"), side=2, line=2.5)
	if(i==2){legend("topleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")}
}


# =============
# = Vegetation =
# =============
AllVegD <- subset(DataAll, Type=="Macrophyte") #& Taxon!="PKS"
summary(lm(dD~as.factor(Year)+Taxon*Week, data=AllVegD))
#No differences between years was detected for any species.  However, there really isn't enough 2012 data for any of the species to detect a seasonal trend, b/c each species only had 2 2012 samples sent in for analysis.  There really isn't anything too striking here.
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
	boxplot(dD~Week*Year, data=VegD, col=BGcols, border=FG2010, names=VegMonths)
	mtext(VegKey[ThisVeg], side=2, line=2)
	if(i==2){legend("topright", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")}
}





solve(a=matrix(data=c(-238,1,-135,1), nrow=2), b=matrix(data=c(-175,1), nrow=2))



