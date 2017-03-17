# RDB 2017

# Ojbective: Get zooplankton data used in Figure 3 (biplot) of Batt et al. 2015 in Ecosphere
# At the request of Amber Bellamy and Sam Oliver in Feb/March 2017
# (Sorry for being slow, this code is old to me!)

# First, pull in code from the script used to generate the figure
# Second, subset to the stuff that the Eco
# From Scripts/figure3_isotope_biplots.R
# ======================================
# = Pull in Beginning of Figure 3 Code =
# ======================================
load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Data+Phyto.RData")
CexScale <- 0.8

col2012 <- rgb(t(col2rgb("blue")), alpha=200, maxColorValue=256)
col2010 <- c(
				rgb(t(col2rgb("white")), alpha=200, maxColorValue=256), 
				rgb(t(col2rgb("red")), alpha=200, maxColorValue=256)
				)

CombYears <- which(is.element(Data[,"Type"], c("Macrophyte", "Terrestrial"))) #I don't want to make separate points for samples types that I don't really expect to change between years, such as the macrophytes or the terrestrial samples.  I just want to assume that these were the same in the different years, while the consumers, DOM, POM etc. may have changed.
Data[CombYears,"Year"] <- 2010
ExcludeTaxa <- c("POM","Nija","YWP","PKS","BHD2","Mesocyclops", "Pickerell Weed")
ExcludeSamples <- c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")
EndMemberTaxa <- c("Periphyton","Macrophyte","Terrestrial")
InvertTaxa <- c("Zooplankton","Snail")
BaseData <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon, ExcludeTaxa) & is.element(Type, c(EndMemberTaxa, InvertTaxa, "Fish"))))#select=c("Year","Trophic","Taxon","d13C","d15N","dD")

# ===================
# = Subset to Zoops =
# ===================
library(data.table)
bd <- as.data.table(BaseData)
zoopData <- bd[Type=="Zooplankton" | Type=="Chaoborous"]
save(zoopData, file="~/Desktop/Batt2015Ecosphere_zoopData.RData")