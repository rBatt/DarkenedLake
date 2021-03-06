#Ryan Batt
#23 April 2011
#What is POM made of?
#Given POM, what is a consumer made of?
#The purpose of this script is first calculate the constituent components of Ward POM from the summer of 2010.
#Next, I will determine the composition of a consumer.
#I begin with the simplifying assumption that POM is made of terrestrial (sedge and alder) and phytoplankton sources.
#I will also assume that the consumer is eating some combination of the following four sources: 1) Epi Phyto 2) Meta Phyto 3) Equal parts of Alder, Sedge, Tree 4) DOC
#Version 5:
	#Intended to be the final version
	#Does not do a massive simulation through possible source combinations
	#Looks at 2 possible source combinations: 1 with the phytoplankton split into Epi and Meta and the macrophytes and periphyton grouped, and the other with the phytos grouped but periphyton in a group separate from the macrophytes
	#Previous analyses had forgotten to remove the Watershield data point that was a "stem" (I think, anyway).  
	#I may need to treat the "Tree" variance difference in the future, because this is actually adding another layer of nesting within a "source"
	#This version will use a larger number of chains and longer chain lengths, and will do cross-validation for the density plots
	#The plots should look better overall
	#There are several samples which will be automatically excluded from analysis:
		#All the Meta B POC samples-- Meta B sampling was thwarted throughout the season by a massive patch of Chara
		#The Hypo DOC sample-- it has an N signature that is quite different from the others
		#The "Nija" sample b/c there was only 1 sample of it
		#The watershield sample that was just a stem-- its deuterium was different from the "leaf" samples
		#The DIC data has not been edited for these weird Meta B and Hypo E samples, but those values are not used in this at all
#Version 5.1: 
	#Commented out that bullshit with the terrestrial variance being copied for the pelagic epilimnion and pelagic metalimnion... zomg.
#Version 7.0: Changed the terrestrial end member to not include DOM (DOC).  Also, I later changed the graphing of the phytoplankton posteriors to round to one less digit for carbon-- this is to only have 3 sig figs, and also to make sure the 13C peak for the epi didn't overlap with the estimate
#Version 8.0: Including a new data file for the isotopes, which now includes additional tree data.  Averages for each species are taken from the Cascade database.  For the trees, there are only one or two samples (if 2, it's just analytical replicates) per species for C/N, whereas there are a ton of dueterium samples typically.  The sample number refers to the sample number for the C/N data.
	#Found an error where the number of consumers in the ConsMix function was calculted as the number of columns, but it should have been the number of rows
#Version 8.1: 
	#Run with DOM as it's own "terrestrial" source
#Version 8.2:
	#Run with the "terrestrial" source as Alder, Sedge, and DOM
#Version 8.3:
	#I am reformatting the Figures according to the editorial marks that I received on 10-May-2012
	
#Vesion 0.0.0 (10-Jan-2013): I am starting the versioning of this program over for the analysis of isotope data post-Aquashade manipulation
	#Made changes in several places so that the analysis would be year-specific
	#Automated the placement of the version and the figure folder name into several strings, given their specification near the beginning of this script
	#Option to specify # of iterations
#Version 0.1.0 (11-Jan-2013): The previous version "worked" fine with the new data, but I kept estimating benthic contribution to the zooplankton, which I don't believe is actually happening.  In an effort to correct this, I am changing this script to allow for consumer-specific groupings of sources.  I don't want to just remove the option for zooplankton (etc.) to eat periphyton, I just think that this would be a less likely diet choice if the other 3 options were more appropriate.  Regardless, the idea is to have the option to tailor the resource groupings to the specific consumer being analyzed.
#Version 0.3.0 (03-Feb-2013): I am changing the diet for the zoops.  Grouping the algae together, all macrohpytes, terrestrial, and periphyton.  This says that they are eating some periphyton and that the chaobs didn't really increase in allochthony.
#Version 0.3.1 (03-Feb-2013): I am changing the diet for the zoops.  Grouping the algae together, doing submerused, floating, terrestrial.
#Version 0.3.3 (11-April-2013): I commented-out the per-month analysis and plotting. Added the working directory to the bugs() functions.
#Version 0.4.0 (11-April-2013): Trying to mess w/ the number of sources

rm(list=ls())
graphics.off()

library(R2WinBUGS)
library(R2jags)

Iterations <- 2000

if(Sys.info()["sysname"]=="Windows"){
	windowsFonts(Times=windowsFont("TT Times New Roman"))
	BUGSWorkDir <- NULL
}else{
	BUGSWorkDir <- "~/.wine/drive_c/temp/Rtmp"
}

#Select the top 2 if on Snow Leopard, the bottom 2 if on Leopard, and the selection doesn't matter if on a PC
WINE="/Applications/Darwine/Wine.bundle/Contents/bin/wine"
WINEPATH="/Applications/Darwine/Wine.bundle/Contents/bin/winepath"
# WINEPATH="/opt/local/bin/winepath"
# WINE="/opt/local/bin/wine"



# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2010Analysis")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/ConsMix.R")
DataRaw <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/IsotopeData2012/WardIsotopes_2010&2012_17Jan2013.csv", header=TRUE)
Data <- subset(DataRaw, Taxon!="Nija" & !is.element(SampleID, c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405", "P-1244")) & is.na(FishID))
Months <- c("May", "Jun", "Jul", "Aug")




#Calculate the algal end member from POM
TSources <- c("Alder", "Sedge", "Tamarack", "Tree")#, "Tamarack") #c("Alder", "Sedge", "DOM")

#Signature of the terrestrial source
nTS <- length(TSources)
TMeans <- data.frame("d13C"=rep(NA,nTS),"d15N"=rep(NA,nTS),"dD"=rep(NA,nTS), row.names=TSources)
TVars <- data.frame("d13C"=rep(NA,nTS),"d15N"=rep(NA,nTS),"dD"=rep(NA,nTS), row.names=TSources)
#Td13C_aov 
Td15NObs <- data.frame()
TdDObs <- data.frame() #matrix(ncol=nTS, dimnames=list(NULL,TSources))
for(i in 1:length(TSources)){
	TMeans[i,] <- apply(subset(Data, Taxon==TSources[i], select=c("d13C","d15N","dD")),2,mean)
}
dCNH_Terr_Mu <- apply(TMeans, 2, mean)
dCNH_Terr_Var <- data.frame("d13C"=NA, "d15N"=NA, "dD"=NA)

if(nTS>1){
	Temp_d13C_aov <- anova(lm(d13C ~ Taxon, data=subset(Data, is.element(Taxon, TSources), select=c("Taxon","d13C"))))
	if(Temp_d13C_aov$Pr[1] <= 0.1){
		dCNH_Terr_Var["d13C"] <- sum(Temp_d13C_aov$Mean)
	}else{
		dCNH_Terr_Var["d13C"] <- Temp_d13C_aov$Mean[2]
	}
		
	Temp_d15N_aov <- anova(lm(d15N ~ Taxon, data=subset(Data, is.element(Taxon, TSources), select=c("Taxon","d15N"))))
	if(Temp_d15N_aov$Pr[1] <= 0.1){
		dCNH_Terr_Var["d15N"] <- sum(Temp_d15N_aov$Mean)
	}else{
		dCNH_Terr_Var["d15N"] <- Temp_d15N_aov$Mean[2]
	}
		
	Temp_dD_aov <- anova(lm(dD ~ Taxon, data=subset(Data, is.element(Taxon, TSources), select=c("Taxon","dD"))))
	if(Temp_dD_aov$Pr[1] <= 0.1){
		dCNH_Terr_Var["dD"] <- sum(Temp_dD_aov$Mean)
	}else{
		dCNH_Terr_Var["dD"] <- Temp_dD_aov$Mean[2]
	}
}else{
	dCNH_Terr_Var <- apply(subset(Data, is.element(Taxon, TSources), select=c("d13C", "d15N", "dD")), 2, var)
}
#Define the Terr objects to be used in the POM Mixture portion of the BUGS model

#**************************************
T_dX <- as.numeric(dCNH_Terr_Mu)
T_dX_Var <- as.numeric(dCNH_Terr_Var)
#**************************************




for(YearMix in c(2010, 2012)){
	# TODO The water will need to be defined by year.  Either stored in a higher dimensional object, or have separate objects for each year.
	Water_dD_Mu <- mean(subset(Data, Type=="Water" & Year==YearMix, select="dD")[,])
	Water_dD_Var <- var(subset(Data, Type=="Water" & Year==YearMix, select="dD")[,])

	#Calculate Epi phyto deuterium prior from water dD
	dD_Water_Epi <- subset(Data, Type=="Water" & Habitat=="Epi" & Year==YearMix, select="dD")[,]
	dD_Water_Adj <- mean(c(-152.8, -172.4))#Fractionation range reported in Solomon et al. 2011
	dD_Phyto_Epi_Mu <- mean(dD_Water_Epi + dD_Water_Adj)
	#From Solomon et al. 2011 Appendix A: alpha phyto-water = mean ± sd = 0.84 ± 0.008; qnorm(p=.025, mean=-231.945, sd=5); var=25 it should have been ~70.. ask Grace.
	dD_Phyto_Epi_Var <- var(dD_Water_Epi) + 25#variance of water + fractionation = variance of Phyto
	dD_Phyto_Epi_Shape <- dD_Phyto_Epi_Var*0.1#dD_Phyto_Var~dgamma(shape,rate); shape when rate==0.1
	#Signature of the Epi POM mixture
	dCNH_POM_Epi <- subset(Data, Type=="POM" & Habitat=="Epi" & Year==YearMix, select=c("d13C","d15N","dD"))
	POM_dX_Epi_Obs <- matrix(data=c(dCNH_POM_Epi[,1], dCNH_POM_Epi[,2], dCNH_POM_Epi[,3]), ncol=3)
	POM_dX_Epi_Var <- apply(dCNH_POM_Epi, 2, var)
	nPOM_Epi <- length(POM_dX_Epi_Obs[,1])

	#Same POM and phyto calcs for Meta
	#Calculate Algal deuterium prior from water dD
	dD_Water_Meta <- subset(Data, Type=="Water" & Habitat=="Meta" & Year==YearMix, select="dD")[,]
	dD_Phyto_Meta_Mu <- mean(dD_Water_Meta + dD_Water_Adj)
	#From Solomon et al. 2011 Appendix A: alpha phyto-water = mean ± sd = 0.84 ± 0.008; qnorm(p=.025, mean=-231.945, sd=5); var=25
	dD_Phyto_Meta_Var <- var(dD_Water_Meta) + 25#variance of water + variance of fractionation = variance of Phyto
	dD_Phyto_Meta_Shape <- dD_Phyto_Meta_Var*0.1#dD_Phyto_Var~dgamma(shape,rate); shape when rate==0.1

	#Signature of the Meta POM mixture
	dCNH_POM_Meta <- subset(Data, Type=="POM" & Habitat=="Meta" & Year==YearMix, select=c("d13C","d15N","dD"))
	POM_dX_Meta_Obs <- matrix(data=c(dCNH_POM_Meta[,1], dCNH_POM_Meta[,2], dCNH_POM_Meta[,3]), ncol=3)
	POM_dX_Meta_Var <- apply(dCNH_POM_Meta, 2, var)
	nPOM_Meta <- length(POM_dX_Meta_Obs[,1])

	#Run BUGS Part 1: Using POM, calculate the isotopic signatures of epilimnetic and metalimnetic phytoplankton
	SupplyBUGS_pt1 <- list(T_dX, T_dX_Var, dD_Phyto_Epi_Mu, dD_Phyto_Epi_Shape, POM_dX_Epi_Obs, nPOM_Epi, dD_Phyto_Meta_Mu, dD_Phyto_Meta_Shape, POM_dX_Meta_Obs, nPOM_Meta)
	names(SupplyBUGS_pt1) <- strsplit(c("T_dX, T_dX_Var, dD_Phyto_Epi_Mu, dD_Phyto_Epi_Shape, POM_dX_Epi_Obs, nPOM_Epi, dD_Phyto_Meta_Mu, dD_Phyto_Meta_Shape, POM_dX_Meta_Obs, nPOM_Meta"), split=", ")[[1]]
	ParamBUGS_pt1 <- c("f", "P_dC_Epi", "P_dN_Epi", "P_dD_Epi", "P_dC_Epi_Var", "P_dN_Epi_Var", "P_dD_Epi_Var",  "P_dC_Meta", "P_dN_Meta", "P_dD_Meta", "P_dC_Meta_Var", "P_dN_Meta_Var", "P_dD_Meta_Var", "residSd")
	BUGSfile_pt1 <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/mix_Cons_Mixture_Ward2010_v2_pt1.bug"
	bugsOut_pt1 <- jags(SupplyBUGS_pt1, inits=NULL, ParamBUGS_pt1, BUGSfile_pt1, n.chains=8, n.iter=Iterations)



	#Extract and name relevant information concerning epilimnetic and metalimnetic phytoplankton
	#**************************************
	P_dX_Epi <- c(bugsOut_pt1$BUGSoutput$mean$P_dC_Epi, bugsOut_pt1$BUGSoutput$mean$P_dN_Epi, bugsOut_pt1$BUGSoutput$mean$P_dD_Epi)
	P_dX_Epi_Var <- c(bugsOut_pt1$BUGSoutput$mean$P_dC_Epi_Var, bugsOut_pt1$BUGSoutput$mean$P_dN_Epi_Var, bugsOut_pt1$BUGSoutput$mean$P_dD_Epi_Var)
	P_dX_Meta <- c(bugsOut_pt1$BUGSoutput$mean$P_dC_Meta, bugsOut_pt1$BUGSoutput$mean$P_dN_Meta, bugsOut_pt1$BUGSoutput$mean$P_dD_Meta)
	P_dX_Meta_Var <- c(bugsOut_pt1$BUGSoutput$mean$P_dC_Meta_Var, bugsOut_pt1$BUGSoutput$mean$P_dN_Meta_Var, bugsOut_pt1$BUGSoutput$mean$P_dD_Meta_Var)

	Sim_P_dX_Epi_Obs <- as.data.frame(matrix(data=rep(rnorm(n=nPOM_Epi),3), ncol=3, byrow=FALSE))

	Sim_P_dX_Epi_Obs[,1] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Epi"], size=nPOM_Epi)
	Sim_P_dX_Epi_Obs[,2] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Epi"], size=nPOM_Epi)
	Sim_P_dX_Epi_Obs[,3] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dD_Epi"], size=nPOM_Epi)

	Sim_P_dX_Epi_Obs <- (Sim_P_dX_Epi_Obs-as.data.frame(matrix(data=rep(apply(Sim_P_dX_Epi_Obs,2,mean), nPOM_Epi), ncol=3, byrow=TRUE)))/as.data.frame(matrix(data=rep(apply(Sim_P_dX_Epi_Obs,2,sd), nPOM_Epi), ncol=3, byrow=TRUE))
	Sim_P_dX_Epi_Obs[,1] <- Sim_P_dX_Epi_Obs[,1]*sqrt(P_dX_Epi_Var[1])+P_dX_Epi[1]
	Sim_P_dX_Epi_Obs[,2] <- Sim_P_dX_Epi_Obs[,2]*sqrt(P_dX_Epi_Var[2])+P_dX_Epi[2]
	Sim_P_dX_Epi_Obs[,3] <- Sim_P_dX_Epi_Obs[,3]*sqrt(P_dX_Epi_Var[3])+P_dX_Epi[3]

	colnames(Sim_P_dX_Epi_Obs) <- c("d13C","d15N","dD")
	Sim_P_dX_Epi_Obs <- cbind("Taxon"=rep("EpiPhyto",nPOM_Epi), Sim_P_dX_Epi_Obs)

	Sim_P_dX_Meta_Obs <- as.data.frame(matrix(data=rep(rnorm(n=nPOM_Meta),3), ncol=3, byrow=FALSE))

	Sim_P_dX_Meta_Obs[,1] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Meta"], size=nPOM_Meta)
	Sim_P_dX_Meta_Obs[,2] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Meta"], size=nPOM_Meta)
	Sim_P_dX_Meta_Obs[,3] <- sample(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dD_Meta"], size=nPOM_Meta)

	Sim_P_dX_Meta_Obs <- (Sim_P_dX_Meta_Obs-as.data.frame(matrix(data=rep(apply(Sim_P_dX_Meta_Obs,2,mean), nPOM_Meta), ncol=3, byrow=TRUE)))/as.data.frame(matrix(data=rep(apply(Sim_P_dX_Meta_Obs,2,sd), nPOM_Meta), ncol=3, byrow=TRUE))
	Sim_P_dX_Meta_Obs <- (Sim_P_dX_Meta_Obs-apply(Sim_P_dX_Meta_Obs,2,mean))/apply(Sim_P_dX_Meta_Obs,2,sd)
	Sim_P_dX_Meta_Obs[,1] <- Sim_P_dX_Meta_Obs[,1]*sqrt(P_dX_Meta_Var[1])+P_dX_Meta[1]
	Sim_P_dX_Meta_Obs[,2] <- Sim_P_dX_Meta_Obs[,2]*sqrt(P_dX_Meta_Var[2])+P_dX_Meta[2]
	Sim_P_dX_Meta_Obs[,3] <- Sim_P_dX_Meta_Obs[,3]*sqrt(P_dX_Meta_Var[3])+P_dX_Meta[3]

	colnames(Sim_P_dX_Meta_Obs) <- c("d13C","d15N","dD")
	Sim_P_dX_Meta_Obs <- cbind("Taxon"=rep("MetaPhyto",nPOM_Meta), Sim_P_dX_Meta_Obs)
	#**************************************


	#*****************************************************
	#Begin for consumers and their respective sources
	#*****************************************************

	if(YearMix==2010){
		Cons <- c("Calanoid", "Chaoborus", "Helisoma trivolvis", "FHM", "DAC", "CMM", "BHD1", "BHD2") #, "PKS", "FHM", "YWP", "CMM", "BHD", "Mesocyclops", "DAC")
		TL <- c(1, 2, 1, 2, 2, 2, 2, 3) #, 2, 2, 3, 2.5, 2.5, 1.5, 1)
		GraphTitle <- c("Skistodiaptomus oregonensis", "Chaoborus spp.", "Helisoma trivolvis", "Pimephales promelas", "Phoxinus spp.", "Umbra limi", "Small Ameiurus melas", "Big Ameiurus melas") #, "Lepomis gibbosus", "Pimephales promelas", "Perca flavescens", "Umbra limi", "Ameiurus melas",  "Mesocyclops spp.", "Phoxinus spp.")
	}else{
		Cons <- c("Calanoid", "Chaoborus", "Helisoma trivolvis", "FHM", "DAC", "CMM", "BHD1", "BHD2") #, "PKS", "FHM", "CMM", "BHD", "Mesocyclops", "DAC")
		TL <- c(1, 2, 1, 2, 2, 2, 2, 3) #, 2, 2, 2.5, 2.5, 1.5, 1)
		GraphTitle <- c("Skistodiaptomus oregonensis", "Chaoborus spp.", "Helisoma trivolvis", "Pimephales promelas", "Phoxinus spp.", "Umbra limi", "Small Ameiurus melas", "Big Ameiurus melas") #, "Lepomis gibbosus", "Pimephales promelas", "Umbra limi", "Ameiurus melas",  "Mesocyclops spp.", "Phoxinus spp.")
	}
	
	
	AllMacs <- c("Brasenia schreberi", "Chara", "Najas flexilis", "Nuphar variegata", "Nymphaea odorata", "Potamogeton amplifolius", "Potamogeton nodosus", "Potamogeton pusillus")
	FloatMacs <- c("Brasenia schreberi", "Nuphar variegata", "Nymphaea odorata", "Potamogeton nodosus")
	SubMacs <- c("Chara", "Najas flexilis", "Potamogeton amplifolius", "Potamogeton pusillus")
	AllTerr <- c("Alder", "Sedge", "Tamarack", "Tree")
	LocalTerr <- c("Alder", "Sedge", "Tamarack")


	SourceOpts <- list("All Macrophytes"=AllMacs, "Floating Macrophytes"=FloatMacs, "Submersed Macrophytes"=SubMacs, "All Terrestrial"=AllTerr, "Local Terrestrial"=LocalTerr, "All Phytoplankton"=c("EpiPhyto", "MetaPhyto"), "Epi. Phytoplankton"="EpiPhyto", "Meta. Phytoplankton"="MetaPhyto", "DOM"="DOM", "Periphyton"="Periphyton")
	
	#insert the combinations code here
	# AllCombinations <- t(combn(c("All Terrestrial", "Epi. Phytoplankton", "Meta. Phytoplankton", "All Phytoplankton", "Submersed Macrophytes", "Floating Macrophytes", "All Macrophytes", "Periphyton", NA, NA), m=4))
	
	AllCombinations <- t(combn(c("All Terrestrial", "All Phytoplankton", "Submersed Macrophytes", "Floating Macrophytes", "All Macrophytes", "Periphyton", NA, NA), m=4))

	# Checks <- function(x){
		# c1 <- is.element("All Terrestrial", x) #YES
		# c2 <- TRUE #any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton", "All Phytoplankton"), x)) #YES
		# c3 <- any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton"), x)) & is.element("All Phytoplankton", x) #NO
		# c4 <- any(is.element(c("Floating Macrophytes", "Submersed Macrophytes"), x)) & is.element("All Macrophytes", x) #NO
		# c5 <- is.element("Epi. Phytoplankton", x) == is.element("Meta. Phytoplankton", x) #YES
		# return(all(c1, c2, !c3, !c4, c5))
	# }
	
		Checks <- function(x){
			c1 <- is.element("All Terrestrial", x) #YES
			# c2 <- TRUE #any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton", "All Phytoplankton"), x)) #YES
			# c3 <- any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton"), x)) & is.element("All Phytoplankton", x) #NO
			c4 <- any(is.element(c("Floating Macrophytes", "Submersed Macrophytes"), x)) & is.element("All Macrophytes", x) #NO
			# c5 <- is.element("Epi. Phytoplankton", x) == is.element("Meta. Phytoplankton", x) #YES
			return(all(c1, !c4))
	}

	Checked <- apply(AllCombinations, MARGIN=1, FUN=Checks)

	FewerCombns <- AllCombinations[Checked,]
	Combinations <- unique(FewerCombns)
	ComboList <- as.list(as.data.frame(t(Combinations)))

	ConsChoices <- list(
					"Calanoid"=ComboList,
					"Chaoborus"=ComboList,
					"Helisoma trivolvis"=ComboList,
					"FHM"=ComboList,
					"DAC"=ComboList,
					"CMM"=ComboList,
					"BHD1"=ComboList,
					"BHD2"=ComboList
					)

	SourceData <- subset(Data, Trophic==0 & Taxon!="POM" & Year==YearMix | is.element(Type, c("Macrophyte", "Terrestrial")))
	SourceTaxa <- as.character(unique(SourceData[,"Taxon"]))
	Source_Means <- matrix(ncol=3, nrow=length(SourceTaxa), dimnames=list(SourceTaxa,NULL))
	Source_Vars <- matrix(ncol=3, nrow=length(SourceTaxa), dimnames=list(SourceTaxa,NULL))
	for(i in 1:length(SourceTaxa)){
		Source_Means[i,] <- apply(subset(SourceData, Taxon==SourceTaxa[i], select=c("d13C","d15N","dD")), 2, mean)
		Source_Vars[i,] <- apply(subset(SourceData, Taxon==SourceTaxa[i], select=c("d13C","d15N","dD")), 2, var)
	}
	Source_Means <- rbind(Source_Means, "EpiPhyto"=P_dX_Epi, "MetaPhyto"=P_dX_Meta)
	Source_Vars <- rbind(Source_Vars, "EpiPhyto"=P_dX_Epi_Var, "MetaPhyto"=P_dX_Meta_Var)

	# nSrcs <- length(SourceNames[[f_Src]])
	SourceData_dX_Obs <- SourceData[,c("Taxon","d13C","d15N","dD")]
	SourceData_dX_Obs <- rbind(SourceData_dX_Obs, Sim_P_dX_Epi_Obs, Sim_P_dX_Meta_Obs)
	
	if(YearMix==2010){
		EndMembers <- data.frame("Year"=YearMix, SourceData_dX_Obs)
	}else{
		EndMembers0 <- data.frame("Year"=YearMix, SourceData_dX_Obs)
		EndMembers <- rbind(EndMembers, EndMembers0)
	}
	for(g_Cons in 1:length(Cons)){
		TempoCons <- Cons[g_Cons]
		SourceNames <- ConsChoices[[TempoCons]]
		# FirstSources <- list(SourceOpts[[SourceNames[[1]][1]]], SourceOpts[[SourceNames[[2]][1]]])
		# SecondSources <- list(SourceOpts[[SourceNames[[1]][2]]], SourceOpts[[SourceNames[[2]][2]]])
		# ThirdSources <- list(SourceOpts[[SourceNames[[1]][3]]], SourceOpts[[SourceNames[[2]][3]]])
		# FourthSources <- list(SourceOpts[[SourceNames[[1]][4]]], SourceOpts[[SourceNames[[2]][4]]])

		for(f_Src in 1:length(ComboList)){
			# ComboList[[23]] <- ComboList[[23]][which(!is.na(ComboList[[23]]))]
			SourceNames[[f_Src]] <- as.character(SourceNames[[f_Src]][which(!is.na(SourceNames[[f_Src]]))])
			nSrcs <- length(SourceNames[[f_Src]])
			
			Source1 <- SourceOpts[[SourceNames[[f_Src]][1]]] #FirstSources[[f_Src]]
			Source2 <- SourceOpts[[SourceNames[[f_Src]][2]]]
			if(nSrcs>=3){Source3 <- SourceOpts[[SourceNames[[f_Src]][3]]]}
			if(nSrcs>=4){Source4 <- SourceOpts[[SourceNames[[f_Src]][4]]]}
		
			for(i in 1:nSrcs){
				TempName_Source <- paste("Source", i, sep="")
				TempName_Mean <- paste(paste("Source", paste(i, "_Mean", sep=""), sep=""))
				TempName_Var <- paste(paste("Source", paste(i, "_Var", sep=""), sep=""))
				if(length(get(TempName_Source))>1){assign(TempName_Mean, apply(Source_Means[get(TempName_Source),], 2, mean))}else{assign(TempName_Mean, Source_Means[get(TempName_Source),])}

				if(length(get(TempName_Source))>1){
		
					assign(TempName_Var, data.frame("d13C"=NA, "d15N"=NA, "dD"=NA))#This is to clear the temporary data frame at the beginning of each loop

					Temp_d13C_aov <- anova(lm(d13C ~ Taxon, data=subset(SourceData_dX_Obs, is.element(Taxon, get(TempName_Source)), select=c("Taxon","d13C"))))
					if(Temp_d13C_aov$Pr[1] <= 0.1){
						Temp_d13C_Var <- sum(Temp_d13C_aov$Mean)
					}else{
							Temp_d13C_Var <- Temp_d13C_aov$Mean[2] 
					}
				
					Temp_d15N_aov <- anova(lm(d15N ~ Taxon, data=subset(SourceData_dX_Obs, is.element(Taxon, get(TempName_Source)), select=c("Taxon","d15N"))))
					if(Temp_d15N_aov$Pr[1] <= 0.1){
						Temp_d15N_Var <- sum(Temp_d15N_aov$Mean)
					}else{
							Temp_d15N_Var <- Temp_d15N_aov$Mean[2]
					}
				
					Temp_dD_aov <- anova(lm(dD ~ Taxon, data=subset(SourceData_dX_Obs, is.element(Taxon, get(TempName_Source)), select=c("Taxon","dD"))))
					if(Temp_dD_aov$Pr[1] <= 0.1){
						Temp_dD_Var <- sum(Temp_dD_aov$Mean) 
					}else{
						Temp_dD_Var <- Temp_dD_aov$Mean[2] 
					}
				
					assign(TempName_Var, c(Temp_d13C_Var, Temp_d15N_Var, Temp_dD_Var))
				}else{
					assign(TempName_Var, Source_Vars[get(TempName_Source),])
				}	
			}#Finish the loop that handles each source (1 through 4) one at a time for this particular set of sources for this consumer
			#Then collect the source means and variances from the previous loop; the following could have been condensed into previous loop.
			Srcs_dX_Ward <- c()
			Srcs_dX_Var_Ward <- c()
			for(i in 1:nSrcs){
				TempName_Mean <- paste(paste("Source", paste(i, "_Mean", sep=""), sep=""))
				TempName_Var <- paste(paste("Source", paste(i, "_Var", sep=""), sep=""))
				Srcs_dX_Ward <- cbind(Srcs_dX_Ward, get(TempName_Mean))
				Srcs_dX_Var_Ward <- cbind(Srcs_dX_Var_Ward, get(TempName_Var))
			}
			if(YearMix==2010 & all(SourceNames[[f_Src]]%in%c("All Terrestrial", "All Macrophytes", "All Phytoplankton","Periphyton"))){
				DOM_MeanSrcSigs2010 <- Srcs_dX_Ward
				DOM_VarSrcSigs2010 <- Srcs_dX_Var_Ward
			}
			if(YearMix==2012 & all(SourceNames[[f_Src]]%in%c("All Terrestrial", "All Macrophytes", "All Phytoplankton","Periphyton"))){
				DOM_MeanSrcSigs2012 <- Srcs_dX_Ward
				DOM_VarSrcSigs2012 <- Srcs_dX_Var_Ward
			}

			Temp_BugOut <- paste("bugsOut_", Cons[g_Cons], "_SrcComb", f_Src, "_","Pooled", sep="")
			Cons_Data <- subset(Data, Taxon==Cons[g_Cons] & Year==YearMix, select=c("Trophic","d13C","d15N","dD"))
			Cons_dX_Obs <- matrix(data=c(Cons_Data[,2], Cons_Data[,3], Cons_Data[,4]), ncol=3)
			ConsName <- as.character(subset(Data, Taxon==Cons[g_Cons] & Year==YearMix, select=Type)[1,1])#There should be a better way to do this...
			assign(Temp_BugOut, ConsMix(Cons_dX_Obs=Cons_dX_Obs, TL=TL[g_Cons], Srcs_dX=Srcs_dX_Ward, Srcs_dX_Var=Srcs_dX_Var_Ward, Water_dD_Mu, Water_dD_Var, FractModel=TRUE, SrcNames=SourceNames[[f_Src]], ConsName=ConsName, GraphTitle=GraphTitle[g_Cons], nChains=8, ChainLength=Iterations, Plot=FALSE))
			
			
			
			#Begin working on changing this code so that I can be more flexible with the number of sources, which will also require referring to columns by names (so that graphing labels can be proper, etc)
			# TempOutput <- get(Temp_BugOut)$sims.matrix[,1:length(SourceNames[[f_Src]])] #remove the deviance column, and rename the output matrix
			# 		colnames(TempOutput) <- SourceNames[[f_Src]] #name the columns of the bugs() output matrix according to the names of the sources
			TempOutput <- get(Temp_BugOut)$BUGSoutput$sims.matrix #remove the deviance column, and rename the output matrix
			colnames(TempOutput) <- c(SourceNames[[f_Src]], "Deviance") #name the columns of the bugs() output matrix according to the names of the sources
			NAsources <- names(SourceOpts)[!is.element(names(SourceOpts), SourceNames[[f_Src]])] #the names of the sources that weren't used in this run
			NA_srcs_df <- matrix(data=rep(NA, length(NAsources)), ncol=length(NAsources), dimnames=list(NULL, NAsources))
			
			if(g_Cons==1 & f_Src==1 & YearMix==2010){
				# ResourceUse <- data.frame("Year"=YearMix, "Month"="Pooled", "Consumer"=Cons[g_Cons], "Grouping"=f_Src, get(Temp_BugOut)$sims.matrix[,1:4])
				ResourceUse <- data.frame("Year"=YearMix, "Month"="Pooled", "Consumer"=Cons[g_Cons], "Grouping"=f_Src, TempOutput, NA_srcs_df) #will have to use make.names() to index by name
			}else{
				TempoResourceUse <- data.frame("Year"=YearMix, "Month"="Pooled", "Consumer"=Cons[g_Cons], "Grouping"=f_Src, TempOutput, NA_srcs_df)
				ResourceUse <- rbind(ResourceUse, TempoResourceUse)
			}
			# TempoResourceUse <- data.frame("Year"=YearMix, "Month"="Pooled", "Consumer"=Cons[g_Cons], "Grouping"=f_Src, get(Temp_BugOut)$sims.matrix[,1:4])
			# ResourceUse <- rbind(ResourceUse, TempoResourceUse)
		}#Finish loop the loop that handles the two sets of sources for this particular consumer
	}#Finish loop that estimates resource use for each consumer under 2 scenarios of available resources/ grouping of resources
	

	#Plot the composition of POM
	LegendTitle <- list(c("A)", "B)", "C)", "D)"), c("E", "F", "G", "H")) #CHANGED added )'s
	PubCex=1
	PanelNameAdj <- c(0.25, 0.33, 0.55, 0.58)
	# setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
	
	#Plot EPILIMNION
	pdf(file=paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/EpiPhyto_Post_", YearMix, ".pdf", sep=""), width=3.5, height=3.5, family="Times", pointsize=9)
	par(mfrow=c(2,2), las=1, mar=c(3,2.5,0.1,1), oma=c(0,0,0.2,0), cex=PubCex)
	
	TerrYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[1]"], from=0, to=1)$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[1]"], from=0, to=1),xlab="", ylab="", main="", bty="l", xaxt="s", zero.line=FALSE, ylim=TerrYLim)
	title(main=LegendTitle[[1]][1], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex) #CHANGED changed the adj from 1 to 0.1, added font.main=1, line from -0.5 to -1
	mtext("Terrestrial", side=3, line=-0.9, outer=FALSE, las=0, font=1, adj=PanelNameAdj[1], cex=PubCex) #CHANGED line from 0 to -1, deleted cex=0.85, changed font=3 to 1
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$f[1]*100, 0), "%", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex) #CHANGED deleted cex.main=0.85,
	
	PdCYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Epi"])$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Epi"]), main="", ylab="", xlab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PdCYLim)
	title(main=LegendTitle[[1]][3], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext(expression(Phytoplankton~phantom()^13*C), side=3, line=-1.1, outer=FALSE, las=0, font=1, adj=PanelNameAdj[3], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$P_dC_Epi, 1), "", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	
	PhytYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[2]"], from=0, to=1)$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[2]"], from=0, to=1),  main="", xlab="", ylab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PhytYLim)
	title(main=LegendTitle[[1]][2], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext("Phytoplankton", side=3, line=-0.9, outer=FALSE, las=0, font=1, adj=PanelNameAdj[2], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$f[2]*100, 0), "%", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	mtext("Fraction of POM", side=1, line=2, cex=PubCex, font=1, outer=FALSE)
	
	PdNYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Epi"])$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Epi"]),  main="", ylab="", xlab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PdNYLim)
	title(main=LegendTitle[[1]][4], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext(expression(Phytoplankton~phantom()^15*N), side=3, line=-1.1, outer=FALSE, las=0, font=1, adj=PanelNameAdj[4], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$P_dN_Epi, 2), "", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	mtext("Isotopic signature", side=1, line=2, cex=PubCex, font=1, outer=FALSE)
	
	mtext("Density", side=2, line=-1, font=1, las=0, outer=TRUE, cex=PubCex)
	dev.off()
	
	#Plot METALIMNION
	pdf(file=paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/MetaPhyto_Post_", YearMix, ".pdf", sep=""), width=3.5, height=3.5, family="Times", pointsize=9)
	par(mfrow=c(2,2), las=1, mar=c(3,2.5,0.1,1), oma=c(0,0,0.2,0), cex=PubCex)
	
	TerrYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[3]"], from=0, to=1)$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[3]"], from=0, to=1),xlab="", ylab="", main="", bty="l", xaxt="s", zero.line=FALSE, ylim=TerrYLim)
	title(main=LegendTitle[[1]][1], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex) #CHANGED changed the adj from 1 to 0.1, added font.main=1, line from -0.5 to -1
	mtext("Terrestrial", side=3, line=-0.9, outer=FALSE, las=0, font=1, adj=PanelNameAdj[1], cex=PubCex) #CHANGED line from 0 to -1, deleted cex=0.85, changed font=3 to 1
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$f[3]*100, 0), "%", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex) #CHANGED deleted cex.main=0.85,
	
	PdCYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Meta"])$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dC_Meta"]), main="", ylab="", xlab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PdCYLim, xlim=c(-75 , 0))
	title(main=LegendTitle[[1]][3], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext(expression(Phytoplankton~phantom()^13*C), side=3, line=-1.1, outer=FALSE, las=0, font=1, adj=PanelNameAdj[3], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$P_dC_Meta, 1), "", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	
	PhytYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[4]"], from=0, to=1)$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"f[4]"], from=0, to=1),  main="", xlab="", ylab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PhytYLim)
	title(main=LegendTitle[[1]][2], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext("Phytoplankton", side=3, line=-0.9, outer=FALSE, las=0, font=1, adj=PanelNameAdj[2], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean[[1]][4]*100, 0), "%", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	mtext("Fraction of POM", side=1, line=2, cex=PubCex, font=1, outer=FALSE)
	
	PdNYLim <- range(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Meta"])$y)*c(1, 1.15)
	plot.density(density(bugsOut_pt1$BUGSoutput$sims.matrix[,"P_dN_Meta"]),  main="", ylab="", xlab="", bty="l", xaxt="s", zero.line=FALSE, ylim=PdNYLim)
	title(main=LegendTitle[[1]][4], adj=0.025, line=-0.7, font.main=1, cex.main=PubCex)
	mtext(expression(Phytoplankton~phantom()^15*N), side=3, line=-1.1, outer=FALSE, las=0, font=1, adj=PanelNameAdj[4], cex=PubCex)
	title(paste(round(bugsOut_pt1$BUGSoutput$mean$P_dN_Meta, 2), "", sep=""),  adj=0.1, line=-1.75, font.main=3, cex.main=PubCex)
	mtext("Isotopic signature", side=1, line=2, cex=PubCex, font=1, outer=FALSE)
	
	mtext("Density", side=2, line=-1, font=1, las=0, outer=TRUE, cex=PubCex)
	
	# setwd("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2010Analysis/Figures_v8.3")
	# dev2bitmap(file="MetaPhyto_Post_v8.3.tif", type="tiffgray",height=3.5, width=3.5, res=200, font="Times", method="pdf", pointsize=12)
	dev.off()
	
	
}#End Year loop

# setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
# GroupChoose <- 2
# 
# pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/ConsSummary_ModelSelction.pdf", height=7, width=8.5, pointsize=8, family="Times")
# # dev.new(width=8, height=8)
# par(mfcol=c(3,3), mar=c(2.5,3.5,1,0.5), cex=1)
# ConsChoicesShort <- c("All Terrestrial"="Terr", "Epi. Phytoplankton"= "Epi Phyt", "Meta. Phytoplankton"="Meta Phyt", "DOM"="DOM", "All Macrophytes"="Macroph", "Periphyton"="Periphy", "Floating Macrophytes"="Float Mac", "Submersed Macrophytes"="Sub Mac", "All Phytoplankton"="Phyto", "Local Terrestrial"="Local Terr")
# ConsNameMedium <- c("Calanoid"="Calanoid", "Chaoborus"="Chaoborus", "Helisoma trivolvis"="Helisoma trivolvis", "FHM"="Fathead", "DAC"="Dace", "BHD1"="Sm Bullhead", "BHD2"="Lg Bullhead", "CMM"= "Mud Minn.", "PKS"="Pumpkinseed", "YWP"="Perch", "Mesocyclops"="Mesocyclops")
# for(i in 1:length(Cons)){
# 	TheseGroupings <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Month=="Pooled"))#[AllNamesThisuse]
# 	ScoreGroupings <- aggregate(TheseGroupings[,"Deviance"], by=list("Year"=TheseGroupings[,"Year"], "Grouping"=TheseGroupings[,"Grouping"]), FUN=mean)
# 	
# 	Scores2010 <- subset(ScoreGroupings, Year==2010)
# 	zScores2010 <- Scores2010
# 	zScores2010[,"x"] <- (Scores2010[,"x"]-mean(Scores2010[,"x"]))/sd(Scores2010[,"x"])
# 	
# 	Scores2012 <- subset(ScoreGroupings, Year==2012)
# 	zScores2012 <- Scores2012
# 	zScores2012[,"x"] <- (Scores2012[,"x"]-mean(Scores2012[,"x"]))/sd(Scores2012[,"x"])
# 	
# 	zScores <- data.frame("Grouping"=Scores2010[,"Grouping"], "zDev"=(zScores2010[,"x"]+zScores2012[,"x"])/2)
# 	GroupChoose <- zScores[which.min(zScores[,"zDev"]), "Grouping"]
# 	# GroupChoose <- ScoreGroupings[which.min(ScoreGroupings[,"x"]), "Grouping"]
# 	
# 	NamesThisUse0 <- ConsChoices[[i]][[GroupChoose]][!is.na(ConsChoices[[i]][[GroupChoose]])]
# 	NamesThisUse <- make.names(NamesThisUse0)
# 	# NamesThisUse <- make.names(ConsChoices[[i]][[GroupChoose]])
# 	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
# 	
# 	ThisRU0 <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
# 	ThisRU00 <- reshape(ThisRU0, varying=list(c(NamesThisUse)), times=NamesThisUse0, ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")
# 	ThisRU <- ThisRU00[,c("Year", "Month", "Consumer", "Grouping", "Source", "Proportion")]
# 	row.names(ThisRU) <- NULL
# 	ResourceNames <- ConsChoicesShort[as.character(NamesThisUse0)]
# 	
# 	ThisRU[,"Source"] <- factor(ThisRU[,"Source"], levels=ConsChoices[[i]][[GroupChoose]], ordered=TRUE)
# 	
# 	# aggregate(ThisRU[,"Proportion"], by=list("Year"=ThisRU[,1], "Source"=ThisRU[,"Source"]), FUN=median)
# 	
# 	ThisAt_Axis <- c(0.5,3,5.5,8)[1:length(ResourceNames)]
# 	boxplot(Proportion~Year+Source, data=ThisRU, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5)
# 	axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
# 	mtext(ConsNameMedium[Cons[i]], side=2, line=2)
# }
# dev.off()


# setwd("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
# source("DOM_Composition_v0.2.0.R")

NeatConsNames <- c("Calanoid"="S. oregonensis", "Mesocyclops"="Mesocyclops", "Chaoborus"="Chaoborus spp.", "Helisoma trivolvis"="H. trivolvis", "FHM"="P. promelas", "DAC"="Phoxinus spp.", "BHD1"="A. melas", "BHD2"="A. melas", "CMM"= "U. limi", "PKS"="L. gibbosus", "YWP"="P. flavescens")


save(list=c("Data", "EndMembers"), file=paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Data+Phyto.RData" ,sep=""))

save(list=c("ResourceUse", "ConsChoices"), file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ResourceUse2010&2012_ModelSelection.RData")

save.image(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012_ModelSelection.RData")





