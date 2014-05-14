
#*****************************************************
#Begin for consumers and their respective sources
#*****************************************************

if(YearMix==2010){
	Cons <- c("Calanoid", "Chaoborus", "Helisoma trivolvis") #, "PKS", "FHM", "YWP", "CMM", "BHD", "Mesocyclops", "DAC")
	TL <- c(1, 2, 1) #, 2, 2, 3, 2.5, 2.5, 1.5, 1)
	GraphTitle <- c("Skistodiaptomus oregonensis", "Chaoborus spp.", "Helisoma trivolvis") #, "Lepomis gibbosus", "Pimephales promelas", "Perca flavescens", "Umbra limi", "Ameiurus melas",  "Mesocyclops spp.", "Phoxinus spp.")
}else{
	Cons <- c("Calanoid", "Chaoborus", "Helisoma trivolvis") #, "PKS", "FHM", "CMM", "BHD", "Mesocyclops", "DAC")
	TL <- c(1, 2, 1) #, 2, 2, 2.5, 2.5, 1.5, 1)
	GraphTitle <- c("Skistodiaptomus oregonensis", "Chaoborus spp.", "Helisoma trivolvis") #, "Lepomis gibbosus", "Pimephales promelas", "Umbra limi", "Ameiurus melas",  "Mesocyclops spp.", "Phoxinus spp.")

}
AllMacs <- c("Brasenia schreberi", "Chara", "Najas flexilis", "Nuphar variegata", "Nymphaea odorata", "Potamogeton amplifolius", "Potamogeton nodosus", "Potamogeton pusillus")
FloatMacs <- c("Brasenia schreberi", "Nuphar variegata", "Nymphaea odorata", "Potamogeton nodosus")
SubMacs <- c("Chara", "Najas flexilis", "Potamogeton amplifolius", "Potamogeton pusillus")
AllTerr <- c("Alder", "Sedge", "Tamarack", "Tree")
LocalTerr <- c("Alder", "Sedge", "Tamarack")


SourceOpts <- list("All Macrophytes"=AllMacs, "Floating Macrophytes"=FloatMacs, "Submersed Macrophytes"=SubMacs, "All Terrestrial"=AllTerr, "Local Terrestrial"=LocalTerr, "All Phytoplankton"=c("EpiPhyto", "MetaPhyto"), "Epi. Phytoplankton"="EpiPhyto", "Meta. Phytoplankton"="MetaPhyto", "DOM"="DOM", "Periphyton"="Periphyton")

ConsChoices <- list(
				"Calanoid"=list(c("All Terrestrial", "Epi. Phytoplankton", "Periphyton", "DOM") , c("Local Terrestrial", "Epi. Phytoplankton", "Meta. Phytoplankton", "Periphyton")),
				"Chaoborus"=list(c("All Terrestrial", "Epi. Phytoplankton", "Meta. Phytoplankton", "DOM") , c("Local Terrestrial", "Epi. Phytoplankton", "Meta. Phytoplankton", "DOM")),
				"Helisoma trivolvis"=list(c("All Terrestrial", "All Macrophytes", "Periphyton", "DOM"), c("Local Terrestrial", "Floating Macrophytes", "Submersed Macrophytes", "Periphyton"))
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
				
				
for(g_Cons in 1:length(Cons)){
	# SourceNames <- list(c("Terrestrial", "Pelagic Epilimnion", "Pelagic Metalimnion", "Benthic/ Littoral"), c("Terrestrial", "Pelagic", "Benthic", "Littoral"))
	# FirstSources <- list(TSources, TSources)
	# SecondSources <- list("EpiPhyto", c("EpiPhyto", "MetaPhyto"))
	# ThirdSources <- list("MetaPhyto", "Periphyton")
	# FourthSources <- list(c("Nuphar variegata", "Nymphaea odorata", "Brasenia schreberi", "Periphyton"), c("Nuphar variegata", "Nymphaea odorata", "Brasenia schreberi"))
	
	TempoCons <- Cons[g_Cons]
	SourceNames <- ConsChoices[[TempoCons]]
	FirstSources <- list(SourceOpts[[SourceNames[[1]][1]]], SourceOpts[[SourceNames[[2]][1]]])
	SecondSources <- list(SourceOpts[[SourceNames[[1]][2]]], SourceOpts[[SourceNames[[2]][2]]])
	ThirdSources <- list(SourceOpts[[SourceNames[[1]][3]]], SourceOpts[[SourceNames[[2]][3]]])
	FourthSources <- list(SourceOpts[[SourceNames[[1]][4]]], SourceOpts[[SourceNames[[2]][4]]])

	for(f_Src in 1:2){
		Source1 <- FirstSources[[f_Src]]
		Source2 <- SecondSources[[f_Src]]
		Source3 <- ThirdSources[[f_Src]]
		Source4 <- FourthSources[[f_Src]]
		# SourceData <- subset(Data, Trophic==0 & Taxon!="POM" & Year==YearMix | is.element(Type, c("Macrophyte", "Terrestrial")))
		# SourceTaxa <- as.character(unique(SourceData[,"Taxon"]))
		# Source_Means <- matrix(ncol=3, nrow=length(SourceTaxa), dimnames=list(SourceTaxa,NULL))
		# Source_Vars <- matrix(ncol=3, nrow=length(SourceTaxa), dimnames=list(SourceTaxa,NULL))
		# for(i in 1:length(SourceTaxa)){
			# Source_Means[i,] <- apply(subset(SourceData, Taxon==SourceTaxa[i], select=c("d13C","d15N","dD")), 2, mean)
			# Source_Vars[i,] <- apply(subset(SourceData, Taxon==SourceTaxa[i], select=c("d13C","d15N","dD")), 2, var)
		# }
		# Source_Means <- rbind(Source_Means, "EpiPhyto"=P_dX_Epi, "MetaPhyto"=P_dX_Meta)
		# Source_Vars <- rbind(Source_Vars, "EpiPhyto"=P_dX_Epi_Var, "MetaPhyto"=P_dX_Meta_Var)

		nSrcs <- length(SourceNames[[f_Src]])
		# SourceData_dX_Obs <- SourceData[,c("Taxon","d13C","d15N","dD")]
		# SourceData_dX_Obs <- rbind(SourceData_dX_Obs, Sim_P_dX_Epi_Obs, Sim_P_dX_Meta_Obs)
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
		
		Temp_BugOut <- paste("bugsOut", paste(Cons[g_Cons], paste("SrcComb", f_Src, sep=""), sep="_"), sep="_")
		Cons_Data <- subset(Data, Taxon==Cons[g_Cons] & Year==YearMix, select=c("Trophic","d13C","d15N","dD"))
		Cons_dX_Obs <- matrix(data=c(Cons_Data[,2], Cons_Data[,3], Cons_Data[,4]), ncol=3)
		ConsName <- as.character(subset(Data, Taxon==Cons[g_Cons] & Year==YearMix, select=Type)[1,1])#There should be a better way to do this...
		assign(Temp_BugOut, ConsMix(Cons_dX_Obs=Cons_dX_Obs, TL=TL[g_Cons], Srcs_dX=Srcs_dX_Ward, Srcs_dX_Var=Srcs_dX_Var_Ward, Water_dD_Mu, Water_dD_Var, FractModel=TRUE, SrcNames=SourceNames[[f_Src]], ConsName=ConsName, GraphTitle=GraphTitle[g_Cons], WINE=WINE, WINEPATH= WINEPATH, nChains=8, ChainLength=Iterations, Plot=FALSE))
	}#Finish loop the loop that handles the two sets of sources for this particular consumer
}#Finish loop that estimates resource use for each consumer under 2 scenarios of available resources/ grouping of resources