
# ======================================
# = Load Data for sensitivity analysis =
# ======================================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012_ModelSelection.RData")
cn <- c("Calanoid", "Chaoborus", "Helisoma trivolvis", "FHM", "DAC", "CMM", "BHD1")
sensit <- list()
sensit2 <- list()
for(i in 1:length(cn)){
	crit <- ResourceUse[ResourceUse[,"Consumer"]==cn[i],]
	crit2 <- aggregate(crit[,c("All.Terrestrial","All.Phytoplankton", "All.Macrophytes", "Submersed.Macrophytes","Floating.Macrophytes","Deviance", "Periphyton")], by=crit[,c("Grouping","Year")], FUN=mean)
	if(cn[i]=="Helisoma trivolvis"){
		crit3 <- crit2[is.na(crit2[,"All.Phytoplankton"]),]
	}else{
		if(is.element(cn[i], c("Calanoid", "Chaoborus"))){
			crit3 <- crit2[!is.na(crit2[,"All.Phytoplankton"]),]
		}else{
			crit3 <- crit2
		}
	}
	
	crit4 <- rbind(colMeans(crit3[crit3[,"Year"]==2010,], na.rm=TRUE), colMeans(crit3[crit3[,"Year"]==2012,], na.rm=TRUE))
	critX0 <- crit2[!is.na(crit2[,"All.Phytoplankton"]),]
	critX <- rbind(colMeans(critX0[critX0[,"Year"]==2010,], na.rm=TRUE), colMeans(critX0[critX0[,"Year"]==2012,], na.rm=TRUE))
	sensit[[cn[i]]] <- crit4
	sensit2[[cn[i]]] <- critX
}
save(sensit, sensit2, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Orgd_ConsMix_Sensit.RData")