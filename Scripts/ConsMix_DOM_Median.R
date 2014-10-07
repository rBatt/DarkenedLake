#ConsMix_v3.R (old version, version log below)
#Ryan Batt (rbatt@wisc.edu)
#28-May-2011

#***NOTE***
#You do not need a .bugs file with this function.  The only dependencies are the input arguments (described below), having the package "R2WinBUGS" installed (if it isn't installed yet, type >install.packages(R2WinBUGS)), and having WinBUGS1.4.3 properly installed.  Will work on a PC or Mac.
#Use WinBUGS1.4.3.  
#Go to <mywebspace.wisc.edu/rbatt/web/ConsMixInstructions.txt> for details on how and where to install WinBUGS on a Mac or PC so that it will work with this function.  
#Go to <mywebspace.wisc.edu/rbatt/web/WinBugsR%20Instructions.rtf> for a detailed account as to how I installed WinBUGS and got it working with R and this function on a Mac (some of the steps are relevant to PC users as well).

#Happy Mixing :)

#Description:
#Calculate a consumers diet given observations of individual consumer 13C, 15N, and 2H, the mean and variance of the same signatures for 2 to 4 sources, and the mean and variance of water deuterium signature.

#Arguments:
#Cons_dX_Obs--	A data frame or matrix where each row is a different individual of the same consumer type, the first column is the d13C, the second column d15N, and the third column d2H.
#TL--			The trophic level of the consumer given as a scalar, where TL=1 is a primary consumer, TL=2 is a secondary consumer, etc.
#Srcs_dX--		A data frame or matrix where each column is a different source, the first row is the mean d13C, the second row the mean d15N, and the third row the mean d2H.
#Srcs_dX_Var--	Same as Srcs_dX, but replace means with variances.
#Water_dD_Mu--	The mean of the water d2H values.
#Water_dD_Var--	The variance of the water d2H values.
#SrcsNames--	This is used simply to label the panels of the final plot produced.  If left blank, the planels will be labeled, "Source1", etc.
#ConsName--		Serves to extract the Omega values (mean and variance) for a given consumer type.  Values are taken from Table 1 of Solomon et al. (2009), and available options include (c(mean, variance)):
					# "Snail"=c(0,0)
					# "Zooplankton"=c(0.20, 0.016)
					# "Chaoborus"=c(0.14, 0.036)
					# "Fish"=c(0.12, 0.04)
#Omega_Info--	A vector of the format c(mean, variance) to describe the proportion of tissue H derived from water (omega).  Needs to be supplied if ConsName is NULL, and will be ignored otherwise.
#TL_Var--		The variance of the supplied estimate of the consumer's trophic level.  The default is 0.1.
#Plot--			Display a density plot of the consumer's percent reliance on each source?
#DispMu--		Display the mean % reliance in each panel of the plot?
#GraphTitle--	The main title of the graph; defaults to ConsName if NULL (the default).
#nChains--		The number of Markov Chains
#ChainLength--	The number of iterations/ length of each Markov Chain
#...--			Arguments to be passed on to plot(); Note: This will throw an error if you specify an argument that I've already included (e.g., you cannot specify bty, because bty="l" already exists).

#Version3: Made changes so to be more compatible with Snow Leopard.  The WINE and WINEPATH defaults lend themselves to Leopoard (OSX 10.5), but for Snow Leopard (OSX 10.6) one should use
	# WINE <- "/Applications/Darwine/Wine.bundle/Contents/bin/wine"
	# WINEPATH <-"/Applications/Darwine/Wine.bundle/Contents/bin/winepath"
#Version4: Same as version 3 (different name to simply stay organized with main wrapper that uses this function); a few minor arguements added from Version3 (added after original V4)
#Version5: 
	#Cross-validated bandwidths used in plots
#Version6:
	#Added a new function that does not do any of the dietary water or fractionation (created as an easy way to analyze DOM)
	#Found an error where the number of consumers in the ConsMix function was calculted as the number of columns, but it should have been the number of rows


ConsMix <- function(Cons_dX_Obs, TL, Srcs_dX, Srcs_dX_Var, Water_dD_Mu, Water_dD_Var, FractModel=TRUE, SrcNames=NULL, ConsName=NULL, Omega_Info=NULL, TL_Var=0.1, Plot=TRUE, NewPlot=TRUE, DispMu=TRUE, GraphTitle=NULL, nChains=5, ChainLength=1000, WINE="/opt/local/bin/wine", WINEPATH="/opt/local/bin/winepath", ...){
require(R2WinBUGS)
nCons <- nrow(Cons_dX_Obs)#correction made from _v5: changed ncol to nrow
nSrcs <- ncol(Srcs_dX)

if(TL>1){#for consumers higher than primary consumers
	Cons_Nrich <- (TL-1)*3.4 + 2.52 #carnivore + herbivore enrichment #changed 2.5 to 2.52 as per vander zanden and rasmussen 2001
	Cons_Nrich_Var <- 3.4^2*(TL_Var) + (TL-1)^2*(0.16) + 6.25 #1.58 #changed the variance on 17-June-2011-- I have no idea where I got theold number from, but these N fracionations are from vander zanden and rasmussen 2001; #0.625 #variance of carniv. enrichment + herbiv. enrichment
	}else{#where the consumer is an herbivore
		Cons_Nrich <- 2.52#herbivore enrichment #changed 2.5 to 2.52 as per vander zanden and rasmussen 2001
		Cons_Nrich_Var <- 6.25 #Made the same change as above #0.625 #variance of herbivore enrichment
		}

#Define omega values
ConsOmegaInfo <- data.frame("Snail"=c(0,0), "Zooplankton"=c(0.20, 0.0016), "Chaoborus"=c(0.14, 0.0036), "Fish"=c(0.12, 0.0004)) #from Solomon paper, c(mean, Variance)
if(!is.null(ConsName)){OmegaInfo <- ConsOmegaInfo[,ConsName]}else{OmegaInfo <- Omega_Info}
Omega <- 1-(1-OmegaInfo[1])^TL #compound omega
Omega_Var <- TL*(1-OmegaInfo[1])^(TL-1)*OmegaInfo[2] - (1-OmegaInfo[1])^TL*log(1-OmegaInfo[1])*TL_Var #variance of compounded omega




#*******************************
#BEGIN writing the WinBUGS model
#*******************************
ConsMixModel <- function(){
	#Ryan Batt
	#Version 2: May 28, 2011
	#Version 2 is compatible with a variable number of sources (_v1 was set to 4)
	#Version 2 was also written to be shorter and possibly slightly faster than _v1
	#Portions based on code by CTS
	#Begin consumer calculations
	for (i in 1:nSrcs){
	  DietF.transform[i] ~ dunif(-3,5)
	  }

	#CLR math (untransform source proportions)
	for (i in 1:nSrcs){
  	expDietFTrans[i] <- exp(DietF.transform[i])
 	 }
	DietF.tot <- sum(expDietFTrans[])
	
	for (Src in 1:(nSrcs-1)){
		DietF[Src] <- exp(DietF.transform[Src])/DietF.tot
		}
	DietF[nSrcs] <- 1 - sum(DietF[1:(nSrcs-1)])
	
	for (i in 1:3){
	  residSd_Cons[i] ~ dunif(0,100)
	  residVar_Cons[i] <- residSd_Cons[i]*residSd_Cons[i]
	  }
	  
	#Calculate consumer means and variances
	#Calculate the signature of the consumer based on the signature of it's diet items and the relative proportion that those items are consumed
	for(i in 1:nSrcs){
		Cons_dX_Term[1,i] <- DietF[i]*Srcs_dX[1,i] #13C
		Cons_dX_Term[2,i] <- DietF[i]*Srcs_dX[2,i] #15N
		Cons_dX_Term[3,i] <- DietF[i]*Srcs_dX[3,i] #2H
		}
	Cons_dX[1] <- sum(Cons_dX_Term[1,]) #13C signature of consumer
	Cons_dX[2] <- sum(Cons_dX_Term[2,]) + Cons_Nrich #15N signature of consumer
	Cons_dX[3] <- sum(Cons_dX_Term[3,])*(1-Omega) + Omega*Water_dD_Mu #2H signature of consumer

	#Calculate the variance of the consumer's signatures
	for(i in 1:nSrcs){
		Cons_dC_Var_Term[i] <- pow(DietF[i],2)*Srcs_dX_Var[1,i] #13C variance contributed by each source
		Cons_dN_Var_Term[i] <- pow(DietF[i],2)*Srcs_dX_Var[2,i] #15N variance contributed by each source
		SourceVariance[i] <- pow(DietF[i],2)* Srcs_dX[3,i] #dD variance contributed by each source; dD variance is more complicated than the others
		#SourceMean[i] <- DietF[i]*Srcs_dx[3,i]# Same as  sum(Cons_dX_Term[3,])
		}
	Cons_dC_Var <- sum(Cons_dC_Var_Term[]) + residVar_Cons[1]
	Cons_dN_Var <- sum(Cons_dN_Var_Term[]) + Cons_Nrich_Var + residVar_Cons[2]
	Cons_dD_Var <- pow(Omega,2)*Water_dD_Var + pow(Water_dD_Mu,2)*Omega_Var + pow((1-Omega),2)*sum(SourceVariance[]) + pow(sum(Cons_dX_Term[3,]),2)*Omega_Var + residVar_Cons[3]
	# f = A*B; var(f)/f^2 = var(A)/A^2 + var(B)/B^2 <assuming cov(AB)==0>
	# var(f) = var(A)*B^2 + var(B)*A^2
	#Omega^2*WaterVariance + WaterMean^2*OmegaVariance +(1-Omega)^2*SourceVariance + SourceMean^2*OmegaVariance --- this is how the dHVar is propagated; BATT

	Cons_dX_Prec[1] <- (1/Cons_dC_Var)
	Cons_dX_Prec[2] <- (1/Cons_dN_Var)
	Cons_dX_Prec[3] <- (1/Cons_dD_Var)
	
	for (i in 1:nCons){
		for (iso in 1:3){
			Cons_dX_Obs[i,iso] ~ dnorm(Cons_dX[iso], Cons_dX_Prec[iso])
			}
		}
	}#End Model Function
#*****************************
#END writing the WinBUGS model
#*****************************



#*******************************
#BEGIN writing the WinBUGS model (without fractionation)
#*******************************
ConsMixModel_NoFract <- function(){
	#Ryan Batt
	#Version 2: May 28, 2011
	#Version 2 is compatible with a variable number of sources (_v1 was set to 4)
	#Version 2 was also written to be shorter and possibly slightly faster than _v1
	#Portions based on code by CTS
	#Begin consumer calculations
	#"DOM_BUGS_Model"~Ryan Batt, 10-Jan-2012
	#Remove all the omega, water, trohpic level, and fractionation stuf
	for (i in 1:nSrcs){
	  DietF.transform[i] ~ dunif(-3,5)
	  }

	#CLR math (untransform source proportions)
	for (i in 1:nSrcs){
  	expDietFTrans[i] <- exp(DietF.transform[i])
 	 }
	DietF.tot <- sum(expDietFTrans[])
	
	for (Src in 1:(nSrcs-1)){
		DietF[Src] <- exp(DietF.transform[Src])/DietF.tot
		}
	DietF[nSrcs] <- 1 - sum(DietF[1:(nSrcs-1)])
	
	for (i in 1:3){
	  residSd_Cons[i] ~ dunif(0,100)
	  residVar_Cons[i] <- residSd_Cons[i]*residSd_Cons[i]
	  }
	  
	#Calculate consumer means and variances
	#Calculate the signature of the consumer based on the signature of it's diet items and the relative proportion that those items are consumed
	for(i in 1:nSrcs){
		Cons_dX_Term[1,i] <- DietF[i]*Srcs_dX[1,i] #13C
		Cons_dX_Term[2,i] <- DietF[i]*Srcs_dX[2,i] #15N
		Cons_dX_Term[3,i] <- DietF[i]*Srcs_dX[3,i] #2H
		}
	Cons_dX[1] <- sum(Cons_dX_Term[1,]) #13C signature of consumer
	Cons_dX[2] <- sum(Cons_dX_Term[2,]) #15N signature of consumer
	Cons_dX[3] <- sum(Cons_dX_Term[3,]) #2H signature of consumer

	#Calculate the variance of the consumer's signatures
	for(i in 1:nSrcs){
		Cons_dC_Var_Term[i] <- pow(DietF[i],2)*Srcs_dX_Var[1,i] #13C variance contributed by each source
		Cons_dN_Var_Term[i] <- pow(DietF[i],2)*Srcs_dX_Var[2,i] #15N variance contributed by each source
		SourceVariance[i] <- pow(DietF[i],2)* Srcs_dX[3,i] #dD variance contributed by each source; dD variance is more complicated than the others
		#SourceMean[i] <- DietF[i]*Srcs_dx[3,i]# Same as  sum(Cons_dX_Term[3,])
		}
	Cons_dC_Var <- sum(Cons_dC_Var_Term[]) + residVar_Cons[1]
	Cons_dN_Var <- sum(Cons_dN_Var_Term[]) + residVar_Cons[2]
	Cons_dD_Var <- sum(SourceVariance[]) + residVar_Cons[3]
	# f = A*B; var(f)/f^2 = var(A)/A^2 + var(B)/B^2 <assuming cov(AB)==0>
	# var(f) = var(A)*B^2 + var(B)*A^2
	#Omega^2*WaterVariance + WaterMean^2*OmegaVariance +(1-Omega)^2*SourceVariance + SourceMean^2*OmegaVariance --- this is how the dHVar is propagated; BATT

	Cons_dX_Prec[1] <- (1/Cons_dC_Var)
	Cons_dX_Prec[2] <- (1/Cons_dN_Var)
	Cons_dX_Prec[3] <- (1/Cons_dD_Var)
	
	for (i in 1:nCons){
		for (iso in 1:3){
			Cons_dX_Obs[i,iso] ~ dnorm(Cons_dX[iso], Cons_dX_Prec[iso])
			}
		}
	}#End Model Function
#*****************************
#END writing the WinBUGS model (no fractionation)
#*****************************
ConsMixModel_File <- file.path(tempdir(), "ConsMixModel.bug")
write.model(if(FractModel==TRUE){ConsMixModel}else{ConsMixModel_NoFract}, ConsMixModel_File)
#file.show(ConsMixModel_File)


SupplyBUGS <- list(nSrcs, Srcs_dX, Srcs_dX_Var, nCons, Cons_dX_Obs, Cons_Nrich, Cons_Nrich_Var, Omega, Omega_Var, Water_dD_Mu, Water_dD_Var)
names(SupplyBUGS) <- strsplit(c("nSrcs, Srcs_dX, Srcs_dX_Var, nCons, Cons_dX_Obs, Cons_Nrich, Cons_Nrich_Var, Omega, Omega_Var, Water_dD_Mu, Water_dD_Var"), split=", ")[[1]]

SupplyBUGS_NoFract <- list(nSrcs, Srcs_dX, Srcs_dX_Var, nCons, Cons_dX_Obs)
names(SupplyBUGS_NoFract) <- strsplit(c("nSrcs, Srcs_dX, Srcs_dX_Var, nCons, Cons_dX_Obs"), split=", ")[[1]]

ParamBUGS <- c("DietF")

# if(.Platform$OS.type=="windows"){
	bugsOut <- jags(if(FractModel==TRUE){SupplyBUGS}else{SupplyBUGS_NoFract}, inits=NULL, ParamBUGS, ConsMixModel_File, n.chains=nChains, n.iter=ChainLength,...)
	# }else{
		# print(WINE)
		# 	print(WINEPATH)
		# bugsOut <- bugs(if(FractModel==TRUE){SupplyBUGS}else{SupplyBUGS_NoFract}, inits=NULL, ParamBUGS, ConsMixModel_File, n.chains=nChains, n.iter=ChainLength, program="winbugs", useWINE=TRUE, newWINE=TRUE, WINEPATH=WINEPATH, WINE=WINE,...)
		# }




if(Plot==TRUE){
	if(is.null(GraphTitle)){GraphTitle=ConsName}
	if(is.null(SrcNames)){SrcNames <- paste("Source",1:nSrcs,sep="")}
	#Plot the consumer diet
	if(nSrcs<=2){
		if(NewPlot==TRUE){dev.new(height=8, width=5)}
		par(mfrow=c(2,1), family="Times", las=1, mar=c(2,2,3,2), oma=c(3,3,3,0))
	}else{
		if(NewPlot==TRUE){dev.new(height=8, width=8)}
		par(mfrow=c(2,2), family="Times", las=1, mar=c(2,2,3,2), oma=c(3,3,3,0))
		}
	plot.density(density(bugsOut$sims.matrix[,1], from=0, to=1), xlab="", ylab="", main="", bty="l",...)
	mtext("Density", side=2, line=3, outer=FALSE, las=0)
	mtext(SrcNames[1], side=3, line=0.5, outer=FALSE, las=0, font=3)
	if(DispMu==TRUE){
		if(.Platform$OS.type=="windows"){
			legend("top", legend=paste("mu=", round(bugsOut$BUGSoutput$median[[1]][1], 2), sep=""), bty="n")
				}else{
					legend("top", legend=paste("µ=", round(bugsOut$BUGSoutput$median[[1]][1], 2), sep=""), bty="n")
					}	
		}
	
	if(nSrcs>=2){
		plot.density(density(bugsOut$BUGSoutput$sims.matrix[,2], from=0, to=1), main="", ylab="", xlab="", bty="l",...)
		if(nSrcs<=3){mtext("Proportion of Composition", side=1, line=3, outer=FALSE, las=0)}
		if(nSrcs==2){mtext("Density", side=2, line=3, outer=FALSE, las=0)}
		mtext(SrcNames[2], side=3, line=0.5, outer=FALSE, las=0, font=3)
			if(DispMu==TRUE){
			if(.Platform$OS.type=="windows"){
				legend("top", legend=paste("mu=", round(bugsOut$BUGSoutput$median[[1]][2], 2), sep=""), bty="n")
					}else{
						legend("top", legend=paste("µ=", round(bugsOut$BUGSoutput$median[[1]][2], 2), sep=""), bty="n")
						}	
			}
		}
	if(nSrcs>=3){
		plot.density(density(bugsOut$BUGSoutput$sims.matrix[,3], from=0, to=1), main="", xlab="Proportion of Composition", ylab="", bty="l",...)
		mtext("Proportion of Composition", side=1, line=3, outer=FALSE, las=0)
		mtext("Density", side=2, line=3, outer=FALSE, las=0)
		mtext(SrcNames[3], side=3, line=0.5, outer=FALSE, las=0, font=3)
			if(DispMu==TRUE){
			if(.Platform$OS.type=="windows"){
				legend("top", legend=paste("mu=", round(bugsOut$BUGSoutput$median[[1]][3], 2), sep=""), bty="n")
					}else{
						legend("top", legend=paste("µ=", round(bugsOut$BUGSoutput$median[[1]][3], 2), sep=""), bty="n")
						}	
			}
	}
	if(nSrcs==4){
		plot.density(density(bugsOut$BUGSoutput$sims.matrix[,4], from=0, to=1), main="", ylab="", xlab="Proportion of Composition", bty="l",...)
		mtext("Proportion of Composition", side=1, line=3, outer=FALSE, las=0)
		mtext(paste(SrcNames[4], collapse=", "), side=3, line=0.5, outer=FALSE, las=0, font=3)
			if(DispMu==TRUE){
			if(.Platform$OS.type=="windows"){
				legend("top", legend=paste("mu=", round(bugsOut$BUGSoutput$median[[1]][4], 2), sep=""), bty="n")
					}else{
						legend("top", legend=paste("µ=", round(bugsOut$BUGSoutput$median[[1]][4], 2), sep=""), bty="n")
						}	
			}
	}
	mtext(GraphTitle, side=3, line=0, font=2, cex=1.5, outer=TRUE)
	#*************************
	}

	
return(bugsOut)



}#End Function