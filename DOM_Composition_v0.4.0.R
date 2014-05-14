#Ryan Batt
#12-Jan-2012
#Figure out the composition of DOM in Ward Lake in 2010
#Starting at _v8 to mimic the current version of the IsoWard analysis
#DOM_Composition_v0.2.0 Repurposing the old DOM composition script for use with the 2010+2012 analysis.  Starting at this version when Cons_Mixture_Ward2010&2012_v0.2.2R is the latest version of the script.
	#Intended to be source()'d within the parent script.
#_v0.4.0 (31-Oct-2013) I've realized that the mass of aquashade per liter is about 1/12 of what I thought it was ... so I'm seeing how/ if this changes anything at all

OrigWD <- getwd()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2010Analysis")
source("ConsMix_v6.DOM.Median.R")


#Aquashade dN = -0.02,	dC= -26.82,	 dD = -65.4

#http://chem.sis.nlm.nih.gov/chemidplus/jsp/common/ChemInfo.jsp?calledFrom=lite&type=formulas
#Aquashade formula =   C37-H66-N2-O9-S3   C16-H12-N4-O9-S2   2(NH3)  3Na
#aCarbon = (53)*(12) #(atoms)*(g/mol) 
#aHydrogen = (84)*(1)#(atoms)*(g/mol) 
#aNitrogen = (8)*(14)#(atoms)*(g/mol) 
#aOxygen = (18)*(16)#(atoms)*(g/mol) 
#aSulfur = (5)*(32)#(atoms)*(g/mol) 
#aSodium = (3)*(23)#(atoms)*(g/mol) 
#aMW = 1349 #computed from the above
#My count of the first component is that it has 36 H, not 66; the second has 9; 3H on the two ammonia. My total count of the hydrogen is 51, which doesn't match this reference, but matches the one below. I think that this reference is wrong.

#http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=135020425&viewopt=PubChem
#Molecular Formula: C53H51N8Na3O18S5   
#Molecular Weight: 1317.309148   
#InChIKey: CNXKZKXMIIEMNW-UHFFFAOYSA-K
#ok, given this formula, carbon forms (612/1317.3)  of the molecular weight, or 46.46% (CPSIL calculated 49.31%)
	#nitrogen is (112/1317.3) of the molecular weight, or 8.5% (CPSIL calculated 4.39%)
	#hydrogen is (51/1317.3) of the molecular weight, or 3.9% (CPSIL calculated 4.2%)
#Therefore, 1.5 mg/L of aquashade, is (1.5 * 0.4646) = 0.67mg/L carbon
#In 2012, Ward had an average of 10.65 mg/L of DOC
#Aquashade formed (0.67/10.65) = 6.3% of the dissolved organic carbon mass in Ward
#For 2012, CPSIL is saying that the DOM had ~10% C weight, 0.4% N weight, and 2.3% H weight.
#For 2010, DOM had 2.3% C, 0.6% N (but that's only 1 sample; if I use the other two, the C is 3.15%, and N is 0.43%)


(0.1 - .46*(1.5/106.5))*100 #to get the percent DOM (sans aquashade) that is carbon
(0.004 - .085*(1.5/106.5))*100#to get the percent DOM (sans aquashade) that is nitrogen
(0.023 - .085*(1.5/106.5))*100#to get the percent DOM (sans aquashade) that is hydrogen


#Aquashade #this also needed to be fixed in 0.4.0 for what's need a few steps down
#1) 0.258825 mg		1.5 mg
#2) 46% C (0.12 mg)		(0.67 mg)
#3) 8.5% N (0.022 mg)		(0.13 mg)
#4) 3.9% H (0.01 mg)		(0.06 mg)


#DOM (total, including aquashade, taken from 2012 CPSIL data that came w/ isotopes)
#1) 106.5 mg (in 2012 there was 10.65mg C, and C was 10%, so 10.65/0.1=106.5)
#2) 10% C (10.65 mg)
#3) 0.4% N (0.426 mg)
#4) 2.3% H (2.415 mg)


#DOM (NOT including aquashade)
#1) 105 mg (in 2012 there was 10.65mg C, and C was 10%, and there was 0.258825 mg from Aquashade, so (10.65/0.1)-1.5=105)
#2) 9.35% C (9.96 mg)	#(0.1 - .46*(0.258825/106.5))*100 #to get the percent DOM (sans aquashade) that is carbon
#3) 0.28% N (0.2985 mg)	#(0.004 - .085*(0.258825/106.5))*100 #to get the percent DOM (sans aquashade) that is nitrogen
#4) 2.18% H (2.322 mg)	#(0.023 - .039*(0.258825/106.5))*100 #to get the percent DOM (sans aquashade) that is hydrogen


#(Aquashade/(DOM-Aquashade))
#1) 1.4%
#2) 6.73% C
#3) 43.55% N
#4) 2.58% H


#(Aquashade/(TotalDOM)) #this is what i'm fixing in 0.4.0 b/c it's what is used later ... other bits don't really matter
#total DOM is the 10.65 mg, calculated from the [DOC] and the %C in DOm samples sent in to CPSIL
#1) 0.243%		1.4%
#2) 1.12%	 (0.258825*0.46)/10.65	#6.29% C
#3) 5.164%	(0.258825*0.085)/0.426	30.52% N
#4) 0.418%	(0.258825*0.039)/2.415	2.48% H

#Aquashade dN = -0.02,	dC= -26.82,	 dD = -65.4


-1*0.0629*-26.82

#dDOM[total] = fNormal*dDOM[normal] + fAquashade*dDOM[aquashade]
#dDOM[total] - fAquashade*dDOM[aquashade] = fNormal*dDOM[normal]
#(dDOM[total] - fAquashade*dDOM[aquashade])/fNormal = dDOM[normal]

# CorrectDOM <- function(dTot, dAqua=c("dC"=-26.82, "dN"=-0.02, "dD"=-65.4), fAqua=c("fC"=0.0629, "fN"=0.3052, "fD"=0.0248)){
CorrectDOM <- function(dTot, dAqua=c("dC"=-26.82, "dN"=-0.02, "dD"=-65.4), fAqua=c("fC"=0.0111793, "fN"=0.05164349, "fD"=0.004179783)){
	fNorm <- 1- fAqua
	dNorm <- matrix(data=rep(NA, 3*nrow(dTot)), nrow=nrow(dTot))
	for(i in 1:nrow(dTot)){
		dNorm[i,] <- (dTot[i,]-fAqua*dAqua)/fNorm
	}
	return(dNorm)
}



#Select the top 2 if on Snow Leopard, the bottom 2 if on Leopard, and the selection doesn't matter if on a PC
# WINE="/Applications/Darwine/Wine.bundle/Contents/bin/wine"
# WINEPATH="/Applications/Darwine/Wine.bundle/Contents/bin/winepath"
# WINEPATH="/opt/local/bin/winepath"
# WINE="/opt/local/bin/wine"

DOM2010 <- as.matrix(subset(Data, Taxon=="DOM" & Year==2010, select=c("d13C", "d15N", "dD")))
dimnames(DOM2010) <- NULL
DOM2012_0 <- as.matrix(subset(Data, Taxon=="DOM" & Year==2012, select=c("d13C", "d15N", "dD")))
dimnames(DOM2012_0) <- NULL
DOM2012 <- CorrectDOM(dTot=DOM2012_0)

domOut2010 <- ConsMix(Cons_dX_Obs=DOM2010, TL=0, Srcs_dX=DOM_MeanSrcSigs2010, Srcs_dX_Var=DOM_VarSrcSigs2010, Water_dD_Mu=0, Water_dD_Var=0, FractModel=FALSE, SrcNames=c("All Terr.", "Macroph.", "Phytos.", "Periphyton"), ConsName=NULL, Omega_Info=c(0,0), TL_Var=0, Plot=FALSE, NewPlot=TRUE, DispMu=FALSE, GraphTitle=NULL, nChains=5, ChainLength=1000, debug=FALSE, WINE=WINE, WINEPATH= WINEPATH)
domOut2012 <- ConsMix(Cons_dX_Obs=DOM2012, TL=0, Srcs_dX=DOM_MeanSrcSigs2012, Srcs_dX_Var=DOM_VarSrcSigs2012, Water_dD_Mu=0, Water_dD_Var=0, FractModel=FALSE, SrcNames=c("All Terr.", "Macroph.", "Phytos.", "Periphyton"), ConsName=NULL, Omega_Info=c(0,0), TL_Var=0, Plot=FALSE, NewPlot=TRUE, DispMu=FALSE, GraphTitle=NULL, nChains=5, ChainLength=1000, debug=FALSE, WINE=WINE, WINEPATH= WINEPATH)

DOM_Comp00 <- rbind(data.frame("Year"=2010, domOut2010$sims.matrix[,1:4]), data.frame("Year"=2012, domOut2012$sims.matrix[,1:4]))
DOM_Comp0 <- reshape(DOM_Comp00, varying=list(c("DietF.1.", "DietF.2.", "DietF.3.", "DietF.4.")), times=1:4, ids=1:nrow(DOM_Comp00), timevar="Source", v.names="Proportion", direction="long")
DOM_Comp <- DOM_Comp0[,c("Year","Source", "Proportion")]
row.names(DOM_Comp) <- NULL
domResourceNames <- c("Terr", "Macroph", "Phyto", "Periphy")

setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
# pdf(file=paste(paste("DOM_Comp_", Version, sep=""), ".pdf", sep=""), height=3.5, width=3.5, pointsize=9)
png(file=paste(paste("DOM_Comp_", Version, sep=""), ".png", sep=""), units="in", res=300, height=3.5, width=3.5, pointsize=10)
	# dev.new(width=8, height=8)
par(mar=c(2.5,3.5,1,0.5), cex=1, ps=10)

boxplot(Proportion~Year+Source, data=DOM_Comp, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(c(0.5,3,5.5,8),each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5, cex=1)
axis(side=1, at=c(0.5,3,5.5,8), labels=domResourceNames, cex.axis=1)
mtext("DOM", side=2, line=2, cex=1)

dev.off()

setwd(OrigWD)
