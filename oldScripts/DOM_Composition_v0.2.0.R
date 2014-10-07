#Ryan Batt
#12-Jan-2012
#Figure out the composition of DOM in Ward Lake in 2010
#Starting at _v8 to mimic the current version of the IsoWard analysis
#DOM_Composition_v0.2.0 Repurposing the old DOM composition script for use with the 2010+2012 analysis.  Starting at this version when Cons_Mixture_Ward2010&2012_v0.2.2R is the latest version of the script.
	#Intended to be source()'d within the parent script.

OrigWD <- getwd()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2010Analysis")
source("ConsMix_v6.DOM.Median.R")

#Select the top 2 if on Snow Leopard, the bottom 2 if on Leopard, and the selection doesn't matter if on a PC
# WINE="/Applications/Darwine/Wine.bundle/Contents/bin/wine"
# WINEPATH="/Applications/Darwine/Wine.bundle/Contents/bin/winepath"
# WINEPATH="/opt/local/bin/winepath"
# WINE="/opt/local/bin/wine"

DOM2010 <- as.matrix(subset(Data, Taxon=="DOM" & Year==2010, select=c("d13C", "d15N", "dD")))
dimnames(DOM2010) <- NULL
DOM2012 <- as.matrix(subset(Data, Taxon=="DOM" & Year==2012, select=c("d13C", "d15N", "dD")))
dimnames(DOM2012) <- NULL

domOut2010 <- ConsMix(Cons_dX_Obs=DOM2010, TL=0, Srcs_dX=DOM_MeanSrcSigs2010, Srcs_dX_Var=DOM_VarSrcSigs2010, Water_dD_Mu=0, Water_dD_Var=0, FractModel=FALSE, SrcNames=c("All Terr.", "Macroph.", "Phytos.", "Periphyton"), ConsName=NULL, Omega_Info=c(0,0), TL_Var=0, Plot=TRUE, NewPlot=TRUE, DispMu=FALSE, GraphTitle=NULL, nChains=5, ChainLength=1000, debug=FALSE, WINE=WINE, WINEPATH= WINEPATH)
domOut2012 <- ConsMix(Cons_dX_Obs=DOM2012, TL=0, Srcs_dX=DOM_MeanSrcSigs2012, Srcs_dX_Var=DOM_VarSrcSigs2012, Water_dD_Mu=0, Water_dD_Var=0, FractModel=FALSE, SrcNames=c("All Terr.", "Macroph.", "Phytos.", "Periphyton"), ConsName=NULL, Omega_Info=c(0,0), TL_Var=0, Plot=TRUE, NewPlot=TRUE, DispMu=FALSE, GraphTitle=NULL, nChains=5, ChainLength=1000, debug=FALSE, WINE=WINE, WINEPATH= WINEPATH)

DOM_Comp00 <- rbind(data.frame("Year"=2010, domOut2010$sims.matrix[,1:4]), data.frame("Year"=2012, domOut2012$sims.matrix[,1:4]))
DOM_Comp0 <- reshape(DOM_Comp00, varying=list(c("DietF.1.", "DietF.2.", "DietF.3.", "DietF.4.")), times=1:4, ids=1:nrow(DOM_Comp00), timevar="Source", v.names="Proportion", direction="long")
DOM_Comp <- DOM_Comp0[,c("Year","Source", "Proportion")]
row.names(DOM_Comp) <- NULL
domResourceNames <- c("All Terr.", "Macroph.", "Phytos.", "Periphyton")

setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
pdf(file=paste(paste("DOM_Comp_", Version, sep=""), ".pdf", sep=""), height=3.5, width=3.5, pointsize=9)
	# dev.new(width=8, height=8)
par(mar=c(2.5,3.5,1,0.5), cex=1)

boxplot(Proportion~Year+Source, data=DOM_Comp, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(c(0.5,3,5.5,8),each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5, cex=1)
axis(side=1, at=c(0.5,3,5.5,8), labels=domResourceNames, cex.axis=1)
mtext("DOM", side=2, line=2, cex=1)

dev.off()

setwd(OrigWD)
