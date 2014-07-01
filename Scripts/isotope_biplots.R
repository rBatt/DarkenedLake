
Save <- FALSE
library("cluster")
library("tikzDevice")

# ===================
# = Isotope biplots =
# ===================
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
TaxID <- as.numeric(BaseData[,"Taxon"])
BaseData <- cbind(BaseData, "TaxID"=TaxID)
XY_ErrorBars <- function(X, Y, sdX, sdY, Plot=FALSE){#This is a function that creates error bars given means and sd's
	#Horizontal Error Bars
	XnaughtH <- X-sdX
	XoneH <- X+sdX
	YnaughtH <- Y
	YoneH <- Y
	if(Plot==TRUE){segments(XnaughtH, YnaughtH, XoneH, YoneH, col="black")}
	#Vertical Error Bars
	XnaughtV <- X
	XoneV <- X
	YnaughtV <- Y-sdY
	YoneV <- Y+sdY
	if(Plot==TRUE){segments(XnaughtV,YnaughtV,XoneV,YoneV, col="black")}

	ErrorXlim <- c(min(XnaughtH), max(XoneH))
	ErrorYlim <- c(min(YnaughtV), max(YoneV))
	return(list("ErrorXlim"= ErrorXlim, "ErrorYlim"= ErrorYlim))
}#End function


encircle <- function(ind, cx="d13C", cy=NULL, ...){
	muMat <- as.matrix(mu1[ind,c(cx,cy)])
	sdMat <- as.matrix(sd1[ind,c(cx,cy)])/2
	# f1 <- rbind((muMat-sdMat), (muMat+sdMat))
	f01 <- cbind((muMat[,1]-sdMat[,1]),muMat[,2])
	f02 <- cbind((muMat[,1]+sdMat[,1]),muMat[,2])
	f03 <- cbind(muMat[,1], (muMat[,2]-sdMat[,2]))
	f04 <- cbind(muMat[,1], (muMat[,2]+sdMat[,2]))
	f1 <- rbind(f01, f02, f03, f04)
	# f1 <- muMat
	elo <- ellipsoidhull(f1)
	lines(predict(elo), ...)
	# polygon(predict(elo))
	# f2 <- rbind((muMat*1.1), (muMat-muMat*-0.1))
	# lines(predict(ellipsoidhull(f1)), ...)
	# lines(predict(ellipsoidhull(f2)), col="forestgreen", lty="dashed", lwd=2)
}


MetaPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="MetaPhyto"), "TaxID"=max(unique(TaxID))+2)
EpiPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="EpiPhyto"), "TaxID"=max(unique(TaxID))+1)
Sources_Cons <- rbind(BaseData[,c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")], MetaPhytos, EpiPhytos)
Data <- Sources_Cons
Grouping= c("Taxon", "Year")
IsoNames=c("d13C", "d15N", "dD")
nPlots <- combn(length(IsoNames),2)[,-3]
nPlots <- nPlots[,c(2,1)]
muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), mean)
muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))
sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), sd)
sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]])
Lims2 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]])


GraphLayers <- c("EndMembers", "Inverts2010&2012", "Fish2010&2012")
for(i in 1:length(GraphLayers)){

	if(GraphLayers[i]=="EndMembers"){
			Data <- droplevels(subset(BaseData, Trophic==0 & is.element(Type, EndMemberTaxa), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	if(GraphLayers[i]=="Inverts2010&2012"){
			Data <- droplevels(subset(BaseData, Trophic>=0 & is.element(Type, InvertTaxa), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	if(GraphLayers[i]=="Fish2010&2012"){
			Data <- droplevels(subset(BaseData, Trophic>=0 & is.element(Type, "Fish"), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	
	# Data <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")) & Taxon!="POM" & Taxon!="Nija" & Taxon!="YWP" & Taxon!="PKS" &Taxon!="BHD2" & Taxon!="Mesocyclops", select=c("Year","Trophic","Taxon","d13C","d15N","dD")))



	DataNew <- Data
	UniqueTaxa <- as.character(unique(Data[,"Taxon"]))
	RealTaxa <- as.character(UniqueTaxa) #c("Helisoma trivolvis", "DOM", "Perca flavescens", "Lepomis gibbosus", "Pimaphales promelas", "Periphyton", "Potamogeton pussillus", "Chara sp.", "Nymphaea odorata", "Carex sp.", "Alnus incana subsp. Rugosa", "Skistodiaptomus oregonensis", "Chaoborus spp.", "Nuphar lutea", "Brasenia schreberi", "Trees") #"DOM", 

	if(GraphLayers[i]=="EndMembers"){
		MetaPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="MetaPhyto"), "TaxID"=max(unique(TaxID))+2)
		EpiPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="EpiPhyto"), "TaxID"=max(unique(TaxID))+1)
		Sources_Cons <- rbind(DataNew, MetaPhytos, EpiPhytos)
	}else{
		Sources_Cons <- DataNew
	}
	
	#BiPlot <- function(Grouping="Taxon", IsoNames=c("d13C", "d15N", "dD"), Data=NULL, SDs=NULL, PlotLabels=TRUE, Legend=TRUE, LegendLoc="topleft", ...){
		#"Grouping" is the name of the column that contains the names of the different varieties of your group
		#"IsoNames" is the name of the columns that contain the isotope values
		#"Data" is the data frame that contains your means (and sd's if you have them)
		#"SD's" is the data frame that contains your standard deviations (see below for more info on making this)
		#"PlotLabels" asks if you want the name of each group variety plotted with it's data point.  If TRUE (the default), names will often overlap
		#"LegendLoc" describes the location of the legend in the event that you don't want the labels to be plotted.  Options include, "topleft"(the default), "bottomleft", "bottomright", "topright", "center", etc. ?legend for more options
		#"..." pass on other arguments to plot(); for example, type 'main="FUN!"' as an argument to the biplot function (i.e., BiPlot(main="FUN!")) to make "FUN!" the title of each plot.

	PubCex <- 1
	Grouping= c("Taxon", "Year")
	IsoNames=c("d13C", "d15N", "dD")
	SDs=NULL
	LegendLoc="topleft"
	Data=Sources_Cons
	PlotLabels=FALSE
	Legend=FALSE

	nPlots <- combn(length(IsoNames),2)[,-3]#Leave this line alone
	# LabelIso <- list(expression(phantom()^13*C), expression(phantom()^15*N), expression(phantom()^2*H))
	GroupNames <- paste("Group",1:length(Grouping),sep=".")#Leave this line alone

	muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), mean)
	muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))

	if(is.null(SDs)){
		sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), sd)
		sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
		}else{
			sdDataGroup <- SDs
	}

	nPlots <- nPlots[,c(2,1)]

	LegendTitle <- list(c("A", "B"), c("C", "D"), c("E", "F"))[[i]]
	if(Save & i==1){
		
		options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
		    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
		    "\\usepackage{amssymb}","\\usepackage{tensor}"))
		tikz(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/BiPlots/Ward2010&2012_Full_IsoBiPlots.tex", width=6.1, height=5.25, standAlone=TRUE, 
		packages = c("\\usepackage{tikz}",
		                 "\\usepackage[active,tightpage,psfixbb]{preview}",
		                 "\\PreviewEnvironment{pgfpicture}",
		                 "\\setlength\\PreviewBorder{0pt}",
		                 "\\usepackage{amssymb}",
						"\\usepackage{tensor}")
		)
	}else if(i==1){
		dev.new(width=6.1, height=4.25)
	}
	if(i==1){
		par(family="Times", las=0, mfcol=c(2,3), mar=c(1.5, 1.75, 0.1, 0.1), oma=c(0.5,0,0,0), cex=PubCex, ps=8, mgp=c(1.5, 0.15, 0), tcl=-0.15)
	}else{
		par(mar=c(1.5, 0.85, 0.1, 0.1))
	}
	
	# Create first plot
	plot(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", xaxt="s", cex.axis=1)
	# axis(side=1, labels=FALSE)
	if(i == 1){
		mtext("$\\delta^2\\mathrm{H}$", side=2, line=1.1, las=0, cex=PubCex)
	}
	
	
	# Add Circles
	if(GraphLayers[i]=="EndMembers"){
		mu1 <- muDataGroup #being lazy -- when looping through the groups, I need to use the muDataGroup (and sd) for the end members when drawing circles (called in function)
		sd1 <- sdDataGroup
		floatI <- muDataGroup[,"Taxon"] %in% c("Brasenia schreberi", "Nuphar variegata", "Nymphaea odorata", "Potamogeton nodosus")
		subI <- muDataGroup[,"Taxon"] %in% c("Chara","Najas flexilis", "Potamogeton amplifolius", "Potamogeton pusillus")
		terrI <- muDataGroup[,"Taxon"] %in% c("Alder", "Sedge", "Tamarack", "Tree")
		phytoI <- muDataGroup[,"Taxon"] %in% c("EpiPhyto", "MetaPhyto")
		periI <- muDataGroup[,"Taxon"] %in% c("Periphyton")
	}
	encircle(floatI, cy="dD", col="springgreen4", lty="dashed", lwd=2)
	encircle(subI, cy="dD", col="springgreen4", lty="dashed", lwd=2)
	encircle(terrI, cy="dD", col="burlywood4", lty="solid", lwd=2)
	encircle(phytoI, cy="dD", col="springgreen4", lty="solid", lwd=2)
	encircle(periI, cy="dD", col="springgreen4", lty="solid", lwd=2)
	
	#Add error bars
	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]], Plot=T)
	
	# Add points
	colGroups <- 1L + (muDataGroup[,"Year"]==2012)*2L + (muDataGroup[,"Year"]==2010 & (muDataGroup[,"Taxon"]%in%c("EpiPhyto", "MetaPhyto", "Periphyton")|(GraphLayers[i]!="EndMembers")))
	points(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=c(col2010, col2012)[colGroups], cex=2.7, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]])
	
	# Get Ready to add #'s inside points, need to adjust overlapping
	TaxID <- as.numeric(muDataGroup[,"TaxID"])
	TXTmuDataGroup <- muDataGroup
	
	if(GraphLayers[i]=="EndMembers"){
		Move3_Index <- TXTmuDataGroup[,"TaxID"]==3
		Move3_d_CND <- c(2.5, 0, 0)
		TXTmuDataGroup[Move3_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move3_Index,c("d13C","d15N","dD")] + Move3_d_CND
		
		# Move12_Index <- TXTmuDataGroup[,"TaxID"]==12
		# Move12_d_CND <- c(0,0,-0.0)
		# TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] + Move12_d_CND
		
		Move15_Index <- TXTmuDataGroup[,"TaxID"]==15
		Move15_d_CND <- c(0,0,10)
		TXTmuDataGroup[Move15_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move15_Index,c("d13C","d15N","dD")] + Move15_d_CND

		Move19_Index <- TXTmuDataGroup[,"TaxID"]==19
		Move19_d_CND <- c(0,-1.95,-0.0)
		TXTmuDataGroup[Move19_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move19_Index,c("d13C","d15N","dD")] + Move19_d_CND
		
		Move2112_Index <- TXTmuDataGroup[,"TaxID"]==21 & TXTmuDataGroup[,"Year"]==2012
		Move2112_d_CND <- c(2.9, -0.6, 6)
		TXTmuDataGroup[Move2112_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move2112_Index,c("d13C","d15N","dD")] + Move2112_d_CND
		
		txtGrps <- colGroups
		txtGrps[Move3_Index | Move15_Index | Move19_Index | Move2112_Index] <- 1
		txtGrps2 <- txtGrps
		
	}else if(GraphLayers[i]=="Fish2010&2012"){
		
		
		Move710_Index <- TXTmuDataGroup[,"TaxID"]==7 & TXTmuDataGroup[,"Year"]==2010
		Move710_d_CND <- c(2.1, 0.5, 5.0)
		TXTmuDataGroup[Move710_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move710_Index,c("d13C","d15N","dD")] + Move710_d_CND
		
		Move712_Index <- TXTmuDataGroup[,"TaxID"]==7 & TXTmuDataGroup[,"Year"]==2012
		Move712_d_CND <- c(-2, 0.2, 0.0)
		TXTmuDataGroup[Move712_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move712_Index,c("d13C","d15N","dD")] + Move712_d_CND
		
		Move910_Index <- TXTmuDataGroup[,"TaxID"]==9 & TXTmuDataGroup[,"Year"]==2010
		Move910_d_CND <- c(0.45, -1.75, -0.0)
		TXTmuDataGroup[Move910_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move910_Index,c("d13C","d15N","dD")] + Move910_d_CND
		
		txtGrps <- colGroups
		txtGrps[Move710_Index | Move712_Index] <- 1
		txtGrps2 <- txtGrps + Move910_Index*2
	}else if(GraphLayers[i] == "Inverts2010&2012"){
		
		txtGrps <- colGroups
		txtGrps2 <- txtGrps
		
	}
	

	
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,1]]], labels=as.character(TaxID), cex=PubCex, col=c("black", "white", "white")[txtGrps], font=2)
	title(main=LegendTitle[1], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)

	# Create 2nd plot
	plot(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims2[[1]], ylim=Lims2[[2]], bty="l", xaxt=ifelse(2==1, "n", "s"), cex.axis=1)
	if(i == 1){
		mtext("$\\delta^{15}\\mathrm{N}$", side=2, line=1, las=0, cex=PubCex)	
	}
	if(i ==2){
		mtext("$\\delta^{13}\\mathrm{C}$", side=1, line=1, outer=FALSE, cex=PubCex)
	}
	
	
	
	# Add circles
	encircle(floatI, cy="d15N", col="springgreen4", lty="dashed", lwd=2)
	encircle(subI, cy="d15N", col="springgreen4", lty="dashed", lwd=2)
	encircle(terrI, cy="d15N", col="burlywood4", lty="solid", lwd=2)
	encircle(phytoI, cy="d15N", col="springgreen4", lty="solid", lwd=2)
	encircle(periI, cy="d15N", col="springgreen4", lty="solid", lwd=2)
	
	# Add error bars
	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]], Plot=T)

	# Add points
	points(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=c(col2010, col2012)[colGroups], cex=2.7, xlab="", ylab="", xlim=Lims2[[1]], ylim=Lims2[[2]])
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,2]]], labels=as.character(TaxID), cex=PubCex, col=c("black", "white", "white")[txtGrps2], font=2)
	#}

	title(main=LegendTitle[2], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)
	
	
	if(Save){dev.off()}
}