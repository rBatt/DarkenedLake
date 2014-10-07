#Ryan Batt
#21-May-2011
#Create Bi-Plots from isotope data
#Data can be entered in 2 formats, 1) Data that does not have mean and sd calculated, or 2)The argument "Data" is just the means and the groups, and the "SDs" argument is the standard deviations entered as a data.frame(). See examples below with "PI Data" ;)
 
rm(list=ls())
graphics.off()
Version <- "v0.4.2"
FigureFolder <- paste("Figures/Figures_", Version, sep="")
# YearMix <- 2012
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load(paste("Data+Phyto_", Version, ".RData",sep=""))
ThisVersion <- "v0.3.0"
CexScale <- 0.8


CombYears <- which(is.element(Data[,"Type"], c("Macrophyte", "Terrestrial"))) #I don't want to make separate points for samples types that I don't really expect to change between years, such as the macrophytes or the terrestrial samples.  I just want to assume that these were the same in the different years, while the consumers, DOM, POM etc. may have changed.
Data[CombYears,"Year"] <- 2010
ExcludeTaxa <- c("POM","Nija","YWP","PKS","BHD2","Mesocyclops", "Pickerell Weed")
ExcludeSamples <- c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")
EndMemberTaxa <- c("Periphyton","DOM","Macrophyte","Terrestrial")
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
	if(Plot==TRUE){segments(XnaughtH, YnaughtH, XoneH, YoneH, col="gray")}
	#Vertical Error Bars
	XnaughtV <- X
	XoneV <- X
	YnaughtV <- Y-sdY
	YoneV <- Y+sdY
	if(Plot==TRUE){segments(XnaughtV,YnaughtV,XoneV,YoneV, col="gray")}

	ErrorXlim <- c(min(XnaughtH), max(XoneH))
	ErrorYlim <- c(min(YnaughtV), max(YoneV))
	return(list("ErrorXlim"= ErrorXlim, "ErrorYlim"= ErrorYlim))
}#End function





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
# Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]])
Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]])
Lims2 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]])
# rm(list=c("MetaPhytos", "EpiPhytos","Sources_Cons", "Data", "Grouping", "IsoNames", "nPlots", "muDataGroup", "sdDataGroup"))

GraphLayers <- c("EndMembers", "Inverts2010", "Inverts2010&2012", "Fish2010&2012")
for(i in 1:length(GraphLayers)){

	# if(GraphLayers[i]=="EndMembers"){
	# 		Data <- droplevels(subset(Data, Trophic==0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon, ExcludeTaxa) & is.element(Type, EndMemberTaxa), select=c("Year","Trophic","Taxon","d13C","d15N","dD")))
	# }
	# if(GraphLayers[i]=="Inverts2010"){
	# 		Data <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon, ExcludeTaxa) & is.element(Type, InvertTaxa) & Year==2010, select=c("Year","Trophic","Taxon","d13C","d15N","dD")))
	# }
	# if(GraphLayers[i]=="Inverts2010&2012"){
	# 		Data <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon,ExcludeTaxa) & is.element(Type, InvertTaxa), select=c("Year","Trophic","Taxon","d13C","d15N","dD")))
	# }
	# if(GraphLayers[i]=="Fish2010&2012"){
	# 		Data <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon, ExcludeTaxa) & is.element(Type, "Fish"), select=c("Year","Trophic","Taxon","d13C","d15N","dD")))
	# }
	
	if(GraphLayers[i]=="EndMembers"){
			Data <- droplevels(subset(BaseData, Trophic==0 & is.element(Type, EndMemberTaxa), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	if(GraphLayers[i]=="Inverts2010"){
			Data <- droplevels(subset(BaseData, Trophic>=0 & is.element(Type, InvertTaxa) & Year==2010, select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
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
	# for(i in 1:length(UniqueTaxa)){
	
		# DataNew[,"Taxon"] <- gsub(UniqueTaxa[i], RealTaxa[i], as.character(DataNew[,"Taxon"]))
	
	# }

	#gsub(UniqueTaxa, RealTaxa, as.character(Data[,"Taxon"]))

	# MetaPhytos <- data.frame("Trophic"=rep(0, nrow(Sim_P_dX_Meta_Obs)),"Taxon"=rep("MetaPhyto", nrow(Sim_P_dX_Meta_Obs)),  "d13C"=Sim_P_dX_Meta_Obs[,2], "d15N"=Sim_P_dX_Meta_Obs[,3], "dD"=Sim_P_dX_Meta_Obs[,4], "Year"=Sim_P_dX_Meta_Obs[,"Year"])
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
	LabelIso <- list(expression(phantom()^13*C), expression(phantom()^15*N), expression(phantom()^2*H))
	GroupNames <- paste("Group",1:length(Grouping),sep=".")#Leave this line alone

	muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), mean)
	muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))
	# print(muDataGroup)

	if(is.null(SDs)){
		sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), sd)
		sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
		}else{
			sdDataGroup <- SDs
	}

	nPlots <- nPlots[,c(2,1)]

	LegendTitle <- c("A", "B")
	# dev.new(height=7.25, width=7.25) #dD vs. d13C
	# pdf(file=paste(FigureFolder,"/Ward2010&2012_IsoBiPlot_", GraphLayers[i], "_", ThisVersion, ".pdf", sep=""), width=3.5, height=5.25, pointsize=9)
	png(file=paste(FigureFolder,"/Ward2010&2012_IsoBiPlot_", GraphLayers[i], "_", ThisVersion, ".png", sep=""), units="in", res=300, width=3.5, height=5.25, pointsize=9)
	par(family="Times", las=1, mfrow=c(2,1), mar=c(2.1,4,0.1,0.1), oma=c(1.5,0,0,0), cex=PubCex)
	#Plot the Raw Isotope Data in 3 (scratch that, 2) Combinations

	plot(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", xaxt=ifelse(1==1, "n", "s"), cex.axis=PubCex)#
	mtext(expression(phantom()^2*H), side=2, line=3, las=0, cex=PubCex)

	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]], Plot=T)
	
	points(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*CexScale, xlab=LabelIso[nPlots[1,1]], ylab=LabelIso[nPlots[2,1]], xlim=Lims1[[1]], ylim=Lims1[[2]])
	TaxID <- as.numeric(muDataGroup[,"TaxID"])
	Move17_Index <- which(muDataGroup[,"TaxID"]==17)
	Move17_d_CND <- c(-0.5,0,8)
	TXTmuDataGroup <- muDataGroup
	if(is.element(16, muDataGroup[,"TaxID"])){
		Move16_Index <- which(TXTmuDataGroup[,"TaxID"]==16)
		Move16_d_CND <- c(0,0,5)
		TXTmuDataGroup[Move16_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move16_Index,c("d13C","d15N","dD")] + Move16_d_CND
		
		Move12_Index <- which(TXTmuDataGroup[,"TaxID"]==12)
		Move12_d_CND <- c(0,0,-0.0)
		TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] + Move12_d_CND
		
	}
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,1]]], labels=as.character(TaxID), cex=PubCex*CexScale)
	title(main=LegendTitle[1], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)




	plot(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims2[[1]], ylim=Lims2[[2]], bty="l", xaxt=ifelse(2==1, "n", "s"), cex.axis=PubCex)#
	mtext(expression(phantom()^15*N), side=2, line=3, las=0, cex=PubCex)

	if(2==2){mtext(expression(phantom()^13*C), side=1, line=2.5, outer=FALSE, cex=PubCex)}
	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]], Plot=T)

	points(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*CexScale, xlab=LabelIso[nPlots[1,2]], ylab=LabelIso[nPlots[2,2]], xlim=Lims2[[1]], ylim=Lims2[[2]])
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,2]]], labels=as.character(TaxID), cex=PubCex*CexScale)
	#}

	title(main=LegendTitle[2], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)


	dev.off()
	print(cbind(muDataGroup, "ID"=TaxID)[order(TaxID),])
}


# write.csv(cbind(muDataGroup, "ID"=TaxID)[order(TaxID),], paste(FigureFolder,"/Ward2010&2012_IsoBiPlot_DataTable_",Version,".csv",sep=""))




