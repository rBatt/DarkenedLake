#Ryan Batt
#21-May-2011
#Create Bi-Plots from isotope data
#Data can be entered in 2 formats, 1) Data that does not have mean and sd calculated, or 2)The argument "Data" is just the means and the groups, and the "SDs" argument is the standard deviations entered as a data.frame(). See examples below with "PI Data" ;)
 
rm(list=ls())
graphics.off()
Version <- "v0.1.0"
FigureFolder <- paste("Figures_", Version, sep="")
# YearMix <- 2012
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")

#Load the 
load(paste("Data+Phyto_NoTree_", 2012, "_", Version, ".RData",sep=""))
Sim_P_dX_Epi_Obs[,"Year"] <- 2012
Sim_P_dX_Meta_Obs[,"Year"] <- 2012
Sim_P_dX_Epi_Obs_2012 <- Sim_P_dX_Epi_Obs
Sim_P_dX_Meta_Obs_2012 <- Sim_P_dX_Meta_Obs
rm(list=c("Sim_P_dX_Epi_Obs", "Sim_P_dX_Meta_Obs"))

load(paste("Data+Phyto_NoTree_", 2010, "_", Version, ".RData",sep=""))
Sim_P_dX_Epi_Obs[,"Year"] <- 2010
Sim_P_dX_Meta_Obs[,"Year"] <- 2010
Sim_P_dX_Epi_Obs_2010 <- Sim_P_dX_Epi_Obs
Sim_P_dX_Meta_Obs_2010 <- Sim_P_dX_Meta_Obs
rm(list=c("Sim_P_dX_Epi_Obs", "Sim_P_dX_Meta_Obs"))
# Data <- subset(DataRaw, Trophic>=0 & SampleID!="O-0362" & Taxon!="POC" & Taxon!="Nija" & Taxon!="Tree"  & Taxon!="DOC", select=c("Trophic","Taxon","d13C","d15N","dD"))

CombYears <- which(is.element(DataRaw[,"Type"], c("Macrophyte", "Terrestrial"))) #I don't want to make separate points for samples types that I don't really expect to change between years, such as the macrophytes or the terrestrial samples.  I just want to assume that these were the same in the different years, while the consumers, DOM, POM etc. may have changed.
DataRaw[CombYears,"Year"] <- 2010

Data <- droplevels(subset(DataRaw, Trophic>=0 & !is.element(SampleID, c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")) & Taxon!="POM" & Taxon!="Nija" & Taxon!="YWP", select=c("Year","Trophic","Taxon","d13C","d15N","dD")))

Sim_P_dX_Epi_Obs <- rbind(Sim_P_dX_Epi_Obs_2010, Sim_P_dX_Epi_Obs_2012)
Sim_P_dX_Meta_Obs <- rbind(Sim_P_dX_Meta_Obs_2010, Sim_P_dX_Meta_Obs_2012)

DataNew <- Data
UniqueTaxa <- as.character(unique(Data[,"Taxon"]))
RealTaxa <- as.character(UniqueTaxa) #c("Helisoma trivolvis", "DOM", "Perca flavescens", "Lepomis gibbosus", "Pimaphales promelas", "Periphyton", "Potamogeton pussillus", "Chara sp.", "Nymphaea odorata", "Carex sp.", "Alnus incana subsp. Rugosa", "Skistodiaptomus oregonensis", "Chaoborus spp.", "Nuphar lutea", "Brasenia schreberi", "Trees") #"DOM", 
# for(i in 1:length(UniqueTaxa)){
	
	# DataNew[,"Taxon"] <- gsub(UniqueTaxa[i], RealTaxa[i], as.character(DataNew[,"Taxon"]))
	
# }

#gsub(UniqueTaxa, RealTaxa, as.character(Data[,"Taxon"]))

MetaPhytos <- data.frame("Trophic"=rep(0, nrow(Sim_P_dX_Meta_Obs)),"Taxon"=rep("MetaPhyto", nrow(Sim_P_dX_Meta_Obs)),  "d13C"=Sim_P_dX_Meta_Obs[,2], "d15N"=Sim_P_dX_Meta_Obs[,3], "dD"=Sim_P_dX_Meta_Obs[,4], "Year"=Sim_P_dX_Meta_Obs[,"Year"]) 
EpiPhytos <- data.frame("Trophic"=rep(0, nrow(Sim_P_dX_Epi_Obs)),"Taxon"=rep("EpiPhyto", nrow(Sim_P_dX_Epi_Obs)),  "d13C"=Sim_P_dX_Epi_Obs[,2], "d15N"=Sim_P_dX_Epi_Obs[,3], "dD"=Sim_P_dX_Epi_Obs[,4], "Year"=Sim_P_dX_Epi_Obs[,"Year"]) 


Sources_Cons <- rbind(DataNew, MetaPhytos, EpiPhytos)


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

muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic")], by=as.list(Data[,Grouping]), mean)
muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))
print(muDataGroup)

if(is.null(SDs)){
	sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic")], by=as.list(Data[,Grouping]), sd)
	sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
	}else{
		sdDataGroup <- SDs
}

Consumers <- c("Calanoid", "Chaoborus", "FHM", "PKS", "Snail", "YWP")
LegendTitle <- c("A", "B")
# dev.new(height=7.25, width=7.25) #dD vs. d13C
pdf(file=paste(FigureFolder,"/Ward2010&2012_IsoBiPlot_",Version,".pdf",sep=""), width=5.25, height=5.25, pointsize=9)
# dev.new(width=5.25, height=5.25, pointsize=9)
par(family="Times", las=1, mfrow=c(2,1), mar=c(2.1,4,0.1,0.1), oma=c(1.5,0,0,0), cex=PubCex)
#Plot the Raw Isotope Data in 3 (scratch that, 2) Combinations

print(nPlots)
nPlots <- nPlots[,c(2,1)]
print(nPlots)

# rm(list="i")
#for(i in 1:length(nPlots[1,])){
Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]])
plot(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch="", col=c("darkred","darkblue")[as.factor(is.element(muDataGroup[,1], Consumers))], cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", xaxt=ifelse(1==1, "n", "s"), cex.axis=PubCex)#

mtext(expression(phantom()^2*H), side=2, line=3, las=0, cex=PubCex)

# if(Legend==TRUE){legend(LegendLoc, legend=paste(pch=as.character(1:nrow(muDataGroup)), muDataGroup[,GroupNames], sep=" "), bty="n")}

XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]], Plot=T)

points(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*0.8, xlab=LabelIso[nPlots[1,1]], ylab=LabelIso[nPlots[2,1]], xlim=Lims1[[1]], ylim=Lims1[[2]])
TaxID <- as.numeric(muDataGroup[,"Taxon"])
# text(muDataGroup[,IsoNames[nPlots[1:2,1]]], labels=as.character(1:nrow(muDataGroup)), cex=PubCex)
text(muDataGroup[,IsoNames[nPlots[1:2,1]]], labels=as.character(TaxID), cex=PubCex*0.8)
#}

title(main=LegendTitle[1], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)








Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]])
plot(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch="", col=c("darkred","darkblue")[as.factor(is.element(muDataGroup[,1], Consumers))], cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", xaxt=ifelse(2==1, "n", "s"), cex.axis=PubCex)#

mtext(expression(phantom()^15*N), side=2, line=3, las=0, cex=PubCex)

if(2==2){mtext(expression(phantom()^13*C), side=1, line=2.5, outer=FALSE, cex=PubCex)}

# if(Legend==TRUE){legend(LegendLoc, legend=paste(pch=as.character(1:nrow(muDataGroup)), muDataGroup[,GroupNames], sep=" "), bty="n")}

XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]], Plot=T)

points(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7, xlab=LabelIso[nPlots[1,2]], ylab=LabelIso[nPlots[2,2]], xlim=Lims1[[1]], ylim=Lims1[[2]])
# text(muDataGroup[,IsoNames[nPlots[1:2,2]]], labels=as.character(1:nrow(muDataGroup)), cex=PubCex)
text(muDataGroup[,IsoNames[nPlots[1:2,2]]], labels=as.character(TaxID), cex=PubCex)
#}

title(main=LegendTitle[2], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)


dev.off()


# setwd("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2010Analysis/Figures_v8.3")
# dev2bitmap(file="Ward2010_IsoBiPlot_v8.3.tif", type="tiffgray",height=7.25, width=7.25, res=200, font="Times", method="pdf", pointsize=12)










#Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,i]]], muDataGroup[,IsoNames[nPlots[2,i]]], sdDataGroup[,IsoNames[nPlots[1,i]]], sdDataGroup[,IsoNames[nPlots[2,i]]])
#plot(muDataGroup[,IsoNames[nPlots[1:2,i]]], pch="", col=c("darkred","darkblue")[as.factor(is.element(muDataGroup[,1], Consumers))], cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", #xaxt=ifelse(i==1, "n", "s"), cex.axis=1/0.85)#
#
#mtext(list(expression(delta^15*N), expression(delta^2*H))[[i]], side=2, line=3.5, las=0, cex=1/0.85)
#
#if(i==1){mtext(expression(delta^13*C), side=1, line=2.5, outer=FALSE, cex=1/0.85)}
#
#if(Legend==TRUE){legend(LegendLoc, legend=paste(pch=as.character(1:nrow(muDataGroup)), muDataGroup[,GroupNames], sep=" "), bty="n")}
#
#XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,i]]], muDataGroup[,IsoNames[nPlots[2,i]]], sdDataGroup[,IsoNames[nPlots[1,i]]], sdDataGroup[,IsoNames[nPlots[2,i]]], Plot=T)
#
#points(muDataGroup[,IsoNames[nPlots[1:2,i]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Trophic"]==0, "white", "gray85"), cex=2.7, xlab=LabelIso[nPlots[1,i]], #ylab=LabelIso[nPlots[2,i]], xlim=Lims1[[1]], ylim=Lims1[[2]])
#text(muDataGroup[,IsoNames[nPlots[1:2,i]]], labels=as.character(1:nrow(muDataGroup)), cex=1)
#}
#
#title(main=LegendTitle[i], adj=0, line=-3, cex.main=1/0.85*1.2)
	

#}#End function


#BiPlot(Data=Sources_Cons, PlotLabels=FALSE, Legend=FALSE)

write.csv(cbind(muDataGroup, "ID"=TaxID)[order(TaxID),], paste(FigureFolder,"/Ward2010&2012_IsoBiPlot_DataTable_",Version,".csv",sep=""))




