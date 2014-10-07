#Compare the routine limnology of Paul and Ward Lakes in 2010 and 2012.
#_v0 (05-Feb-2013)

rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/")
Data0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
POC0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/POC_Ward2010&2012.csv")

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")
PlotNames <- c("L 10", "L 12", "W 10", "W 12")
PML_Data0 <- subset(Data0, Layer=="PML")[,]
PML_Data <- merge(PML_Data0, POC0, all.x=TRUE)
PML_Data[,"C:Chl"] <- PML_Data[,"POC"]/(PML_Data[,"Chla"]/1000)
dev.new(width=8, height=8)
par(mar=c(3,4,0.5,0.5), ps=10, cex=1, mfrow=c(4,4))
for(i in 1:length(VarAnalyze)){
	TempoDat <- subset(PML_Data, select=c("Lake", "Year", "Date", "Layer", "Depth", VarAnalyze[i]))[,]
	names(TempoDat) <- c("Lake", "Year", "Date", "Layer", "Depth", "TempoVar")
	boxplot(TempoVar~Year*Lake, data=TempoDat, names=PlotNames)
	mtext(VarAnalyze[i], side=2, line=2.5)
	if(i==1){legend("topleft", legend="Epi")}

}


VarAnalyze <- c("Color", "Chla", "POC","PON","C:Chl", "Temp", "DOC", "DIC", "DO_Conc")
# Meta_Data <- subset(Data0, Layer=="Meta")[,]
Meta_Data0 <- subset(Data0, Layer=="Meta")[,]
Meta_Data <- merge(Meta_Data0, POC0, all.x=TRUE)
Meta_Data[,"C:Chl"] <- Meta_Data[,"POC"]/(Meta_Data[,"Chla"]/1000)
dev.new(width=8, height=6)
par(mar=c(3,4,0.5,0.5), ps=10, cex=1, mfrow=c(3,3))
for(i in 1:length(VarAnalyze)){
	TempoDat <- subset(Meta_Data, select=c("Lake", "Year", "Date", "Layer", "Depth", VarAnalyze[i]))[,]
	names(TempoDat) <- c("Lake", "Year", "Date", "Layer", "Depth", "TempoVar")
	boxplot(TempoVar~Year*Lake, data=TempoDat, names=PlotNames)
	mtext(VarAnalyze[i], side=2, line=2.5)
	if(i==1){legend("topleft", legend="Meta")}
}


VarAnalyze <- c("Color", "Chla", "Temp", "DOC", "DIC", "DO_Conc")
Hypo_Data <- subset(Data0, Layer=="Hypo")[,]
dev.new(width=8, height=5)
par(mar=c(3,4,0.5,0.5), ps=10, cex=1, mfrow=c(2,3))
for(i in 1:length(VarAnalyze)){
	TempoDat <- subset(Hypo_Data, select=c("Lake", "Year", "Date", "Layer", "Depth", VarAnalyze[i]))[,]
	names(TempoDat) <- c("Lake", "Year", "Date", "Layer", "Depth", "TempoVar")
	boxplot(TempoVar~Year*Lake, data=TempoDat, names=PlotNames)
	mtext(VarAnalyze[i], side=2, line=2.5)
	if(i==1){legend("topleft", legend="Hypo")}

}
aggregate(Hypo_Data[,"DOC"], by=list(Hypo_Data[,"Year"],Hypo_Data[,"Lake"]), FUN=mean, na.rm=TRUE)

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")
# write.table(t(aggregate(PML_Data[,VarAnalyze], by=list(PML_Data[,"Lake"], PML_Data[,"Year"]), FUN=mean, na.rm=TRUE)), file="WardPaulWeeklyTable.csv", sep=",")

aggregate(Meta_Data[,"DOC"],by=list(Meta_Data[,"Lake"], Meta_Data[,"Year"]), FUN=mean)