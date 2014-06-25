#Compare the routine limnology of Paul and Ward Lakes in 2010 and 2012.
#_v0 (05-Feb-2013)

rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/")
Data0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
POC0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/POC_PAULWard2010&2012.csv")

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

Conc2Sat <- function(O2conc, temp){  
	O2mgL=O2conc
	O2sat=(O2mgL*100)/(-0.00000002057759*temp^5 + 0.000002672016*temp^4 + -0.0001884085*temp^3 + 0.009778012*temp^2 + -0.4147241*temp + 14.621)
	return(O2sat)
}
Meta_Data[,"DO_sat"] <- Conc2Sat(Meta_Data[,"DO_Conc"], Meta_Data[,"Temp"])

aggregate(Meta_Data[,"DO_Conc"],by=list(Meta_Data[,"Lake"], Meta_Data[,"Year"]), FUN=mean)
aggregate(Meta_Data[,"DO_sat"],by=list(Meta_Data[,"Lake"], Meta_Data[,"Year"]), FUN=mean)
aggregate(Meta_Data[,"DO_sat"],by=list(Meta_Data[,"Lake"], Meta_Data[,"Year"]), FUN=sd)


aggregate(PML_Data[,"C:Chl"],by=list(PML_Data[,"Lake"], PML_Data[,"Year"]), FUN=mean, na.rm=TRUE)

png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/DOM_PaulWard_2010&2012.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
par(mar=c(2.5,3.5,0.5,0.5), ps=10)
pEpiDOM <- subset(Data0, Layer=="PML", select=c("Lake","Year","DOC"))
boxplot(DOC~Year+Lake, data=pEpiDOM, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
legend("topleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext("DOC (mg/L)", side=2, line=2.25)
dev.off()


png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/pCO2_PaulWard_2010&2012.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
par(mar=c(2.5,3.5,0.5,0.5), ps=10)
ppCO2 <- subset(Data0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
legend("topleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext(expression(Surface~water~italic(p)*CO[2]~(mu*atm)), side=2, line=2.25)
dev.off()