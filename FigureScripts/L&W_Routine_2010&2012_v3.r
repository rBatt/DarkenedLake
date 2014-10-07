#Compare the routine limnology of Paul and Ward Lakes in 2010 and 2012.
#_v0 (05-Feb-2013)

rm(list=ls())
graphics.off()
library("plyr")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/")
Data0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
POC0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/POC_PaulWard2010&2012.csv")
ChlaP0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_ChlaProfile_2010&2012.csv")
ChlaP0 <- subset(ChlaP0, Zid!="Meta")
areal <- function(x){
	mean(x[,"Chla"])*max(x[,"Z"])
}
aChl0 <- ddply(ChlaP0, .variables=c("Lake", "Year", "Date"), .fun=areal)
names(aChl0) <- c("Lake", "Year", "Date", "aChla")
# aChl0 <-aggregate(ChlaP0[,"Chla"], by=list(ChlaP0[,"Lake"], ChlaP0[,"Year"], ChlaP0[,"Date"]), FUN=sum)

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

aggregate(PML_Data[,"C:Chl"],by=list(PML_Data[,"Lake"], PML_Data[,"Year"]), FUN=mean, na.rm=TRUE)

aggregate(PML_Data[,"POC"],by=list(PML_Data[,"Lake"], PML_Data[,"Year"]), FUN=mean, na.rm=TRUE)

aggregate(aChl0[,"aChla"], by=list(aChl0[,"Lake"], aChl0[,"Year"]), FUN=mean, na.rm=TRUE)
(2600/6615)*60


# =============================================
# = 4-panel, areal chlorophyll, pCO2, DO, NEP =
# =============================================

# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/NEP_PaulWard_2010&2012_Metabolism_v0.1.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)

# dev.new(width=3, height=7, pointsize=10)
# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Chl_DO_pCO2_NEP_PaulWard_2010&2012.png", units="in", res=150, height=6, width=3, pointsize=10)
# jpeg(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Chl_DO_pCO2_NEP_PaulWard_2010&2012.jpg", units="in", quality=100, res=200, height=6, width=3, pointsize=10)
setEPS()
postscript(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Chl_DO_pCO2_NEP_PaulWard_2010&2012.eps", width=3, height=6, pointsize=10)
par(mfrow=c(4,1), mar=c(1,4,0.5,0.5), oma=c(1, 0, 0.5, 0), ps=10, cex=1)

#Chlorophyll
# boxplot(aChla~Year+Lake, data=aChl0, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(aChla~Year+Lake, data=aChl0, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(expression(Chlorophyll~(mg/m^2)), side=2, line=2.25)
# legend("top", legend="A", text.font=2)

#pCO2
ppCO2 <- subset(Data0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
# boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)

axis(side=1, at=c(1,4), labels=FALSE)
mtext(expression(italic(p)*CO[2]~(mu*atm)), side=2, line=2.25)
# legend("top", legend="B", text.font=2)

#DO
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_2010&2012_Metabolism_v0.2.RData")
# boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)

axis(side=1, at=c(1,4), labels=FALSE)
mtext(expression(DO~("%"*saturation)), side=2, line=2.25)
# legend("top", legend="C", text.font=2)

#NEP
# load(file="/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_2010&2012_Metabolism_v0.2.RData")
pWardPaul_Metabolism <- subset(WardPaul_Metabolism, DoY>=143)
pWardPaul_Metabolism[,"Lake"] <- relevel(pWardPaul_Metabolism[,"Lake"], ref="Paul")
# boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)

axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(expression(NEP~(mmol~O[2]~m^-3~d^-1)), side=2, line=2.25)
# legend("top", legend="D", text.font=2)
dev.off()

aggregate(pWardPaul_Metabolism[,"NEP"], by=list(pWardPaul_Metabolism[,"Lake"], pWardPaul_Metabolism[,"Year"]), FUN=mean, na.rm=TRUE)

DOsummary <- aggregate(AllDO[,"MeanDO"], by=list(AllDO[,"Lake"], AllDO[,"Year"]), FUN=mean, na.rm=TRUE)
cbind(DOsummary, "localDO"=DOsummary[,"x"]+(100-94.6))

# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/DOpercent_PaulWard_2010&2012_Metabolism_v0.1.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
# par(mar=c(2.5,3.5,0.5,0.5), ps=10)



# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/DOM_PaulWard_2010&2012.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
# par(mar=c(2.5,3.5,0.5,0.5), ps=10)
# pEpiDOM <- subset(Data0, Layer=="PML", select=c("Lake","Year","DOC"))
# boxplot(DOC~Year+Lake, data=pEpiDOM, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
# axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
# legend("topleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
# mtext("DOC (mg/L)", side=2, line=2.25)
# dev.off()


# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/pCO2_PaulWard_2010&2012.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
# par(mar=c(2.5,3.5,0.5,0.5), ps=10)

# dev.off()

# png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/ArealChla_PaulWard_2010&2012.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
# par(mar=c(2.5,3.5,0.5,0.5), ps=10)

