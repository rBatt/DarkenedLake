
# ===============
# = ZOOPLANKTON =
# ===============
zData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardZoopMass2010&2012.csv")
cData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardChaobMass2010&2012.csv")
zData <- reshape(zData0, varying=list(names(zData0[,4:17])), times=names(zData0[,4:17]), ids=1:nrow(zData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(zData) <- NULL
zData <- subset(zData, DoY>=143)

zData[,"Year"] <- as.factor(zData[,"Year"])

cData <- reshape(cData0, varying=list(names(cData0[,4:7])), times=names(cData0[,4:7]), ids=1:nrow(cData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(cData) <- NULL
cData <- subset(cData, DoY>=143)
cData[,"Year"] <- as.factor(cData[,"Year"])

sumzData <- aggregate(zData[,"Mass"], by=list(zData[,"Lake"], zData[,"Year"], zData[,"DoY"]), sum)
names(sumzData) <- c("Lake", "Year", "DoY", "Mass")
zYearMean <- aggregate(sumzData[,"Mass"], by=list(sumzData[,"Lake"], sumzData[,"Year"]), mean)
names(zYearMean) <- c("Lake", "Year", "Mass")
sumcData <- aggregate(cData[,"Mass"], by=list(cData[,"Lake"], cData[,"Year"], cData[,"DoY"]), sum)
names(sumcData) <- c("Lake", "Year", "DoY", "Mass")
cYearMean <- aggregate(sumcData[,"Mass"], by=list(sumcData[,"Lake"], sumcData[,"Year"]), mean)
names(cYearMean) <- c("Lake", "Year", "Mass")


if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/ZoopChaob_PaulWard_2010&2012.pdf", height=1.75, width=4.33)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/ZoopChaob_PaulWard_2010&2012.png", units="in", res=200, height=1.75, width=4.33)}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/ZoopChaob_PaulWard_2010&2012.eps", height=1.75, width=4.33)}
}else{
	dev.new(height=1.75, width=4.33)
}
par(mfrow=c(1,2), mar=c(1,2,0,0.25), oma=c(0.3, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(Zooplankton~(g/m^2)), side=2, line=1)

par(mar=c(1,2.25,0,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(italic(Chaoborus)~spp.~(g/m^2)), side=2, line=1)
dev.off()



# ============
# = ROUTINES =
# ============
#Taken from L&W_Routine_2010&2012_v3.r and then revised
Data0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
POC0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/POC_PaulWard2010&2012.csv")
ChlaP0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_ChlaProfile_2010&2012.csv")
ChlaP0 <- subset(ChlaP0, Zid!="Meta")
areal <- function(x){
	mean(x[,"Chla"])*max(x[,"Z"])
}
aChl0 <- ddply(ChlaP0, .variables=c("Lake", "Year", "Date"), .fun=areal)
names(aChl0) <- c("Lake", "Year", "Date", "aChla")

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")
PlotNames <- c("L 10", "L 12", "W 10", "W 12")
PML_Data0 <- subset(Data0, Layer=="PML")[,]
PML_Data <- merge(PML_Data0, POC0, all.x=TRUE)
PML_Data[,"C:Chl"] <- PML_Data[,"POC"]/(PML_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "POC","PON","C:Chl", "Temp", "DOC", "DIC", "DO_Conc")
Meta_Data0 <- subset(Data0, Layer=="Meta")[,]
Meta_Data <- merge(Meta_Data0, POC0, all.x=TRUE)
Meta_Data[,"C:Chl"] <- Meta_Data[,"POC"]/(Meta_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "Temp", "DOC", "DIC", "DO_Conc")
Hypo_Data <- subset(Data0, Layer=="Hypo")[,]

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")

# =============================================
# = 4-panel, areal chlorophyll, pCO2, DO, NEP =
# =============================================

if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/LimnoMetab_PaulWard_2010&2012.pdf", height=4, width=4.33)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/LimnoMetab_PaulWard_2010&2012.png", units="in", res=200, height=3.5, width=4.33)}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/LimnoMetab_PaulWard_2010&2012.eps", height=4, width=4.33)}
}else{
	dev.new(height=4, width=4.33)
}
par(mfrow=c(2,2), mar=c(1,2,0.2,0.25), oma=c(0.2, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)

#Chlorophyll
boxplot(aChla~Year+Lake, data=aChl0, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(Chlorophyll~(mg/m^2)), side=2, line=1)


#pCO2
par(mar=c(1,2.25,0.2,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
ppCO2 <- subset(Data0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(italic(p)*CO[2]~(mu*atm)), side=2, line=1)

#DO
par(mar=c(1,2,0.2,0.25), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_2010&2012_Metabolism_v0.2.RData")
boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(DO~("%"*saturation)), side=2, line=1)

#NEP
par(mar=c(1,2.25,0.2,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
pWardPaul_Metabolism <- subset(WardPaul_Metabolism, DoY>=143)
pWardPaul_Metabolism[,"Lake"] <- relevel(pWardPaul_Metabolism[,"Lake"], ref="Paul")
boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(NEP~(mmol~O[2]~m^-3~d^-1)), side=2, line=1)
dev.off()
