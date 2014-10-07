
rm(list=ls())
graphics.off()

zData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardZoopMass2010&2012.csv")
cData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardChaobMass2010&2012.csv")
zData <- reshape(zData0, varying=list(names(zData0[,4:17])), times=names(zData0[,4:17]), ids=1:nrow(zData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(zData) <- NULL
# zData <- subset(zData, DoY>=143 & Taxon!="Nauplii")
zData <- subset(zData, DoY>=143)
# zSameDays <- which(zData[,"DoY"] >= 143)
# zData <- zData[zSameDays,]
zData[,"Year"] <- as.factor(zData[,"Year"])


cData <- reshape(cData0, varying=list(names(cData0[,4:7])), times=names(cData0[,4:7]), ids=1:nrow(cData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(cData) <- NULL
# cData <- subset(cData, DoY>=143 & Taxon!="FirstInstars")
cData <- subset(cData, DoY>=143)
# cSameDays <- which(cData[,"DoY"] >= 143)
# cData <- cData[cSameDays,]
cData[,"Year"] <- as.factor(cData[,"Year"])



sumzData <- aggregate(zData[,"Mass"], by=list(zData[,"Lake"], zData[,"Year"], zData[,"DoY"]), sum)
names(sumzData) <- c("Lake", "Year", "DoY", "Mass")
zYearMean <- aggregate(sumzData[,"Mass"], by=list(sumzData[,"Lake"], sumzData[,"Year"]), mean)
names(zYearMean) <- c("Lake", "Year", "Mass")
sumcData <- aggregate(cData[,"Mass"], by=list(cData[,"Lake"], cData[,"Year"], cData[,"DoY"]), sum)
names(sumcData) <- c("Lake", "Year", "DoY", "Mass")
cYearMean <- aggregate(sumcData[,"Mass"], by=list(sumcData[,"Lake"], sumcData[,"Year"]), mean)
names(cYearMean) <- c("Lake", "Year", "Mass")


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures")
# dev.new(height=3.5, width=6.5, pointsize=10)
# pdf(file="ZoopChaob_PaulWard_2010&2012.pdf", height=3.5, width=6.5, pointsize=10)
# png(file="ZoopChaob_PaulWard_2010&2012.png", units="in", res=200, height=3.5, width=6.5, pointsize=10)
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("NeatSummary_", Version, sep=""), ".pdf", sep=""), height=7, width=6.81)}
	if(SaveType==".png"){png(file=paste(paste("NeatSummary_", Version, sep=""), ".png", sep=""), units="in", res=200, height=7, width=6.81)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("NeatSummary_", Version, sep=""), ".eps", sep=""), height=7, width=6.81)}
}else{
	dev.new(height=7, width=6.811)
}
setEPS()
postscript(file="ZoopChaob_PaulWard_2010&2012.eps", width=3, height=4.5, pointsize=10)
par(mfrow=c(2,1), mar=c(1,3,0,0),, oma=c(1.5, 0, 0, 0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
# boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
# axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
axis(side=1, at=c(1,4), labels=FALSE)
# legend("topleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext(expression(Zooplankton~Biomass~(g/m^2)), side=2, line=2)
# summary(lm(Mass~Year*Lake, data=sumzData))

# boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
# legend("topleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext(expression(italic(Chaoborus)~spp.~Biomass~(g/m^2)), side=2, line=2)
# summary(lm(Mass~Year*Lake, data=sumcData))
dev.off()


# sumzData2 <- subset(sumzData, Lake=="Ward")
# meanzData2 <- aggregate(sumzData2[,"Mass"], by=list(sumzData2[,"Year"]), mean)
# 
# sumcData2 <- subset(sumcData, Lake=="Ward")
# meancData2 <- aggregate(sumcData2[,"Mass"], by=list(sumcData2[,"Year"]), mean)
# 
# 
# barplot(height=meanzData2[,2])
# barplot(height=meancData2[,2])