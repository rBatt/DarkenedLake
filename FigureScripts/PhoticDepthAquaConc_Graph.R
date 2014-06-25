#Code copied (except some adjustments to the graph aesthetics) from Isotopes2012PreliminarAnalysis_v2.R

rm(list=ls())
graphics.off()

FigureFolder <- paste("Figures/Figures_", "v0.5.2", sep="")
SaveType <- c(".pdf", ".png", ".eps")[2]
Save <- TRUE
Version <- "_v3"


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load("Photic_WardMetab2012_v1.RData")
load("Ward2012_AquaConc.RData")
library("plyr")

Perc1 <- function(x){x[which.min(abs(x[,"PercSurf"]-0.01)),]}
PhoticD <- ddply(Photic, .variables=c("Year","DoY"), .fun=Perc1)
# dev.new(width=7, height=7)
# png(filename="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/PhoticDepthAquaConc_v0.png", res=150, units="in", width=3, height=3, pointsize=10)


setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("PhoticDepthAquaConc_Full_", Version, sep=""), ".pdf", sep=""), height=3, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("PhoticDepthAquaConc_Full_", Version, sep=""), ".png", sep=""), units="in", res=600, height=3, width=3.23)}
	if(SaveType==".eps"){setEPS();postscript(paste(paste("PhoticDepthAquaConc_Full_", Version, sep=""), ".eps", sep=""), width=3.23, height=3)}
}

par(mar=c(2,2,0,2), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")
with(PhoticD[which(PhoticD[,"Year"]==2010),], plot(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="red", ylab="", xlab="", lwd=3, pch=19, las=0))
with(PhoticD[which(PhoticD[,"Year"]==2012),], lines(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="blue", lwd=3, pch=19))
mtext("Day of year", side=1, line=1.2)
mtext("Depth (m)", side=2, line=1.2)
text(x=167, y=2.35, "Photic Depth 2012", col="blue", cex=1)
text(x=167, y=3.75, "Photic Depth 2010", col="red", cex=1)
# legend("topleft", c("Photic Depth (2010)", "Photic Depth (2012)", "[Aquashade]")[c(3,1,2)], text.col=c("red", "blue", "black")[c(3,1,2)], bty="n", inset=c(0.25,-0.03), text.font=1, cex=1)
par(new=TRUE)
AquaDoY <- as.numeric(format.Date(OrderedFileDates, format="%j"))
plot(AquaDoY, EstAquashadeConc, xlim=c(105,240), ylim=c(0,2), type="o", col="black", xaxt="n", yaxt="n", xlab="", ylab="", lwd=3, pch=19)
axis(side=4, las=0)
mtext("Aquashade (ppm)", side=4, line=1.2)
text(x=167, y=1.9, "[Aquashade]", col="black", cex=1)
dev.off()


# =================================================
# = Make graph with Depth profiles of Temp and DO =
# =================================================
#Use day 160 from 2010, and day 159 from 2012
profile0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Ward_TempDO_Profiles_2010&2012.csv")
pLog2010 <- profile0[,"Year"]==2010 & profile0[,"DoY"]==160
pLog2012 <- profile0[,"Year"]==2012 & profile0[,"DoY"]==159
profile <- profile0[pLog2010|pLog2012,]

p10 <- profile[profile[,"Year"]==2010,c("Temp", "DO", "Depth")]
p12 <- profile[profile[,"Year"]==2012,c("Temp", "DO", "Depth")]

pYlim <- c(7,0)
pXlim <- c(0,max(profile[,c("Temp","DO")]))
pXlim10 <- c(0,max(p10[,c("Temp","DO")]))
pXlim12 <- c(0,max(p12[,c("Temp","DO")]))

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("PhoticAquaProfiles_", Version, sep=""), ".pdf", sep=""), height=4.5, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("PhoticAquaProfiles_", Version, sep=""), ".png", sep=""), units="in", res=600, height=4.5, width=3.23)}
	if(SaveType==".eps"){setEPS();postscript(paste(paste("PhoticAquaProfiles_", Version, sep=""), ".eps", sep=""), width=3.23, height=4.5)}
}
layout(mat=matrix(c(1,1,1,2,2,  1,1,1,2,2, 1,1,1,3,3, 1,1,1,3,3), ncol=4))
par(mar=c(3,2,0.2,2), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")
with(PhoticD[which(PhoticD[,"Year"]==2010),], plot(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="red", ylab="", xlab="", lwd=3, pch=19, las=0))
with(PhoticD[which(PhoticD[,"Year"]==2012),], lines(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="blue", lwd=3, pch=19))
mtext("Day of year", side=1, line=1.25)
mtext("Depth (m)", side=2, line=1.15)
text(x=167, y=2.35, "Photic Depth 2012", col="blue", cex=1)
text(x=167, y=3.75, "Photic Depth 2010", col="red", cex=1)

par(new=TRUE)
plot(AquaDoY, EstAquashadeConc, xlim=c(105,240), ylim=c(0,2), type="o", col="black", xaxt="n", yaxt="n", xlab="", ylab="", lwd=3, pch=19)
axis(side=4, las=0, mgp=c(3,0.2,0), tcl=-0.25)
mtext("Aquashade (ppm)", side=4, line=0.9)
text(x=167, y=1.9, "[Aquashade]", col="black", cex=1)

par(mar=c(2.5,2,0,0.2), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")
plot(p10[,c("Temp","Depth")], ylim=pYlim, xlim=pXlim10, xlab="", ylab="", type="l", col="red")
lines(p10[,c("DO","Depth")], lwd=2, col="red")
mtext("Depth (m)", side=2, line=1.15)
text(x=13.5, y=6.5, "Day 160", col="red", cex=1)

par(mar=c(2.5,0.2,0,2), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")
plot(p12[,c("Temp","Depth")], ylim=pYlim, xlim=pXlim12, yaxt="n", xlab="", ylab="", type="l", col="blue")
axis(side=4, mgp=c(3,0.3,0), tcl=-0.25)
lines(p12[,c("DO","Depth")], lwd=2, col="blue")
text(x=16, y=6.5, "Day 159", col="blue", cex=1)
mtext(quote(Dissolved~Oxygen~(mg~L^-1)~and~Temperature~(degree*C)), side=1, line=-1, outer=TRUE)

dev.off()



