#Code copied (except some adjustments to the graph aesthetics) from Isotopes2012PreliminarAnalysis_v2.R

rm(list=ls())
graphics.off()


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load("Photic_WardMetab2012_v1.RData")
load("Ward2012_AquaConc.RData")


Perc1 <- function(x){x[which.min(abs(x[,"PercSurf"]-0.01)),]}
PhoticD <- ddply(Photic, .variables=c("Year","DoY"), .fun=Perc1)
dev.new(width=7, height=7)
par(mar=c(4,4,1,4))
with(PhoticD[which(PhoticD[,"Year"]==2010),], plot(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="red", ylab="", xlab=""))
with(PhoticD[which(PhoticD[,"Year"]==2012),], lines(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="blue"))
mtext("Day of year", side=1, line=2.5)
mtext("Depth (m)", side=2, line=2.5)
legend("bottomleft", c("Photic Depth (2010)", "Photic Depth (2012)", "[Aquashade]"), text.col=c("red", "blue", "darkturquoise"), bty="n", inset=c(-0.03,0))
par(new=TRUE)
AquaDoY <- as.numeric(format.Date(OrderedFileDates, format="%j"))
plot(AquaDoY, EstAquashadeConc, xlim=c(105,240), ylim=c(0,2), type="o", col="darkturquoise", xaxt="n", yaxt="n", xlab="", ylab="")
axis(side=4)
mtext("Aquashade (ppm)", side=4, line=2.5)