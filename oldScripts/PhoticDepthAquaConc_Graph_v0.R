#Code copied (except some adjustments to the graph aesthetics) from Isotopes2012PreliminarAnalysis_v2.R

rm(list=ls())
graphics.off()


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load("Photic_WardMetab2012_v1.RData")
load("Ward2012_AquaConc.RData")


Perc1 <- function(x){x[which.min(abs(x[,"PercSurf"]-0.01)),]}
PhoticD <- ddply(Photic, .variables=c("Year","DoY"), .fun=Perc1)
# dev.new(width=7, height=7)
png(filename="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/PhoticDepthAquaConc_v0.png", res=150, units="in", width=3, height=3, pointsize=10)
par(mar=c(3,3,1,3), ps=10)
with(PhoticD[which(PhoticD[,"Year"]==2010),], plot(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="red", ylab="", xlab="", lwd=3, pch=19, las=0))
with(PhoticD[which(PhoticD[,"Year"]==2012),], lines(DoY, Depth, type="o", ylim=c(5,0), xlim=c(105,240), col="blue", lwd=3, pch=19))
mtext("Day of year", side=1, line=2)
mtext("Depth (m)", side=2, line=2)
text(x=167, y=2.35, "Photic Depth 2012", col="blue", cex=1)
text(x=167, y=3.75, "Photic Depth 2010", col="red", cex=1)
# legend("topleft", c("Photic Depth (2010)", "Photic Depth (2012)", "[Aquashade]")[c(3,1,2)], text.col=c("red", "blue", "black")[c(3,1,2)], bty="n", inset=c(0.25,-0.03), text.font=1, cex=1)
par(new=TRUE)
AquaDoY <- as.numeric(format.Date(OrderedFileDates, format="%j"))
plot(AquaDoY, EstAquashadeConc, xlim=c(105,240), ylim=c(0,2), type="o", col="black", xaxt="n", yaxt="n", xlab="", ylab="", lwd=3, pch=19)
axis(side=4, las=0)
mtext("Aquashade (ppm)", side=4, line=2)
text(x=167, y=1.9, "[Aquashade]", col="black", cex=1)
dev.off()