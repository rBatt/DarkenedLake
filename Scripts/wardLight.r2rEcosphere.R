library(plyr)
library(LakeMetabolizer)

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/light.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/lightFunctions.R")



irrInt(1500, 0, 4, -log(0.01)/4)

irrInt(1500, 0, 2, -log(0.01)/2)

(irrInt(1500, 0, 3, -log(0.01)/3)) / (irrInt(1500, 0, 9, -log(0.01)/9)) # again, exactly 1/3


light.r2r.1 <- light[light[,"lake"]=="Ward",]
ddply(light.r2r.1, "year", function(x)data.frame(mean=mean(x[,"irr.sonde"], na.rm=TRUE), median=median(x[,"irr.sonde"], na.rm=TRUE), sd=sd(x[,"irr.sonde"], na.rm=TRUE)))




light.r2r.2 <- light[light[,"lake"]=="Ward"&light[,"z.mix"]>=0.5,]
ddply(light.r2r.2, "year", function(x)data.frame(mean=mean(x[,"irr.sonde"], na.rm=TRUE), median=median(x[,"irr.sonde"], na.rm=TRUE), sd=sd(x[,"irr.sonde"], na.rm=TRUE)))


