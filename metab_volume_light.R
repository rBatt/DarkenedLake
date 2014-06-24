

library(plyr)
library(LakeMetabolizer)

source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lightFunctions.R")


# ======================
# = Read in sonde data =
# ======================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/sondes_ward2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/sondes_paul2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/sondes_paul2012.RData")

# ==============================
# = Read in metabolism results =
# ==============================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward10.epi.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward12.epi.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward10.meta.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward12.meta.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/paul10.epi.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/paul12.epi.kal.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/kf.epi.good.RData")

# ===================
# = Calculate light =
# ===================

# Read in light profiles
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
light.prof00 <- read.table(paste(datDir,"PaulWard_Weekly_2010&2012/PaulWard_Light_2010&2012.csv", sep=""), sep=",", header=TRUE)
light.prof00[,"datetime"] <- as.POSIXct(light.prof00[,"datetime"], tz="GMT")
light.prof00 <- light.prof00[complete.cases(light.prof00),]
light.prof0  <- reshape(light.prof00, v.names=c("fracLight"), timevar="depth", idvar=c("lake","datetime"), direction="wide")
names(light.prof0) <- gsub("^fracLight\\.", "irr_", names(light.prof0))


# Calculate the average zmix for each day
w10.light000 <- cbind("lake"="Ward", LakeMetabolizer:::addNAs(ward10.epi.full)[,c("datetime","doy","z.mix", "irr")])
w12.light000 <- cbind("lake"="Ward", LakeMetabolizer:::addNAs(ward12.epi.full)[,c("datetime","doy","z.mix", "irr")])
l10.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul10.epi.full)[,c("datetime","doy","z.mix", "irr")])
l12.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul12.epi.full)[,c("datetime","doy","z.mix", "irr")])
light000 <- rbind(w10.light000, w12.light000, l10.light000, l12.light000)


light000[,"doy"] <- trunc(light000[,"doy"])

light000[,"year"] <- format.Date(light000[,"datetime"], format="%Y")

light00 <- light000[complete.cases(light000),]


light0 <- merge(light00, light.prof0, all=TRUE)



light <- ddply(light0, c("lake","year"), doLight)

light.vals00 <- data.frame(light[,1:4], "val.sonde"=light[,"dz.sonde"], "val.bot"=light[,"dz.bot"]*light[,"irr.bot"]/light[,"irr.sonde"])
light.vals0 <- light.vals00
light.vals0[,"datetime"] <- as.POSIXct(trunc.POSIXt(light.vals0[,"datetime"], "days"), tz="GMT")
light.vals <- ddply(.data=light.vals0, c("lake","year","doy","datetime"), function(x)colMeans(x[,c("val.sonde","val.bot")], na.rm=TRUE))

metab.light0 <- merge(kf.epi.good, light.vals, all.x=TRUE)

metab.light <- metab.light0
metab.light[,"gpp.sonde"] <- metab.light[,"GPP"]*metab.light[,"val.sonde"]
metab.light[,"gpp.bot"] <- metab.light[,"GPP"]*metab.light[,"val.bot"]


metab.light[,"gpp.tot"] <- metab.light[,"gpp.sonde"] + metab.light[,"gpp.bot"]

# ddply(metab.light, c("lake","year"), function(x)colMeans(x[,-c(1:4)]))


save(light, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/light.RData")
save(metab.light, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/metab.light.RData")


