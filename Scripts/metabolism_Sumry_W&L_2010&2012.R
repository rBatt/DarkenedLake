library(plyr)
library(LakeMetabolizer)

# ======================
# = Read in sonde data =
# ======================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2012.RData")

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




# ===============
# = Light stuff =
# ===============
# Read in light profiles
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
light.prof00 <- read.table(paste(datDir,"PaulWard_Weekly_2010&2012/PaulWard_Light_2010&2012.csv", sep=""), sep=",", header=TRUE)
light.prof00[,"datetime"] <- as.POSIXct(light.prof00[,"datetime"], tz="GMT")
light.prof00 <- light.prof00[complete.cases(light.prof00),]
light.prof0  <- reshape(light.prof00, v.names=c("fracLight"), timevar="depth", idvar=c("lake","datetime"), direction="wide")
names(light.prof0) <- gsub("^fracLight\\.", "irr_", names(light.prof0))

grab1perc <- function(x0){
	x <- x0[,-c(1,2)]
	x <- abs(x-0.01)
	x <- apply(x, 1, function(z)names(which.min(z)))
	data.frame(x0[,1:2], "z1perc"=as.numeric(gsub("[^ 0-9 \\.]", "", x)))
}

aprxProf <- function(x0){
	seqStart <- min(x0[,"datetime"])
	seqStop <- max(x0[,"datetime"])
	t.Seq <- seq(seqStart, seqStop, by=60*60*24)
	x <- approx(x=x0[,"datetime"], y=x0[,"z1perc"], xout=t.Seq, rule=2, method="constant", f=0)
	data.frame("lake"=x0[1,"lake"], "datetime"=x$x, "z1perc"=x$y)
}

z1perc0 <- grab1perc(light.prof0)
z1perc0[,"year"] <- as.integer(format.Date(z1perc0[,"datetime"], format="%Y"))


kf.epi.good.z1 <- merge(kf.epi.good, z1perc0, all=TRUE)
z1perc <- ddply(kf.epi.good.z1, c("lake", "year"), aprxProf)[,c("lake","datetime", "z1perc")]
z1perc[,"lake"] <- as.character(z1perc[,"lake"])
kf.epi.good.z1 <- merge(kf.epi.good, z1perc, all.x=TRUE)
kf.epi.good.z1[,c("GPP.areal", "R.areal", "NEP.areal")] <- kf.epi.good.z1[,c("GPP", "R", "NEP")]* kf.epi.good.z1[,"z1perc"]



ddply(kf.epi.good.z1, c("lake","year"), function(x)colMeans(x[,c("GPP","GPP.areal","R","R.areal","NEP","NEP.areal","z1perc")])/(15.999*2)*1000)
