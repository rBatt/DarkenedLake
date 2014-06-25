

# =================
# = Load packages =
# =================
library(rLakeAnalyzer)

# ===================
# = Write Functions =
# ===================
Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}

# =====================
# = Load Weather Data =
# =====================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")


# ================================
# = Read in Paul 2010 Thermistor =
# ================================
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
paul10.therm00 <- read.table(paste(datDir,"PaulSondes_2010/Paul_tChain_2010.csv", sep=""), sep=",", header=TRUE)
paul10.therm00[,"datetime"] <- as.POSIXct(paste(paul10.therm00[,"Date"], paul10.therm00[,"Time"]), tz="GMT")

paul10.therm0 <- paul10.therm00[,c("datetime", names(paul10.therm00)[grepl("wtr_[0-9]\\.[0-9]", names(paul10.therm00))])]

paul10.zmix <- ts.meta.depths(paul10.therm0)[,c("datetime", "top")]
names(paul10.zmix)[2] <- "z.mix"

paul10.therm <- merge(paul10.therm0, paul10.zmix, all=TRUE)

# ============================
# = Read in Paul 2010 Sondes =
# ============================
paul10.epi000 <- read.table(paste(datDir,"PaulSondes_2010/PaulSonde2008to10_LinearDOpH_RawFluoro.csv", sep=""), sep=",", header=TRUE)

paul10.epi00 <- paul10.epi000[paul10.epi000[,"Year"]==2010L,]
paul10.epi00 <- paul10.epi00[,c("dayfrac", "Temp", "fixDO", "RawChl", "fixpH", "SpCond")]

paul10.epi0 <- paul10.epi00
names(paul10.epi0) <- c("datetime", "wtr", "DOsat", "RawChl", "fixpH", "SpCond")
paul10.epi0[,"datetime"] <- as.POSIXct(paul10.epi00[,"dayfrac"]*24*60*60, origin="2009-12-31 00:00:00", tz="GMT")
paul10.epi0[,"datetime"] <- LakeMetabolizer:::round.time(paul10.epi0[,"datetime"], "5 minutes")
paul10.epi0[,"do.obs"] <- Sat2Conc(paul10.epi0[,"DOsat"], paul10.epi0[,"wtr"])

paul10.epi.full00 <- paul10.epi0
paul10.epi.full00[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(paul10.epi.full00[,"wtr"], 960)


paul10.epi.full0 <- merge(paul10.epi.full00, paul10.zmix, all.x=TRUE)

# fix na's in zmix after merge
z.mix.isna <- is.na(paul10.epi.full0[,"z.mix"])
paul10.epi.full0[,"z.mix"] <- approx(x=paul10.epi.full0[!z.mix.isna,"datetime"], y=paul10.epi.full0[!z.mix.isna,"z.mix"], xout=paul10.epi.full0[,"datetime"], method="constant", f=0)$y

# add in the wind data
paul10.epi.full0 <- merge(paul10.epi.full0, irr_wnd_2010, all.x=TRUE)

# calculate other values needed for metabolism
paul10.epi.full0[,"wnd"] <- LakeMetabolizer:::wind.scale.base(paul10.epi.full0[,"wnd"], 2) # note that this step is already done for 2012 weather data
paul10.epi.full0[,"k600.cole"] <- LakeMetabolizer:::k.cole.base(paul10.epi.full0[,"wnd"])
paul10.epi.full0[,"k.gas"] <- LakeMetabolizer:::k600.2.kGAS.base(paul10.epi.full0[,"k600.cole"], paul10.epi.full0[,"wtr"])

# cut off gas exchange when zmix is shallow
paul.10.0k <- paul10.epi.full0[,"z.mix"] < 0.7
paul10.epi.full0[paul.10.0k, "k.gas"] <- 0


paul10.epi.full <- paul10.epi.full0
paul10.epi <- paul10.epi.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr", "wnd")]

# =======================
# = Save organized data =
# =======================
save(paul10.therm, paul10.epi.full, paul10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")

write.table(paul10.therm, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.therm.txt", row.names=FALSE)
write.table(paul10.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.epi.full.txt", row.names=FALSE)
write.table(paul10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.epi.txt", row.names=FALSE)



