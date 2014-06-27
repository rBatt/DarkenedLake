
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/lightFunctions.R")
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
# load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/weather.full.2010.RData")


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


# ======================
# = Read in Light Data =
# ======================
# Read in light profiles
datDir.light <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
light.prof00 <- read.table(paste(datDir.light,"PaulWard_Weekly_2010&2012/PaulWard_Light_2010&2012.csv", sep=""), sep=",", header=TRUE)
light.prof00[,"datetime"] <- as.POSIXct(light.prof00[,"datetime"], tz="GMT")
light.prof00 <- light.prof00[complete.cases(light.prof00),]
light.prof  <- reshape(light.prof00, v.names=c("fracLight"), timevar="depth", idvar=c("lake","datetime"), direction="wide")
names(light.prof) <- gsub("^fracLight\\.", "irr_", names(light.prof))
light.prof[,"year"] <- as.integer(format.Date(light.prof[,"datetime"], format="%Y"))
paul10.light.prof <- light.prof[light.prof[,"lake"]=="Paul"&light.prof[,"year"]==2010L,]
paul10.light.prof <- paul10.light.prof[,!names(paul10.light.prof)%in%c("lake","year")]


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



paul10.epi.full0 <- merge(paul10.epi.full00, paul10.zmix, all.x=TRUE)

# fix na's in zmix after merge
z.mix.isna <- is.na(paul10.epi.full0[,"z.mix"])
paul10.epi.full0[,"z.mix"] <- approx(x=paul10.epi.full0[!z.mix.isna,"datetime"], y=paul10.epi.full0[!z.mix.isna,"z.mix"], xout=paul10.epi.full0[,"datetime"], method="constant", f=0)$y

# add in the wind data
paul10.epi.full0 <- merge(paul10.epi.full0, weather.full.2010, all.x=TRUE)



# ==========================================
# = Calculate Light data for Paul Epi 2010 =
# ==========================================
# Calculate the average zmix for each day
paul10.light000 <- LakeMetabolizer:::addNAs(paul10.epi.full0)[,c("datetime","doy","z.mix", "irr")]
# paul10.epi.full0[,"doy"] <- LakeMetabolizer:::date2doy(paul10.epi.full0[,"datetime"])
# paul10.light000 <- paul10.epi.full0[,c("datetime","doy","z.mix", "irr")]


paul10.light000[,"doy"] <- trunc(paul10.light000[,"doy"])
paul10.light000[,"year"] <- format.Date(paul10.light000[,"datetime"], format="%Y")
paul10.light00 <- paul10.light000[complete.cases(paul10.light000),]
light0 <- merge(paul10.light00, paul10.light.prof, all=TRUE)
paul10.light0 <- doLight(light0)
paul10.light <- paul10.light0[,c("datetime", "irr.sonde", "irr.bot", "dz.sonde", "dz.bot", "kd")]

paul10.epi.full0 <- merge(paul10.epi.full0, paul10.light, all.x=TRUE)


# =========================================
# = Calculate other other params for Paul =
# =========================================
# calculate other values needed for metabolism
paul10.baro <- paul10.epi.full0[,"baro"]
paul10.baro[is.na(paul10.epi.full0[,"baro"])] <- 960
paul10.epi.full0[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(paul10.epi.full0[,"wtr"], paul10.epi.full0[,"baro"])
paul10.epi.full0[,"wnd"] <- LakeMetabolizer:::wind.scale.base(paul10.epi.full0[,"wnd"], 2) # note that this step is already done for 2012 weather data
paul10.epi.full0[,"k600.cole"] <- LakeMetabolizer:::k.cole.base(paul10.epi.full0[,"wnd"])
paul10.epi.full0[,"kgas.cole"] <- LakeMetabolizer:::k600.2.kGAS.base(paul10.epi.full0[,"k600.cole"], paul10.epi.full0[,"wtr"])





# calculate k.read
paul10.sw <- par.to.sw.base(paul10.epi.full0[,"irr"])
paul10.lw <- calc.lw.net.base(
	dateTime=paul10.epi.full0[,"datetime"], 
	sw=paul10.sw, Ts=paul10.epi.full0[,"wtr"], 
	lat=46.28, 
	atm.press=paul10.epi.full0[,"baro"], 
	airT=paul10.epi.full0[,"airTemp"],
	RH=paul10.epi.full0[,"RH"]
	)
	
miss.read <- is.na(paul10.epi.full0[,"wtr"]) | is.na(paul10.epi.full0[,"kd"]) | is.na(paul10.epi.full0[,"baro"]) | is.na(paul10.epi.full0[,"airTemp"]) | is.na(paul10.epi.full0[,"RH"]) | is.na(paul10.epi.full0[,"z.mix"])

paul10.epi.full0[!miss.read,"k600.read"] <- k.read.base(
	wnd.z=10,
	Kd=paul10.epi.full0[!miss.read,"kd"],
	lat=46.28,
	lake.area=17000,
	atm.pres=paul10.epi.full0[!miss.read,"baro"],
	dateTime=paul10.epi.full0[!miss.read,"datetime"],
	Ts=paul10.epi.full0[!miss.read,"wtr"],
	z.aml=paul10.epi.full0[!miss.read,"z.mix"],
	airT=paul10.epi.full0[!miss.read,"airTemp"],
	wnd=paul10.epi.full0[!miss.read,"wnd"],
	RH=paul10.epi.full0[!miss.read,"RH"],
	sw=paul10.sw[!miss.read],
	lwnet=paul10.lw[!miss.read]
	)
	
miss.read2 <- is.na(paul10.epi.full0[,"k600.read"])

dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(2.25, 2.25, 0.5, 0.5), mgp=c(1.25, 0.35, 0), tcl=-0.25, ps=10, family="Times")
plot(paul10.epi.full0[!miss.read2,"datetime"], paul10.epi.full0[!miss.read2,"k600.read"], type="l", xlim=range(paul10.epi.full0[,"datetime"]), xlab="", ylab=bquote(k600~~(m~d^-1)))
lines(paul10.epi.full0[,"datetime"], paul10.epi.full0[,"k600.cole"], type="l", col="red")
legend("topright", legend=c("k.cole", "k.read"), lty=1, lwd=1.5, col=c("red","black"))
plot(log(paul10.epi.full0[!miss.read2,"k600.cole"]), paul10.epi.full0[!miss.read2,"k600.read"], xlab=bquote(log[e]~k600.cole~~(m~d^-1)), ylab=bquote(k600.read~~(m~d^-1)))

paul10.log.cole <- log(paul10.epi.full0[,"k600.cole"])
paul10.lm.read.cole <- lm(paul10.epi.full0[,"k600.read"]~paul10.log.cole)
paul10.read.from.cole <- predict(paul10.lm.read.cole, newdata=data.frame(paul10.log.cole=paul10.log.cole[miss.read2]))
paul10.epi.full0[miss.read2,"k600.read"] <- paul10.read.from.cole
paul10.epi.full0[,"kgas.read"] <- LakeMetabolizer:::k600.2.kGAS.base(paul10.epi.full0[,"k600.read"], paul10.epi.full0[,"wtr"], gas="O2")




# cut off gas exchange when zmix is shallow
paul.10.0k <- paul10.epi.full0[,"z.mix"] < 0.7 & !is.na(paul10.epi.full0[,"z.mix"])
paul10.epi.full0[paul.10.0k, "kgas.cole"] <- 0
paul10.epi.full0[paul.10.0k,"kgas.read"] <- 0L

paul10.epi.full <- paul10.epi.full0
paul10.epi <- paul10.epi.full[,c("datetime", "do.obs", "do.sat", "kgas.cole", "kgas.read", "z.mix",  "irr", "wtr", "wnd")]

# =======================
# = Save organized data =
# =======================
save(paul10.therm, paul10.epi.full, paul10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")

write.table(paul10.therm, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.therm.txt", row.names=FALSE)
write.table(paul10.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.epi.full.txt", row.names=FALSE)
write.table(paul10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul10.epi.txt", row.names=FALSE)



