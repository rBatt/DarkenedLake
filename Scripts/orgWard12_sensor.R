
# Weather, sonde, and thermistor data from Ward Lake in 2012
library("plyr")
library("rLakeAnalyzer")

# detach(package:LakeMetabolizer, unload=TRUE)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")


source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/lightFunctions.R")

# =====================
# = Load Weather Data =
# =====================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/weather.full.2012.RData")


# ===============
# = Manual Zmix =
# ===============
w12.manZ.date0 <- c(105, 110, 117, 124, 131, 138, 145, 152, 159, 166, 173, 180, 187, 194, 201, 208, 215, 222, 229, 236, 241)
w12.manZ.date <- as.POSIXct((w12.manZ.date0)*24*60*60, origin="2011-12-31", tz="GMT") # confirm: as.POSIXct("105-2012", format="%j-%Y")
w12.manZ00 <- c(1.0,2.0,2.0,0.5,0.5, 1.0, 1.0, 1.5, 0.25, 1.0, 1.0, 0.5, 0.5, 1, 1, 1, 1, 1.5, 1.5, 1, 0.25)
w12.manZ.interpDate <- seq(from=as.POSIXct("2012-04-14 00:00:00", tz="GMT"), to=as.POSIXct("2012-08-28 23:55:00", tz="GMT"), by=60*5)
w12.manZ0 <- approx(x=w12.manZ.date, y=w12.manZ00, xout=w12.manZ.interpDate, method="constant", f=0, rule=2)
w12.manZ <- data.frame("datetime"=w12.manZ0$x, "manZ"=w12.manZ0$y)



# =======================
# = Temp Profile/ z.mix =
# =======================
Therms0 <- data.frame("Year"=2012, read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2012/CompiledWardTherms2012_v2.txt", header=TRUE))
ThermDepths <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3)
t.wtr.names <- paste("wtr", as.character(ThermDepths), sep="_")
names(Therms0) <- c("year", "doy", t.wtr.names)
Therms0[,"datetime"] <- as.POSIXct(Therms0[,"doy"]*24*60*60, origin="2012-01-01", tz="GMT")

Therms <- Therms0[,c("datetime",t.wtr.names)]
Therms[,"datetime"] <- LakeMetabolizer:::round.time(Therms[,"datetime"], "5 minutes")
Therms <- Therms[!duplicated(Therms[,"datetime"]),]
ward12.therm <- Therms

ward12.zmix000 <- ts.meta.depths(Therms)
ward12.zmix00 <- ward12.zmix000
ward12.zmix00_is0 <- ward12.zmix00[,"top"] < 0.25 & !is.na(ward12.zmix00[,"top"])
ward12.zmix00[ward12.zmix00_is0, "top"] <- 0.25

ward12.zmix0 <- ward12.zmix00[,c("datetime", "top")]
names(ward12.zmix0) <- c("datetime", "z.mix")

ward12.therm <- merge(Therms, ward12.zmix0, all=TRUE)

w12.bothZ <- merge(ward12.zmix0, w12.manZ, all=TRUE)
w12z.isNA <- is.na(w12.bothZ[,"z.mix"])
ward12.zmix <- w12.bothZ[,c("datetime", "z.mix")]
ward12.zmix[w12z.isNA, "z.mix"] <- w12.bothZ[w12z.isNA, "manZ"]


# ======================
# = Read in Light Data =
# ======================
# Read in light profiles
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
light.prof00 <- read.table(paste(datDir,"PaulWard_Weekly_2010&2012/PaulWard_Light_2010&2012.csv", sep=""), sep=",", header=TRUE)
light.prof00[,"datetime"] <- as.POSIXct(light.prof00[,"datetime"], tz="GMT")
light.prof00 <- light.prof00[complete.cases(light.prof00),]
light.prof  <- reshape(light.prof00, v.names=c("fracLight"), timevar="depth", idvar=c("lake","datetime"), direction="wide")
names(light.prof) <- gsub("^fracLight\\.", "irr_", names(light.prof))
light.prof[,"year"] <- as.integer(format.Date(light.prof[,"datetime"], format="%Y"))
w12.light.prof <- light.prof[light.prof[,"lake"]=="Ward"&light.prof[,"year"]==2012L,]
w12.light.prof <- w12.light.prof[,!names(w12.light.prof)%in%c("lake","year")]



# ==============
# = Sonde data =
# ==============

dat.path <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2012/"
Key2012 <- read.csv(paste(dat.path, "WardSondes2012_FileKey.csv", sep=""))
ShalKey2012 <- Key2012[Key2012[,"Depth"]=="A",]
DeepKey2012 <- Key2012[Key2012[,"Depth"]=="B",]

Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}



# =======================
# = Read Ward Epi 2012 =
# =======================
for(i in 1:nrow(ShalKey2012)){
	tName <- as.character(ShalKey2012[i,"Name"])
	tFormat <- as.character(ShalKey2012[i,"Format"])
	
	if(tFormat=="CDF"){
		if(tName=="A14APR12.CDF"){
			ReadClassCDF <- c("character", "character", rep("numeric", 8))
			tDat <- read.table(paste(dat.path, tName, sep=""), sep=",", header=TRUE, skip=1, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "Battery")
			tDat[,"do.obs"] <- Sat2Conc(tDat[,"DOsat"], tDat[,"wtr"])
		}
	
	# reformat datetime for .CDF files
	tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%m/%d/%y %T", tz="GMT")
	}
	
	
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste(dat.path, tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
		names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "do.obs", "Chl_conc", "Chl_RFU", "Battery")
		tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%Y/%m/%d %T", tz="GMT")
	}
	
	
	tDat <- tDat[,c("datetime", "wtr", "do.obs", "DOsat", "SpCond", "pH", "Chl_conc")]
	if(i == 1){
		ward12.epi0 <- tDat
	}else{
		ward12.epi0 <- rbind(ward12.epi0, tDat)
	}
}


# =============================
# = Format and merge Ward Epi =
# =============================
ward12.epi0[,"datetime"] <- LakeMetabolizer:::round.time(ward12.epi0[,"datetime"], "5 minutes")
ward12.epi.full0 <- merge(ward12.epi0, ward12.zmix, all.x=TRUE)
# ward12.epi.full <- merge(ward12.epi.full0, irr_wnd_2012, all.x=TRUE)
ward12.epi.full <- merge(ward12.epi.full0, weather.full.2012, all.x=TRUE)


# ==========================================
# = Calculate Light data for Ward Epi 2012 =
# ==========================================
# Calculate the average zmix for each day
w12.light000 <- LakeMetabolizer:::addNAs(ward12.epi.full)[,c("datetime","doy","z.mix", "irr")]

w12.light000[,"doy"] <- trunc(w12.light000[,"doy"])
w12.light000[,"year"] <- format.Date(w12.light000[,"datetime"], format="%Y")
w12.light00 <- w12.light000[complete.cases(w12.light000),]
light0 <- merge(w12.light00, w12.light.prof, all=TRUE)
w12.light0 <- doLight(light0)
w12.light <- w12.light0[,c("datetime", "irr.sonde", "irr.bot", "dz.sonde", "dz.bot", "kd")]

ward12.epi.full <- merge(ward12.epi.full, w12.light, all.x=TRUE)


# ==========================================
# = Calculate other params for Ward Epi 10 =
# ==========================================
# wind is already scaled above (when combining UNDERC and Peter)

# ward12.epi.k.cole <- k.cole.base(ward12.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.epi.full[,"k600.cole"] <- k.cole.base(ward12.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.epi.full[,"kgas.cole"] <- k600.2.kGAS.base(ward12.epi.full[,"k600.cole"], ward12.epi.full[,"wtr"], gas="O2")

w12.baro <- ward12.epi.full[,"baro"]
w12.baro[is.na(ward12.epi.full[,"baro"])] <- 960
ward12.epi.full[,"do.sat"] <- o2.at.sat.base(ward12.epi.full[,"wtr"], baro=ward12.epi.full[,"baro"]) # caulate oxygen concentration at saturation



# calculate k.read
ward12.sw <- par.to.sw.base(ward12.epi.full[,"irr"])
ward12.lw <- calc.lw.net.base(
	dateTime=ward12.epi.full[,"datetime"], 
	sw=ward12.sw, Ts=ward12.epi.full[,"wtr"], 
	lat=46.28, 
	atm.press=ward12.epi.full[,"baro"], 
	airT=ward12.epi.full[,"airTemp"],
	RH=ward12.epi.full[,"RH"]
	)
	
miss.read <- is.na(ward12.epi.full[,"wtr"]) | is.na(ward12.epi.full[,"kd"]) | is.na(ward12.epi.full[,"baro"]) | is.na(ward12.epi.full[,"airTemp"]) | is.na(ward12.epi.full[,"RH"]) | is.na(ward12.epi.full[,"z.mix"])

ward12.epi.full[!miss.read,"k600.read"] <- k.read.base(
	wnd.z=10,
	Kd=ward12.epi.full[!miss.read,"kd"],
	lat=46.28,
	lake.area=19000,
	atm.pres=ward12.epi.full[!miss.read,"baro"],
	dateTime=ward12.epi.full[!miss.read,"datetime"],
	Ts=ward12.epi.full[!miss.read,"wtr"],
	z.aml=ward12.epi.full[!miss.read,"z.mix"],
	airT=ward12.epi.full[!miss.read,"airTemp"],
	wnd=ward12.epi.full[!miss.read,"wnd"],
	RH=ward12.epi.full[!miss.read,"RH"],
	sw=ward12.sw[!miss.read],
	lwnet=ward12.lw[!miss.read]
	)
	
miss.read2 <- is.na(ward12.epi.full[,"k600.read"])

# dev.new(width=3.5, height=5)
# par(mfrow=c(2,1), mar=c(2.25, 2.25, 0.5, 0.5), mgp=c(1.25, 0.35, 0), tcl=-0.25, ps=10, family="Times")
# plot(ward12.epi.full[!miss.read2,"datetime"], ward12.epi.full[!miss.read2,"k600.read"], type="l", xlim=range(ward12.epi.full[,"datetime"]), xlab="", ylab=bquote(k600~~(m~d^-1)))
# lines(ward12.epi.full[,"datetime"], ward12.epi.full[,"k600.cole"], type="l", col="red")
# legend("topright", legend=c("k.cole", "k.read"), lty=1, lwd=1.5, col=c("red","black"))
# plot(log(ward12.epi.full[!miss.read2,"k600.cole"]), ward12.epi.full[!miss.read2,"k600.read"], xlab=bquote(log[e]~k600.cole~~(m~d^-1)), ylab=bquote(k600.read~~(m~d^-1)))

ward12.log.cole <- log(ward12.epi.full[,"k600.cole"])
ward12.lm.read.cole <- lm(ward12.epi.full[,"k600.read"]~ward12.log.cole)
ward12.read.from.cole <- predict(ward12.lm.read.cole, newdata=data.frame(ward12.log.cole=ward12.log.cole[miss.read2]))
ward12.epi.full[miss.read2,"k600.read"] <- ward12.read.from.cole
ward12.epi.full[,"kgas.read"] <- LakeMetabolizer:::k600.2.kGAS.base(ward12.epi.full[,"k600.read"], ward12.epi.full[,"wtr"], gas="O2")
w12.noK <- ward12.epi.full[,"z.mix"] < 0.7 & !is.na(ward12.epi.full[,"z.mix"])
ward12.epi.full[w12.noK,"kgas.cole"] <- 0L
ward12.epi.full[w12.noK,"kgas.read"] <- 0L








w12.chl.bad <- ward12.epi.full[,"Chl_conc"] > 100 | is.na(ward12.epi.full[,"Chl_conc"])
ward12.epi.full[w12.chl.bad,"Chl_conc"] <- NA

ward12.epi <- ward12.epi.full[,c("datetime", "do.obs", "do.sat", "kgas.cole", "kgas.read", "z.mix",  "irr", "wtr")] # subset for use with LakeMetabolizer



# ======================
# = Read Ward Meta 2012 =
# ======================
extra.meta12 <- list(
	"B14APR12.CDF"=data.frame("top"=2.25, "bot"=2.5, "z1perc"=2.5),
	"B03MAY12.txt"= data.frame("top"=2.25, "bot"=2, "z1perc"=2), # this is when the phytoflash went missing lol
	"B03JUN12.txt"= data.frame("top"=0.75, "bot"=1.25, "z1perc"=2),
	"B23JUN12.txt"= data.frame("top"=0.75, "bot"=1.25, "z1perc"=2),
	"B14JUL12.txt"= data.frame("top"=0.75, "bot"=1.25, "z1perc"=2),
	"B06AUG12.CDF"= data.frame("top"=0.75, "bot"=1.25, "z1perc"=2)
	)

for(i in 1:nrow(DeepKey2012)){
	tName <- as.character(DeepKey2012[i,"Name"])
	tFormat <- as.character(DeepKey2012[i,"Format"])
	
	if(tFormat=="CDF"){
		if(tName=="B14APR12.CDF"){
			ReadClassCDF <- c("character", "character", rep("numeric", 8))
			tDat <- read.table(paste(dat.path, tName, sep=""), sep=",", header=TRUE, skip=1, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU")
			tDat[,"do.obs"] <- Sat2Conc(tDat[,"DOsat"], tDat[,"wtr"])
		}
		
		if(all(tName!="B14APR12.CDF" & DeepKey2012[,"Depth"]=="B")){
			ReadClassCDF <- c("character", "character", rep("numeric", 9))
			tDat <- read.table(paste(dat.path, tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "do.obs")
		}
	
	# reformat datetime for .CDF files
	tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%m/%d/%y %T", tz="GMT")
	}
	
	
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste(dat.path, tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt)
		names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "do.obs", "Chl_conc", "Chl_RFU", "Battery")
		tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%Y/%m/%d %T", tz="GMT")
	}
	
	
	
	tDat <- data.frame(tDat[,c("datetime", "wtr", "do.obs", "DOsat", "SpCond", "pH", "Chl_conc")], extra.meta12[[tName]])
	if(i == 1){
		ward12.meta0 <- tDat
	}else{
		ward12.meta0 <- rbind(ward12.meta0, tDat)
	}
}



# =============================
# = Format and merge Ward Meta =
# =============================
ward12.meta0[,"datetime"] <- LakeMetabolizer:::round.time(ward12.meta0[,"datetime"], "5 minutes")
ward12.meta.full0 <- merge(ward12.meta0, ward12.zmix, all.x=TRUE)
# ward12.meta.full <- merge(ward12.meta.full0, irr_wnd_2012, all.x=TRUE)
ward12.meta.full <- merge(ward12.meta.full0, weather.full.2012, all.x=TRUE)

# wind is already scaled above where UNDERC and peter were being combined

ward12.meta.full[,"k600.cole"] <- k.cole.base(ward12.meta.full[,"wnd"]) # calculate k600 using Cole & Caraco method
# ward12.meta.full[,"kgas.cole"] <- k600.2.kGAS.base(ward12.meta.full[,"k600.cole"], ward12.meta.full[,"wtr"], gas="O2")
ward12.meta.full[,"k.gas"] <- 0L

ward12.meta.full[,"doy"] <- LakeMetabolizer:::date2doy(ward12.meta0[,"datetime"])


# Smooth metalimnetic temperature time series
ward12.meta.full[,"watts"] <- watts.in(ward12.meta.full[,"top"], ward12.meta.full[,"bot"], ward12.meta.full[,"irr"], ward12.meta.full[,"z1perc"])
ward12.meta.full[,"roundDoy"] <- trunc(ward12.meta.full[,"doy"])
w12m.tsmooth.func <- function(x){ # write this function to handle the short day (42nd day, maybe more after)
	if(any(is.na(x[,c("wtr","watts")])) | length(x[,1])<2){ # note the logic after | ... fails if only 1 obs
		cbind(x, "smooth.wtr"=NA)
	}else{
		cbind(x, "smooth.wtr"=temp.kalman(x[,"wtr"], x[,"watts"], ampH=1))
	}
}
ward12.meta.wtrSmooth <- ddply(ward12.meta.full, .variables="roundDoy", .fun=w12m.tsmooth.func) # smooth the temperature time series
ward12.meta.full[,"smooth.wtr"] <- ward12.meta.wtrSmooth[,"smooth.wtr"] # store the smoothed ts in original data.frame


ward12.meta.full[,"do.sat"] <- o2.at.sat.base(ward12.meta.full[,"wtr"], baro=960) # caulate oxygen concentration at saturation

# w12.chl.bad <- ward12.meta.full[,"Chl_conc"] > 100
# ward12.meta.full[w12.chl.bad,"Chl_conc"] <- NA

ward12.meta <- ward12.meta.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr")] # subset for use with LakeMetabolizer


# ==============================
# = Save 2012 Ward Sensor Data =
# ==============================
save(ward12.therm, ward12.epi.full, ward12.epi, ward12.meta.full, ward12.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")

write.table(ward12.therm, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward12.therm.txt", row.names=FALSE)
write.table(ward12.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward12.epi.full.txt", row.names=FALSE)
write.table(ward12.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward12.epi.txt", row.names=FALSE)
write.table(ward12.meta.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward12.meta.full.txt", row.names=FALSE)
write.table(ward12.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward12.meta.txt", row.names=FALSE)

# dev.new();plot(ward12.epi0[,"datetime"], ward12.epi0[,"DOsat"], type="o", ylim=c(0, 120)); lines(ward12.meta0[,"datetime"], ward12.meta0[,"DOsat"], type="o", col="red")
