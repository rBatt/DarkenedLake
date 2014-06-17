
# Weather, sonde, and thermistor data from Ward Lake in 2012
library("plyr")
library("rLakeAnalyzer")

detach(package:LakeMetabolizer, unload=TRUE)
install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Alme.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/ByeShort.R")

# ===============
# = Manual Zmix =
# ===============
w12.manZ.date0 <- c(105, 110, 117, 124, 131, 138, 145, 152, 159, 166, 173, 180, 187, 194, 201, 208, 215, 222, 229, 236, 241)
w12.manZ.date <- as.POSIXct((w12.manZ.date0)*24*60*60, origin="2011-12-31", tz="GMT") # confirm: as.POSIXct("105-2012", format="%j-%Y")
w12.manZ00 <- c(1.0,2.0,2.0,0.5,0.5, 1.0, 1.0, 1.5, 0.25, 1.0, 1.0, 0.5, 0.5, 1, 1, 1, 1, 1.5, 1.5, 1, 0.25)
w12.manZ.interpDate <- seq(from=as.POSIXct("2012-04-14 00:00:00", tz="GMT"), to=as.POSIXct("2012-08-28 23:55:00", tz="GMT"), by=60*5)
w12.manZ0 <- approx(x=w12.manZ.date, y=w12.manZ00, xout=w12.manZ.interpDate, method="constant", f=0, rule=2)
w12.manZ <- data.frame("datetime"=w12.manZ0$x, "manZ"=w12.manZ0$y)


# ===========
# = Weather =
# ===========

# UNDERC Weather (some Peter weather is missing)
UNDERC_Weather <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSensorData2012/UNDERC_Weather_2012.csv")
names(UNDERC_Weather) <- c("Year", "DoY", "Time", "Wind", "PAR")
UNDERC_Weather[,"Time"] <- UNDERC_Weather[,"Time"]-100
UNDERC_Weather[,"Frac"] <- UNDERC_Weather[,"Time"]/2400
HourChar <- as.character(ifelse(UNDERC_Weather[,"Time"]/100==24, "0",UNDERC_Weather[,"Time"]/100))
WhichHourSingleDigit <- which(nchar(HourChar)==1)
HourChar[WhichHourSingleDigit] <- paste("0",HourChar[WhichHourSingleDigit],sep="")
UNDERC_Weather[,"Time"] <- paste(HourChar, "00", sep=":")
UNDERC_Weather[,"DoY"] <- UNDERC_Weather[,"DoY"] + UNDERC_Weather[,"Frac"]
UNDERC_Interp_xout <- seq(min(UNDERC_Weather[,"DoY"]), max(UNDERC_Weather[,"DoY"]), by=(1/288))
UNDERC_PAR_WIND <- data.frame("Year"=2012, "DoY"=UNDERC_Interp_xout, "PAR0"=approx(UNDERC_Weather[,"DoY"], UNDERC_Weather[,"PAR"], xout=UNDERC_Interp_xout)$y, "Wind"=approx(UNDERC_Weather[,"DoY"], UNDERC_Weather[,"Wind"], xout=UNDERC_Interp_xout)$y)

# New 15-June-2014: Investigating issue with first part of underc wind speed never getting below
lowest.underc.wind <- min(UNDERC_PAR_WIND[,"Wind"])
UNDERC_PAR_WIND[,"Wind"] <- UNDERC_PAR_WIND[,"Wind"] - min(UNDERC_PAR_WIND[,"Wind"])
UNDERC_PAR_WIND[,"datetime"] <- as.POSIXct(UNDERC_PAR_WIND[,"DoY"]*24*60*60, origin="2011-12-31", tz="GMT")


# Peter Weather
Peter_PAR_Wind0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/Peter_PAR_Wind_2012.csv")
Peter_PAR_Wind0[,"datetime"] <- as.POSIXct(paste(Peter_PAR_Wind0[,"Date"], Peter_PAR_Wind0[,"Time"]), tz="GMT")
# Peter_PAR_Wind0[,"doy"] <- LakeMetabolizer:::date2doy(Peter_PAR_Wind0[,"datetime"])
# Peter_PAR_Wind[,"year"] <- 2012
Peter_PAR_Wind <- Peter_PAR_Wind0[,c("datetime", "PAR", "WindSpeed")]
names(Peter_PAR_Wind) <- c("datetime", "irr", "wnd")

# Combine Peter and UNDERC weather
CombinePeterUNDERC <- merge(UNDERC_PAR_WIND, Peter_PAR_Wind, all=TRUE)


# CombinePeterUNDERC <- Alme(X=Peter_PAR_Wind, Y=UNDERC_PAR_WIND, CheckMissDups=TRUE)
# New 15-June-2014: just checking for the hysteresis in the relationship between UNDERC and Peter PAR (it's there, o well)
# ind <- (CombinePeterUNDERC[,"DoY"] - trunc(CombinePeterUNDERC[,"DoY"])) < 0.5
# plot(CombinePeterUNDERC[,"PAR0"], CombinePeterUNDERC[,"PAR"], col=ind+1)
# summary(lm(I(PAR0+0.01)~I(PAR-1.19) + ind, data=CombinePeterUNDERC))
PAR_Conv <- summary(lm(I(PAR0+0.01)~I(irr-1.19) -1, data=CombinePeterUNDERC))$coef[1]
UNDERC_PAR_WIND[,"PAR0"] <- UNDERC_PAR_WIND[,"PAR0"]/PAR_Conv
names(UNDERC_PAR_WIND) <- c("Year", "DoY", "PAR", "Wind", "datetime")
# Peter_PAR_Wind <- ByeShort(Peter_PAR_Wind)
# UNDERC_PAR_WIND <- ByeShort(UNDERC_PAR_WIND)
# UNDERC_2add <- UNDERC_PAR_WIND[,"datetime"] < min(Peter_PAR_Wind[,"datetime"])
# PAR_Wind0 <- rbind(UNDERC_PAR_WIND[UNDERC_2add,], Peter_PAR_Wind)

UNDERC_2add.wnd <- is.na(CombinePeterUNDERC[,"wnd"]) & !is.na(CombinePeterUNDERC[,"Wind"])
UNDERC_2add.irr <- is.na(CombinePeterUNDERC[,"irr"]) & !is.na(CombinePeterUNDERC[,"PAR0"])

CombinePeterUNDERC[UNDERC_2add.wnd, "wnd"] <- CombinePeterUNDERC[UNDERC_2add.wnd, "Wind"]
CombinePeterUNDERC[UNDERC_2add.irr, "irr"] <- CombinePeterUNDERC[UNDERC_2add.irr, "PAR0"]
PAR_Wind <- CombinePeterUNDERC[,c("datetime", "irr", "wnd")]

# New 15-June-2014
# PAR_Wind0[,"datetime"] <- as.POSIXct(PAR_Wind0[,"DoY"]*24*60*60, origin="2012-01-01", tz="GMT")
# PAR_Wind <- PAR_Wind0[,c("datetime", "PAR", "Wind")]
# names(PAR_Wind) <- c("datetime", "irr", "wnd")


# =======================
# = Temp Profile/ z.mix =
# =======================
Therms0 <- data.frame("Year"=2012, read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/CompiledWardTherms2012_v2.txt", header=TRUE))
ThermDepths <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3)
t.wtr.names <- paste("wtr", as.character(ThermDepths), sep="_")
names(Therms0) <- c("year", "doy", t.wtr.names)
Therms0[,"datetime"] <- as.POSIXct(Therms0[,"doy"]*24*60*60, origin="2012-01-01", tz="GMT")

Therms <- Therms0[,c("datetime",t.wtr.names)]
Therms[,"datetime"] <- LakeMetabolizer:::round.time(Therms[,"datetime"], "5 minutes")
Therms <- Therms[!duplicated(Therms[,"datetime"]),]

ward12.zmix000 <- ts.meta.depths(Therms)
ward12.zmix00 <- ward12.zmix000
ward12.zmix00_is0 <- ward12.zmix00[,"top"] < 0.25 & !is.na(ward12.zmix00[,"top"])
ward12.zmix00[ward12.zmix00_is0, "top"] <- 0.25

# just checking some things
# class(ward12.zmix[,"datetime"])
# plot(ward12.zmix[,"datetime"], ward12.zmix[,"top"], type="l")
# LakeMetabolizer:::addNAs(head(ward12.zmix, 288))
# blah <- LakeMetabolizer:::addNAs(ward12.zmix)

ward12.zmix0 <- ward12.zmix00[,c("datetime", "top")]
names(ward12.zmix0) <- c("datetime", "z.mix")

w12.bothZ <- merge(ward12.zmix0, w12.manZ, all=TRUE)
w12z.isNA <- is.na(w12.bothZ[,"z.mix"])
ward12.zmix <- w12.bothZ[,c("datetime", "z.mix")]
ward12.zmix[w12z.isNA, "z.mix"] <- w12.bothZ[w12z.isNA, "manZ"]



# ==============
# = Sonde data =
# ==============

dat.path <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
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
			tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=1, colClasses=ReadClassCDF)
			# head(read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=","))
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "Battery")
			tDat[,"do.obs"] <- Sat2Conc(tDat[,"DOsat"], tDat[,"wtr"])
		}
	
	# reformat datetime for .CDF files
	tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%m/%d/%y %T", tz="GMT")
	}
	
	
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
		# head(read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=","))
		names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "do.obs", "Chl_conc", "Chl_RFU", "Battery")
		tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%Y/%m/%d %T", tz="GMT")
		# tDat <- tDat[,-9]
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
ward12.epi.full <- merge(ward12.epi.full0, PAR_Wind, all.x=TRUE)

ward12.epi.full[!UNDERC_2add, "wnd"] <- scale.exp.wind.base(ward12.epi.full[!UNDERC_2add, "wnd"], 2) # scale wind speed to 10 m

# ward12.epi.k.cole <- k.cole.base(ward12.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.epi.full[,"k600.cole"] <- k.cole.base(ward12.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.epi.full[,"k.gas"] <- k600.2.kGAS.base(ward12.epi.full[,"k600.cole"], ward12.epi.full[,"wtr"], gas="O2")
w12.noK <- ward12.epi.full[,"z.mix"] < 0.7 & !is.na(ward12.epi.full[,"z.mix"])
ward12.epi.full[w12.noK,"k.gas"] <- 0L
# those k values seem high, so just checking (looks good)
# plot(ward12.epi.k.cole)
# par(new=TRUE); plot(ward12.epi.full[,"wnd"], type="l", col="red", xlab="", ylab="", xaxt="n", yaxt="n"); axis(side=4)
# abline(v=min(which(ward12.epi.full[,"wnd"] < 0.1)), col="blue") # note that I adjusted the underc wind speeds down so that the lowest was 0, not 0.2811882

ward12.epi.full[,"do.sat"] <- o2.at.sat.base(ward12.epi.full[,"wtr"], baro=960) # caulate oxygen concentration at saturation
w12.chl.bad <- ward12.epi.full[,"Chl_conc"] > 100
ward12.epi.full[w12.chl.bad,"Chl_conc"] <- NA

ward12.epi <- ward12.epi.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr")] # subset for use with LakeMetabolizer



# ======================
# = Read Ward Meta 2012 =
# ======================
for(i in 1:nrow(DeepKey2012)){
	tName <- as.character(DeepKey2012[i,"Name"])
	tFormat <- as.character(DeepKey2012[i,"Format"])
	
	if(tFormat=="CDF"){
		if(tName=="B14APR12.CDF"){
			ReadClassCDF <- c("character", "character", rep("numeric", 8))
			tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=1, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU")
			tDat[,"do.obs"] <- Sat2Conc(tDat[,"DOsat"], tDat[,"wtr"])
		}
		
		if(all(tName!="B14APR12.CDF" & DeepKey2012[,"Depth"]=="B")){
			ReadClassCDF <- c("character", "character", rep("numeric", 9))
			tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "do.obs")
		}
	
	# reformat datetime for .CDF files
	tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%m/%d/%y %T", tz="GMT")
	}
	
	
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
		names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "do.obs", "Chl_conc", "Chl_RFU", "Battery")
		tDat[,"datetime"] <- as.POSIXct(paste(as.character(tDat[,"Date"]), as.character(tDat[,"Time"])), format="%Y/%m/%d %T", tz="GMT")
		# tDat <- tDat[,-9]
	}
	
	
	tDat <- tDat[,c("datetime", "wtr", "do.obs", "DOsat", "SpCond", "pH", "Chl_conc")]
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
ward12.meta.full <- merge(ward12.meta.full0, PAR_Wind, all.x=TRUE)

ward12.meta.full[!UNDERC_2add, "wnd"] <- scale.exp.wind.base(ward12.meta.full[!UNDERC_2add, "wnd"], 2) # scale wind speed to 10 m

# ward12.meta.k.cole <- k.cole.base(ward12.meta.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.meta.full[,"k600.cole"] <- k.cole.base(ward12.meta.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.meta.full[,"k.gas"] <- k600.2.kGAS.base(ward12.meta.full[,"k600.cole"], ward12.meta.full[,"wtr"], gas="O2")
# w12.noK <- ward12.meta.full[,"z.mix"] < 0.7 & !is.na(ward12.meta.full[,"z.mix"])
# ward12.meta.full[w12.noK,"k.gas"] <- 0L
ward12.meta.full[,"k.gas"] <- 0L
# those k values seem high, so just checking (looks good)
# plot(ward12.meta.k.cole)
# par(new=TRUE); plot(ward12.meta.full[,"wnd"], type="l", col="red", xlab="", ylab="", xaxt="n", yaxt="n"); axis(side=4)
# abline(v=min(which(ward12.meta.full[,"wnd"] < 0.1)), col="blue") # note that I adjusted the underc wind speeds down so that the lowest was 0, not 0.2811882

ward12.meta.full[,"do.sat"] <- o2.at.sat.base(ward12.meta.full[,"wtr"], baro=960) # caulate oxygen concentration at saturation
w12.chl.bad <- ward12.meta.full[,"Chl_conc"] > 100
ward12.meta.full[w12.chl.bad,"Chl_conc"] <- NA

ward12.meta <- ward12.meta.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr")] # subset for use with LakeMetabolizer



# dev.new();plot(ward12.epi0[,"datetime"], ward12.epi0[,"DOsat"], type="o", ylim=c(0, 120)); lines(ward12.meta0[,"datetime"], ward12.meta0[,"DOsat"], type="o", col="red")
