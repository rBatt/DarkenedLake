
# Weather, sonde, and thermistor data from Ward Lake in 2012
library("plyr")
library("rLakeAnalyzer")

detach(package:LakeMetabolizer, unload=TRUE)
install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Alme.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/ByeShort.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Chunks.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Light.R")

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


# Peter Weather
Peter_PAR_Wind <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/Peter_PAR_Wind_2012.csv")
Peter_PAR_Wind[,"DoY"] <- as.numeric(format.Date(as.POSIXct(Peter_PAR_Wind[,1], format="%Y-%m-%d"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(Peter_PAR_Wind[,1], Peter_PAR_Wind[,2]), format="%Y-%m-%d %H:%M"), time2=as.POSIXct(paste(Peter_PAR_Wind[,1], "00:00:00"), format="%Y-%m-%d %H:%M:%S"), units="days"))
Peter_PAR_Wind[,"Year"] <- 2012
Peter_PAR_Wind <- Peter_PAR_Wind[,c("Year", "DoY", "PAR", "WindSpeed")]
names(Peter_PAR_Wind) <- c("Year", "DoY", "PAR", "Wind")

CombinePeterUNDERC <- Alme(X=Peter_PAR_Wind, Y=UNDERC_PAR_WIND, CheckMissDups=TRUE)
# New 15-June-2014: just checking for the hysteresis in the relationship between UNDERC and Peter PAR (it's there, o well)
# ind <- (CombinePeterUNDERC[,"DoY"] - trunc(CombinePeterUNDERC[,"DoY"])) < 0.5
# plot(CombinePeterUNDERC[,"PAR0"], CombinePeterUNDERC[,"PAR"], col=ind+1)
# summary(lm(I(PAR0+0.01)~I(PAR-1.19) + ind, data=CombinePeterUNDERC))
PAR_Conv <- summary(lm(I(PAR0+0.01)~I(PAR-1.19) -1, data=CombinePeterUNDERC))$coef[1]
UNDERC_PAR_WIND[,"PAR0"] <- UNDERC_PAR_WIND[,"PAR0"]/PAR_Conv
names(UNDERC_PAR_WIND) <- c("Year", "DoY", "PAR", "Wind")
Peter_PAR_Wind <- ByeShort(Peter_PAR_Wind)
UNDERC_PAR_WIND <- ByeShort(UNDERC_PAR_WIND)
UNDERC_2add <- UNDERC_PAR_WIND[,"DoY"] < min(Peter_PAR_Wind[,"DoY"])
PAR_Wind0 <- rbind(UNDERC_PAR_WIND[UNDERC_2add,], Peter_PAR_Wind)

# New 15-June-2014
PAR_Wind0[,"datetime"] <- as.POSIXct(PAR_Wind0[,"DoY"]*24*60*60, origin="2012-01-01", tz="GMT")
PAR_Wind <- PAR_Wind0[,c("datetime", "PAR", "Wind")]
names(PAR_Wind) <- c("datetime", "irr", "wnd")


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

ward12.zmix0 <- ts.meta.depths(Therms)
w12.zmix.is0 <- ward12.zmix0[,"top"] < 0.25 & !is.na(ward12.zmix0[,"top"])
ward12.zmix0[w12.zmix.is0, "top"] <- 0.25

# just checking some things
# class(ward12.zmix[,"datetime"])
# plot(ward12.zmix[,"datetime"], ward12.zmix[,"top"], type="l")
# LakeMetabolizer:::addNAs(head(ward12.zmix, 288))
# blah <- LakeMetabolizer:::addNAs(ward12.zmix)

ward12.zmix <- ward12.zmix0[,c("datetime", "top")]
names(ward12.zmix) <- c("datetime", "z.mix")

# ==============
# = Sonde data =
# ==============

dat.path <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
Key2012 <- read.csv(paste(dat.path, "WardSondes2012_FileKey.csv", sep=""))
ShalKey2012 <- Key2012[Key2012[,"Depth"]=="B",]
DeepKey2012 <- Key2012[Key2012[,"Depth"]=="A",]

Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}

# ======================
# = Read Ward Epi 2012 =
# ======================
for(i in 1:nrow(ShalKey2012)){
	tName <- as.character(ShalKey2012[i,"Name"])
	tFormat <- as.character(ShalKey2012[i,"Format"])
	
	if(tFormat=="CDF"){
		if(tName=="B14APR12.CDF"){
			ReadClassCDF <- c("character", "character", rep("numeric", 8))
			tDat <- read.table(paste(dat.path, "WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=1, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU")
			tDat[,"do.obs"] <- Sat2Conc(tDat[,"DOsat"], tDat[,"wtr"])
		}
		
		if(all(tName!="B14APR12.CDF" & ShalKey2012[,"Depth"]=="B")){
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
		ward12.epi0 <- tDat
	}else{
		ward12.epi0 <- rbind(ward12.epi0, tDat)
	}
}


# scaletemp <- ddply(data.frame(ward12.epi0, "week"=format.Date(ward12.epi0[,"datetime"], format="%W")), "week", function(x){abs(scale(x[,"wtr"]))})[,2]
# plot(ward12.epi0[,"wtr"], col=(scaletemp>4)+1)



# =======================
# = Read Ward Meta 2012 =
# =======================
for(i in 1:nrow(DeepKey2012)){
	tName <- as.character(DeepKey2012[i,"Name"])
	tFormat <- as.character(DeepKey2012[i,"Format"])
	
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
		ward12.meta0 <- tDat
	}else{
		ward12.meta0 <- rbind(ward12.epi0, tDat)
	}
}



# =============================
# = Format and merge Ward Epi =
# =============================
ward12.epi0[,"datetime"] <- LakeMetabolizer:::round.time(ward12.epi0[,"datetime"], "5 minutes")
ward12.epi.full0 <- merge(ward12.epi0, ward12.zmix, all.x=TRUE)
ward12.epi.full <- merge(ward12.epi.full0, PAR_Wind, all.x=TRUE)

ward12.epi.full[!UNDERC_2add, "wnd"] <- scale.exp.wind.base(ward12.epi.full[!UNDERC_2add, "wnd"], 2) # scale wind speed to 10 m

ward12.epi.k.cole <- k.cole.base(ward12.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward12.epi.full[,"k.gas"] <- k600.2.kGAS.base(ward12.epi.k.cole, ward12.epi.full[,"wtr"], gas="O2")
w12.noK <- ward12.epi.full[,"z.mix"] < 0.5 & !is.na(ward12.epi.full[,"z.mix"])
ward12.epi.full[w12.noK,"k.gas"] <- 0L
# those k values seem high, so just checking (looks good)
# plot(ward12.epi.k.cole)
# par(new=TRUE); plot(ward12.epi.full[,"wnd"], type="l", col="red", xlab="", ylab="", xaxt="n", yaxt="n"); axis(side=4)
# abline(v=min(which(ward12.epi.full[,"wnd"] < 0.1)), col="blue") # note that I adjusted the underc wind speeds down so that the lowest was 0, not 0.2811882

ward12.epi.full[,"do.sat"] <- o2.at.sat.base(ward12.epi.full[,"wtr"], baro=960) # caulate oxygen concentration at saturation

ward12.epi <- ward12.epi.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr")] # subset for use with LakeMetabolizer



# =======================
# = Test epi metabolism =
# =======================

# MLE
ward12.epi.mle <- metab(ward12.epi, "mle")
ward12.epi.mle.res <- ward12.epi.mle[,c("doy","GPP","R", "NEP")]
ward12.epi.mle.res <- merge(data.frame("doy"=138:239), ward12.epi.mle.res, all=TRUE)

# Kalman
ward12.epi.kal <- metab(ward12.epi, "kalman")
ward12.epi.kal.res <- ward12.epi.kal[,c("doy","GPP","R", "NEP")]
ward12.epi.kal.res <- merge(data.frame("doy"=138:239), ward12.epi.kal.res, all=TRUE)

# Bayesian
ward12.epi.bay <- metab(ward12.epi, "bayesian")
ward12.epi.bay.res <- ward12.epi.bay[,c("doy","GPP","R", "NEP")]
ward12.epi.bay.res <- merge(data.frame("doy"=138:239), ward12.epi.bay.res, all=TRUE)

# Reformat data for bookkeeping, do BK

# ward12.epi3 <- ward12.epi2
# bk.irr <- as.integer(LakeMetabolizer:::is.day(ward12.epi3[,"datetime"], 46.28))
# ward12.epi3[,"irr"] <- bk.irr
# ward12.epi.bk <- metab(ward12.epi3, "bookkeep")
ward12.epi.bk <- metab(ward12.epi, "bookkeep", lake.lat=46.28)
ward12.epi.bk.res <- ward12.epi.bk[,c("doy","GPP","R", "NEP")]
ward12.epi.bk.res <- merge(data.frame("doy"=138:239), ward12.epi.bk.res, all=TRUE)

# OLS
ward12.epi.ols <- metab(ward12.epi, "ols")
ward12.epi.ols.res <- ward12.epi.ols[,c("doy","GPP","R", "NEP")]
ward12.epi.ols.res <- merge(data.frame("doy"=138:239), ward12.epi.ols.res, all=TRUE)

# merge all results
w12e.r1 <- merge(cbind("method"="mle",ward12.epi.mle.res), cbind("method"="kalman",ward12.epi.kal.res), all=TRUE)
w12e.r2 <- merge(w12e.r1, cbind("method"="bayes", ward12.epi.bay.res), all=TRUE)
w12e.r3 <- merge(w12e.r2, cbind("method"="bookkeep", ward12.epi.bk.res), all=TRUE)
w12e.r4 <- merge(w12e.r3, cbind("method"="ols", ward12.epi.ols.res), all=TRUE)





# ==========
# = lower? =
# ==========
good10 <- ward10.epi.kal[,"GPP"] > 0 & ward10.epi.kal[,"R"] < 0
colMeans(ward10.epi.kal[good10,])

good12 <- ward12.epi.kal[,"GPP"] > 0 & ward12.epi.kal[,"R"] < 0
colMeans(ward12.epi.kal[good12,])


# =====================================
# = test for z.mix deepening at night =
# =====================================
w12.day <- is.day(ward12.epi[, "datetime"], 46.28)

w12.doy <- LakeMetabolizer:::date2doy(ward12.epi[,"datetime"])
w12.dZ <- ddply(cbind(doy=trunc(w12.doy), isDay=w12.day, ward12.epi), .variables=c("doy", "isDay"), function(x){data.frame(dZ=sum(diff(x[,"z.mix"]), na.rm=TRUE))})
w12.dZ.day <- w12.dZ[w12.dZ[,"isDay"],"dZ"]
w12.dZ.night <- w12.dZ[!w12.dZ[,"isDay"],"dZ"]
dev.new()
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))
plot(w12.dZ[!w12.dZ[,"isDay"],c("doy","dZ")], type="l")
plot(w12.dZ.night-w12.dZ.day, type="l")
abline(h=0, col="blue")



