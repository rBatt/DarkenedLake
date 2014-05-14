
rm(list=ls())
graphics.off()
library("zoo")
library("plyr")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/PAR_DO_PhaseShift/")
source("Alme.R")
source("ByeShort.R")
source("Chunks.R")
source("Light.R")
source("Metabolism_LM_v2.2.R")

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/WardSensorData2012")
UNDERC_Weather <- read.csv("UNDERC_Weather_2012.csv")
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

Peter_PAR_Wind <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/Peter_PAR_Wind_2012.csv")
Peter_PAR_Wind[,"DoY"] <- as.numeric(format.Date(as.POSIXct(Peter_PAR_Wind[,1], format="%Y-%m-%d"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(Peter_PAR_Wind[,1], Peter_PAR_Wind[,2]), format="%Y-%m-%d %H:%M"), time2=as.POSIXct(paste(Peter_PAR_Wind[,1], "00:00:00"), format="%Y-%m-%d %H:%M:%S"), units="days"))
Peter_PAR_Wind[,"Year"] <- 2012
Peter_PAR_Wind <- Peter_PAR_Wind[,c("Year", "DoY", "PAR", "WindSpeed")]
names(Peter_PAR_Wind) <- c("Year", "DoY", "PAR", "Wind")

CombinePeterUNDERC <- Alme(X=Peter_PAR_Wind, Y=UNDERC_PAR_WIND, CheckMissDups=TRUE)
PAR_Conv <- summary(lm(I(PAR0+0.01)~I(PAR-1.19) -1, data=CombinePeterUNDERC))$coef[1]
UNDERC_PAR_WIND[,"PAR0"] <- UNDERC_PAR_WIND[,"PAR0"]/PAR_Conv
names(UNDERC_PAR_WIND) <- c("Year", "DoY", "PAR", "Wind")
Peter_PAR_Wind <- ByeShort(Peter_PAR_Wind)
UNDERC_PAR_WIND <- ByeShort(UNDERC_PAR_WIND)
UNDERC_2add <- which(UNDERC_PAR_WIND[,"DoY"] < min(Peter_PAR_Wind[,"DoY"]))
PAR_Wind <- rbind(UNDERC_PAR_WIND[UNDERC_2add,], Peter_PAR_Wind)


##****WORKING ON READING IN ZMIX*******
Therms <- data.frame("Year"=2012, read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Data/CompiledWardTherms2012.txt", header=TRUE))
ThermDepths <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3)
zmix <- function(temp, depth){
        tempDiff <- temp[-length(temp)] - temp[-1]
        depthDiff <- depth[-1] - depth[-length(depth)]
        ratio <- tempDiff / depthDiff
		zhat <- depth[which(ratio >= 2)][1]
		Z <- ifelse(zhat==0, 0.25, zhat)
        return(depth[which(ratio >= 2)][1])
}
Zmixs <- apply(Therms[,-c(1,2)], MARGIN=1, FUN=zmix, depth=ThermDepths)
Therms <- data.frame(Therms, "Zmix"=Zmixs)
# dim(Therms[which(is.element(Therms[,"DoY"], PAR_Wind[,"DoY"])),])
# length(which(PAR_Wind[,"DoY"] > 105 & PAR_Wind[,"DoY"] < 226))
PAR_Wind_Therms <- Alme(X=PAR_Wind, Y=Therms, Freq_High_Low="Low", AllX=TRUE)
# plot(1:length(PAR_Wind_Therms[,"Zmix"]), PAR_Wind_Therms[,"Zmix"])
plot(1:length(PAR_Wind_Therms[,"z0.00"]), PAR_Wind_Therms[,"z0.00"])
# abline(v=which(abs(diff(PAR_Wind_Therms[,5]))>1.5), col="blue", lwd=2)
# min(which(!is.na(PAR_Wind_Therms[,"z0.00"])))
# abline(v=c(17051, 19384), col="red")
# abline(v=c(12495, 12519), col="green")
# PAR_Wind_Therms[c(17051:19384, 12495:12519), 5:16] <- NA
abline(v=c(17004, 19330), col="red")
abline(v=c(12464, 12488), col="green")
PAR_Wind_Therms[c(17004:19330, 12454:12498), 5:16] <- NA
plot(1:length(PAR_Wind_Therms[,"z0.00"]), PAR_Wind_Therms[,"z0.00"])
abline(v=c(17004, 19330), col="red")
abline(v=c(12454, 12498), col="green")
# abline(h=c(27, 33), col="red")
PAR_Wind_Therms[, "Zmix"] <- approx(c(PAR_Wind_Therms[, "DoY"], 152, 159, 229, 236, 241), c(PAR_Wind_Therms[, "Zmix"], 1.5, 0.5, 1.5, 1, 0.25), PAR_Wind_Therms[, "DoY"])$y

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/")
Key2012 <- read.csv("WardSondes2012_FileKey.csv")
ShalKey2012 <- subset(Key2012, Depth=="A")


for(i in 1:nrow(ShalKey2012)){
	tName <- ShalKey2012[i,"Name"]
	tFormat <- ShalKey2012[i,"Format"]
	
	if(tFormat=="CDF"){
		ReadClassCDF <- c("character", "character", rep("numeric", 9))
		tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
		names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "Battery")
		tDoY <- as.numeric(format.Date(as.POSIXct(as.character(tDat[,1]), format="%m/%d/%y"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%m/%d/%y %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%m/%d/%y %H:%M:%S"), units="days"))
	}
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
		names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "DOconc", "Chl_conc", "Chl_RFU", "Battery")
		tDoY <- as.numeric(format.Date(as.POSIXct(tDat[,1], format="%Y/%m/%d"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%Y/%m/%d %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%Y/%m/%d %H:%M:%S"), units="days"))
		tDat <- tDat[,-9]
	}
	
	if(i == 1){
		DataSonde <- data.frame("Year"=2012, "DoY"=tDoY, tDat[,-c(1:2)])
	}
	
	if(i !=1){
		DataSonde <- rbind(DataSonde, data.frame("Year"=2012, "DoY"=tDoY, tDat[,-c(1:2)]))
	}
	print(i)
	print(dim(tDat))
	print(dim(DataSonde))
	
	
}
# DataSonde <- as.data.frame(DataSonde)
# PAR_Wind_Therms <- as.data.frame(PAR_Wind_Therms)
sDataSonde <- ByeShort(DataSonde)
# Test <- Alme(X=sDataSonde, Y=PAR_Wind_Therms, AllX=TRUE)
# head(Test)
# length(unique(trunc(Test[,"DoY"])))
#missed 10 days from the batteries dying, and would miss 1 day for each deployment (calibration/first deployment/last day of season when it was removed= 5+1+1). 1 of the missing days was the same for the battery dying/ calibration (i pulled out a dead sonde). That's 16 missing days. First day of measurement was 105, last was 139, so that is a total of 135 days the sonde could have been in the lake (139 - 105 + 1).  My final data set has 119 days.  119 + 16 = 135.  OK, I didn't lose any days of data due to stupid/sloppy programming.



Data0 <- Alme(X=sDataSonde, Y=PAR_Wind_Therms)

Data <- Data0[,c("Year", "DoY", "DOsat", "Temp", "PAR", "Zmix", "Wind")]

Metabolism_LM(Data)
