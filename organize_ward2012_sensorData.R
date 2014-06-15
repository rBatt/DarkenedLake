
# Weather, sonde, and thermistor data from Ward Lake in 2012
library("plyr")

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
UNDERC_2add <- which(UNDERC_PAR_WIND[,"DoY"] < min(Peter_PAR_Wind[,"DoY"]))
PAR_Wind <- rbind(UNDERC_PAR_WIND[UNDERC_2add,], Peter_PAR_Wind)

# New 15-June-2014
PAR_Wind[,"datetime"] <- as.POSIXct(PAR_Wind[,"DoY"]*24*60*60, origin="2012-01-01", tz="GMT")




