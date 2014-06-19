
# ===================================
# = Read in Peter Weather from 2010 =
# ===================================
PeterWeather00 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/PeterWeather.csv", header=TRUE)[-(96344:96345),]
PeterWeather0 <- PeterWeather00[PeterWeather00[,"Year"]==2010L,]
PeterWeather0 <- PeterWeather0[,c("Date", "Time", "Year", "PAR","WindSpd")]

PeterWeather0[,"datetime"] <- paste(PeterWeather0[,"Date"], PeterWeather0[,"Time"])
PeterWeather0[,"datetime"] <- gsub("^(?=[0-9]/)", "0", PeterWeather0[,"datetime"], perl=TRUE)
PeterWeather0[,"datetime"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", PeterWeather0[,"datetime"], perl=TRUE)
PeterWeather0[,"datetime"] <- as.POSIXct(PeterWeather0[,"datetime"], format="%m/%d/%y %I:%M:%S %p")
irr_wnd_2010 <- PeterWeather0[,c("datetime", "PAR", "WindSpd")]
names(irr_wnd_2010) <- c("datetime", "irr", "wnd")

save(irr_wnd_2010, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")




# ======================================
# = Peter and UNDERC Weather from 2012 =
# ======================================
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
Peter_PAR_Wind0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/Peter_PAR_Wind_2012.csv")
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

CombinePeterUNDERC[!UNDERC_2add.wnd, "wnd"] <- scale.exp.wind.base(CombinePeterUNDERC[!UNDERC_2add.wnd, "wnd"], 2)

CombinePeterUNDERC[UNDERC_2add.wnd, "wnd"] <- CombinePeterUNDERC[UNDERC_2add.wnd, "Wind"]
CombinePeterUNDERC[UNDERC_2add.irr, "irr"] <- CombinePeterUNDERC[UNDERC_2add.irr, "PAR0"]
irr_wnd_2012 <- CombinePeterUNDERC[,c("datetime", "irr", "wnd")]

save(irr_wnd_2012, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2012.RData")
