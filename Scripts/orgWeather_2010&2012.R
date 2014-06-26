library(plyr)

# ======================================
# = Read in 2010 & 2012 UNDERC Weather =
# ======================================
UNDERC_Weather0 <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/undercWeather_2010&2012.txt", header=TRUE, sep="\t")
UNDERC_Weather0[,"datetime"] <- as.POSIXct(UNDERC_Weather0[,"datetime"], tz="GMT")

fill.uw <- function(UNDERC_Weather0){
	uw.time.range <- range(UNDERC_Weather0[,"datetime"])
	uw.datetime <- seq(uw.time.range[1], uw.time.range[2]+60*55, by=5*60)
	uw.airTemp <- approx(UNDERC_Weather0[,"datetime"], UNDERC_Weather0[,"airTemp"], xout=uw.datetime, rule=2)$y
	uw.RH <- approx(UNDERC_Weather0[,"datetime"], UNDERC_Weather0[,"RH"], xout=uw.datetime, rule=2)$y
	uw.wnd <- approx(UNDERC_Weather0[,"datetime"], UNDERC_Weather0[,"wnd"], xout=uw.datetime, rule=2)$y
	uw.irr <- approx(UNDERC_Weather0[,"datetime"], UNDERC_Weather0[,"irr"], xout=uw.datetime, rule=2)$y
	uw.baro <- approx(UNDERC_Weather0[,"datetime"], UNDERC_Weather0[,"baro"], xout=uw.datetime, rule=2)$y
	
	data.frame(datetime=uw.datetime, airTemp=uw.airTemp, RH=uw.RH, wnd=uw.wnd, irr=uw.irr, baro=uw.baro)
}

UNDERC_Weather <- ddply(UNDERC_Weather0, c("year", "doy"), fill.uw)


# ===================================
# = Read in Peter Weather from 2010 =
# ===================================
PeterWeather00 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/PeterWeather.csv", header=TRUE)[-(96344:96345),]
PeterWeather0 <- PeterWeather00[PeterWeather00[,"Year"]==2010L,]
PeterWeather0 <- PeterWeather0[,c("Date", "Time", "Year", "PAR","WindSpd")]

PeterWeather0[,"datetime"] <- paste(PeterWeather0[,"Date"], PeterWeather0[,"Time"])
PeterWeather0[,"datetime"] <- gsub("^(?=[0-9]/)", "0", PeterWeather0[,"datetime"], perl=TRUE)
PeterWeather0[,"datetime"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", PeterWeather0[,"datetime"], perl=TRUE)
PeterWeather0[,"datetime"] <- as.POSIXct(PeterWeather0[,"datetime"], format="%m/%d/%y %I:%M:%S %p", tz="GMT")
irr_wnd_2010 <- PeterWeather0[,c("datetime", "PAR", "WindSpd")]
names(irr_wnd_2010) <- c("datetime", "irr", "wnd")

save(irr_wnd_2010, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")

weather.full.2010 <- merge(irr_wnd_2010, UNDERC_Weather[,c("datetime","airTemp","RH","baro")], all.x=TRUE)
save(weather.full.2010, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/weather.full.2010.RData")




# ======================================
# = Peter and UNDERC Weather from 2012 =
# ======================================
# Peter Weather
Peter_PAR_Wind0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/Peter_PAR_Wind_2012.csv")
Peter_PAR_Wind0[,"datetime"] <- as.POSIXct(paste(Peter_PAR_Wind0[,"Date"], Peter_PAR_Wind0[,"Time"]), tz="GMT")
# Peter_PAR_Wind0[,"doy"] <- LakeMetabolizer:::date2doy(Peter_PAR_Wind0[,"datetime"])
# Peter_PAR_Wind[,"year"] <- 2012
Peter_PAR_Wind <- Peter_PAR_Wind0[,c("datetime", "PAR", "WindSpeed")]
# names(Peter_PAR_Wind) <- c("datetime", "irr", "wnd")

# Combine Peter and UNDERC weather
CombinePeterUNDERC0 <- merge(UNDERC_Weather[UNDERC_Weather[,"year"]==2012L,], Peter_PAR_Wind, all=TRUE)
flatIrr <- function(x){
	if(sum(!is.na(x[,"PAR"]))==nrow(x)){
		x[,"PAR"] <- x[,"PAR"] - min(x[,"PAR"])	
	}
	if(sum(!is.na(x[,"irr"]))==nrow(x)){
		x[,"irr"] <- x[,"irr"] - min(x[,"irr"])	
	}	
	
	x
	
}
CombinePeterUNDERC <- ddply(CombinePeterUNDERC0, c("year", "doy"), flatIrr)

# Convert UNDERC PAR to Peter PAR
PAR_Conv <- summary(lm(I(irr)~I(PAR) - 1, data=CombinePeterUNDERC))$coef[1]
CombinePeterUNDERC[,"irr"] <- (CombinePeterUNDERC[,"irr"] - 0) /PAR_Conv[1]


# Fill in missing Peter w/ UNDERC
UNDERC_2add.wnd <- is.na(CombinePeterUNDERC[,"WindSpeed"]) & !is.na(CombinePeterUNDERC[,"wnd"])
UNDERC_2add.irr <- is.na(CombinePeterUNDERC[,"PAR"]) & !is.na(CombinePeterUNDERC[,"irr"])
CombinePeterUNDERC[!UNDERC_2add.wnd, "WindSpeed"] <- wind.scale.base(CombinePeterUNDERC[!UNDERC_2add.wnd, "WindSpeed"], 2) # scale Peter wind
CombinePeterUNDERC[UNDERC_2add.wnd, "WindSpeed"] <- CombinePeterUNDERC[UNDERC_2add.wnd, "wnd"]
CombinePeterUNDERC[UNDERC_2add.irr, "PAR"] <- CombinePeterUNDERC[UNDERC_2add.irr, "irr"]

# Create the irr_wnd data
irr_wnd_2012 <- CombinePeterUNDERC[,c("datetime", "PAR", "WindSpeed")]
names(irr_wnd_2012) <- c("datetime", "irr", "wnd")
save(irr_wnd_2012, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2012.RData")

# Create the full weather data
weather.full.2012 <- CombinePeterUNDERC[,c("datetime", "airTemp", "RH", "WindSpeed", "PAR", "baro")]
names(weather.full.2012) <- c("datetime", "airTemp", "RH", "wnd", "irr", "baro")





