

library(plyr)
# detach(package:LakeMetabolizer, unload=TRUE)
# library("roxygen2")
# library("devtools")
# roxygenise("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", overwrite=FALSE)
# document("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", roclets="rd")
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lib/LakeMetabolizer", type="source", repos=NULL)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
# update.packages("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lib/LakeMetabolizer", type="source", repos=NULL)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")
library("rLakeAnalyzer")
library("zoo")

# Load weather file
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")

ciClasses <- c("character","numeric", "character", rep("numeric",3), rep("integer",5), rep("numeric",2))
# calInfo <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/WardSondes/WardCalData.csv", sep=",", header=TRUE, colClasses=ciClasses)
calInfo <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2010/WardCalData.csv", sep=",", header=TRUE, colClasses=ciClasses)

sites <- unique(calInfo[,"Site"])

# =============
# = Functions =
# =============
Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}
a.nc <- function(x) as.numeric(as.character(x)) # function to convert factor to numeric, but safe if starting as character, numeric, integer


sondes <- list()
for(i in 1:length(sites)){

	t.site <- sites[i]
	t.calInfo <- calInfo[calInfo[,"Site"]==t.site, ]
	t.files <- t.calInfo[,"FileName"]

	for(j in 1:nrow(t.calInfo)){

		# Read in a file
		t.f <- t.files[j]
		# t.fn <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/WardSondes/",t.f,sep="")
		t.fn <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2010/",t.f,sep="")
		t.dat <- read.table(t.fn, header=TRUE, sep=",")[-1,]
		# print(t.dat[1,])
		# Format dates in file
		if(t.f%in%c("Ward09Test.csv")){
			t.date <- as.POSIXct(paste(t.dat[,"Date"], t.dat[,"Time"]), format="%m/%d/%y %H:%M", tz="GMT")
		}else{
			t.date <- as.POSIXct(paste(t.dat[,"Date"], t.dat[,"Time"]), format="%m/%d/%y %H:%M:%S", tz="GMT")
		}
		
		t.year <- as.integer(format.Date(t.date, format="%Y"))
		t.doy <- LakeMetabolizer:::date2doy(t.date)
		
		# Create a data.frame of sonde data
		t.sonde0 <- data.frame(
			year=t.year, 
			date=t.date, 
			doy=t.doy, 
			doobs=a.nc(t.dat[,"DO.Conc"]), 
			dosat=a.nc(t.dat[,"DO."]), 
			wtr=a.nc(t.dat[,"Temp"])
		)
		
		# Compute values needed to correct for drift
		t.preCal <- mean(t.sonde0[1:t.calInfo[j,"Cal1"],"dosat"])
		t.postCal <- mean(t.sonde0[t.calInfo[j,"Cal4"]:nrow(t.sonde0),"dosat"])
		t.calSat <- c(t.preCal, t.postCal) # Average saturation before and after lake deployment
		
		t.lake.start <- t.calInfo[j, "Cal2"]
		t.lake.stop <- t.calInfo[j, "Cal3"]
		t.calTime <- c(t.sonde0[t.lake.start,"doy"], t.sonde0[t.lake.stop,"doy"]) # doy when initial cal ended, doy when final cal began
		
		t.sat.per.time <- diff(t.calSat)/diff(t.calTime) # average change in %sat (due to drift) per time (days)
		
		# Add NA's, drop days with too few observations, clip to in-lake observations, calculate drift
		t.sonde0.na1 <- LakeMetabolizer:::addNAs(t.sonde0[t.lake.start:t.lake.stop, ]) # KEEP! BELOW IS JUST TEST
		# t.sonde0.na1 <- LakeMetabolizer:::addNAs(t.sonde0[t.lake.start:t.lake.stop, ][-c(2240:2250),]) # TEST ONLY!
		t.elapsed.doy <- c(0, cumsum(diff(t.sonde0.na1[,"doy"])))
		t.cumul.drift.sat <- t.sat.per.time*t.elapsed.doy
		t.cumul.drift.conc <- Sat2Conc(SDO2sat=t.cumul.drift.sat, SDtemp=t.sonde0.na1[,"wtr"])
		
		# Create a data frame that has drift-corrected oxygen
		t.sonde0.na2 <- t.sonde0.na1
		t.sonde0.na2[,"dosat"] <- t.sonde0.na1[,"dosat"] - t.cumul.drift.sat
		t.sonde0.na2[,"doobs"] <- t.sonde0.na1[,"doobs"] - t.cumul.drift.conc
		
		t.sonde0.na2[,"datetime"] <- LakeMetabolizer:::round.time(t.sonde0.na2[,"datetime"], "5 minutes")
		
		t.sonde <- t.sonde0.na2[!duplicated(t.sonde0.na2[,"datetime"]),]
		t.sonde[,"zmix"] <- t.calInfo[j,"LayerBot"]
		
		t.sonde[,"top"] <- t.calInfo[j,"LayerTop"]
		t.sonde[,"bot"] <- t.calInfo[j,"LayerBot"]
		t.sonde[,"sensor_depth"] <- t.calInfo[j,"Depth"]
		t.sonde[,"baro"] <- t.calInfo[j,"CalPress"]*1.33322368
		t.sonde[,"z1perc"] <- t.calInfo[j, "X1Perc"]
		
		if(j==1){
			sondes[[t.site]] <- t.sonde
		}else{
			sondes[[t.site]] <- rbind(sondes[[t.site]], t.sonde)
		}
		
	} # end looping through files pertaining to a given site
		
		
} # end looping through sites



# ===========================
# = Read in Thermistor data =
# ===========================
# Can't actually use therm (or won't) b/c 1) don't have it for first few months of 2010, 2) need to remove out-of-lake obs
thermFiles <- list.files("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2010//WardThermComplete_2010/")
thermDepths <- gsub("WardTherm2010_([0-9]\\.[0-9])m.csv", "\\1", thermFiles)

for(i in 1:length(thermFiles)){
	# tFile <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/WardTherm/WardThermComplete_2010/", thermFiles[i], sep="")
	tFile <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/WardSondes_2010/WardThermComplete_2010/", thermFiles[i], sep="")
	tDepth <- as.character(as.numeric(thermDepths[i]))
	tTherm0 <- read.table(tFile, sep=",", header=FALSE)[,-1]
	names(tTherm0) <- c("datetime", paste("wtr",tDepth, sep="_"))
	tTherm0[,"datetime"] <- as.character(tTherm0[,"datetime"])
	tTherm0[,"datetime"] <- gsub("^(?=[0-9]/)", "0", tTherm0[,"datetime"], perl=TRUE)
	tTherm0[,"datetime"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", tTherm0[,"datetime"], perl=TRUE)
	tTherm0[,"datetime"] <- as.POSIXct(tTherm0[,"datetime"], format="%m/%d/%Y %H:%M", tz="GMT")
	tTherm0[,"datetime"] <- LakeMetabolizer:::round.time(tTherm0[,"datetime"], units="5 minutes")
	rmDups <- !duplicated(tTherm0[,"datetime"])
	tTherm <- tTherm0[rmDups,]
	
	# rollapply(tTherm[,2], width=288, by=288, scale)

	if(i!=1){
		tThermCum <- merge(tThermCum, tTherm, all=TRUE)
	}else{
		tThermCum <- tTherm
	}
}

w10.therm0 <- tThermCum[complete.cases(tThermCum),]


ward10.zmix <- ts.meta.depths(w10.therm0)
ward10.zmix.is0 <- ward10.zmix[,"top"] < 0.25 & !is.na(ward10.zmix[,"top"])
ward10.zmix[ward10.zmix.is0, "top"] <- 0.25

# =====================
# = Organize 2010 Epi =
# =====================
# merge weather and sonde
# LakeMetabolizer:::pred.merge(sondes[[2]], PeterWeather2010, all=TRUE)
ward10.epi0 <- merge(sondes[[2]], irr_wnd_2010, all.x=TRUE)

# fix names
names(ward10.epi0) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "baro", "z1perc", "irr", "wnd")
ward10.epi0 <- ward10.epi0[,c("datetime","do.obs", "do.sat", "wtr", "z.mix", "baro", "z1perc", "irr", "wnd")]


# add in thermistor zmix when available
ward10.epi.full <- merge(ward10.epi0, ward10.zmix[,c("datetime", "top")], all.x=TRUE)
w10.have.therm <- !is.na(ward10.epi.full[,"top"])
ward10.epi.full[w10.have.therm, "z.mix"] <- ward10.epi.full[w10.have.therm, "top"]


# format time
ward10.epi.full[,"datetime"] <- as.POSIXct(ward10.epi.full[,"datetime"])

# scale wind, calculate K, convert to K O2, scale K to sampling frequency
ward10.epi.full[,"wnd"] <- LakeMetabolizer:::scale.exp.wind.base(ward10.epi.full[,"wnd"], 2) # scale wind speed to 10 m
ward10.epi.full[,"k600.cole"] <- LakeMetabolizer:::k.cole.base(ward10.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward10.epi.full[,"k.gas"] <- LakeMetabolizer:::k600.2.kGAS.base(ward10.epi.full[,"k600.cole"], ward10.epi.full[,"wtr"], gas="O2")
ward10.epi.full[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(ward10.epi.full[,"wtr"], baro=ward10.epi.full[,"baro"])

ward10.epi <- ward10.epi.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr", "wnd")]


# ======================
# = Organize 2010 Meta =
# ======================
ward10.meta00 <- rbind(sondes[[3]], sondes[[8]], sondes[[10]])
LakeMetabolizer:::pred.merge(ward10.meta00, irr_wnd_2010, all=TRUE)
ward10.meta0 <- merge(ward10.meta00, irr_wnd_2010, all.x=TRUE)

# fix names
names(ward10.meta0) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "baro", "z1perc", "irr", "wnd")

# format time
ward10.meta0[,"datetime"] <- as.POSIXct(ward10.meta0[,"datetime"], tz="GMT")

# scale wind, calculate K, convert to K O2, scale K to sampling frequency
ward10.meta0[,"k.gas"] <- 0
ward10.meta0[,"do.sat"] <- o2.at.sat.base(ward10.meta0[,"wtr"], baro=ward10.meta0[,"baro"])

# pull in values needed to do a daily smoother of temperature
ward10.meta0[,"watts"] <- watts.in(ward10.meta0[,"top"], ward10.meta0[,"bot"], ward10.meta0[,"irr"], ward10.meta0[,"z1perc"])
ward10.meta0[,"roundDoy"] <- trunc(ward10.meta0[,"doy"])
ward10.meta.wtrSmooth <- ddply(ward10.meta0, .variables="roundDoy", function(x)cbind(x, "smooth.wtr"=temp.kalman(x[,"wtr"], x[,"watts"], ampH=50)))

ward10.meta0[,"smooth.wtr"] <- ward10.meta.wtrSmooth[,"smooth.wtr"]
ward10.meta.full <- ward10.meta0
ward10.meta <- ward10.meta0[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix", "irr", "wtr")]


# ==============================
# = Save 2010 Ward Sensor Data =
# ==============================
save(ward10.epi.full, ward10.epi, ward10.meta.full, ward10.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")

write.table(ward10.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.epi.full.txt", row.names=FALSE)
write.table(ward10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.epi.txt", row.names=FALSE)
write.table(ward10.meta.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.meta.full.txt", row.names=FALSE)
write.table(ward10.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.meta.txt", row.names=FALSE)





