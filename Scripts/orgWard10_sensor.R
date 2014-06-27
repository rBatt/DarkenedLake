

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

source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/lightFunctions.R")

# Load weather file
# load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/weather.full.2010.RData")


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
		# t.sonde[,"baro"] <- t.calInfo[j,"CalPress"]*1.33322368
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

ward10.therm <- merge(w10.therm0, ward10.zmix[,c("datetime","top")], all=TRUE)
names(ward10.therm)[names(ward10.therm)=="top"] <- "z.mix"



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
w10.light.prof <- light.prof[light.prof[,"lake"]=="Ward"&light.prof[,"year"]==2010L,]
w10.light.prof <- w10.light.prof[,!names(w10.light.prof)%in%c("lake","year")]





# =====================
# = Organize 2010 Epi =
# =====================
# merge weather and sonde
# LakeMetabolizer:::pred.merge(sondes[[2]], PeterWeather2010, all=TRUE)
# ward10.epi0 <- merge(sondes[[2]], irr_wnd_2010, all.x=TRUE)
ward10.epi0 <- merge(sondes[[2]], weather.full.2010, all.x=TRUE)

# fix names
names(ward10.epi0) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "z1perc", "irr", "wnd", "airTemp", "RH", "baro")
ward10.epi0 <- ward10.epi0[,c("datetime","do.obs", "do.sat", "wtr", "z.mix", "z1perc", "irr", "wnd", "airTemp", "RH", "baro")]


# add in thermistor zmix when available
ward10.epi.full <- merge(ward10.epi0, ward10.zmix[,c("datetime", "top")], all.x=TRUE)
w10.have.therm <- !is.na(ward10.epi.full[,"top"])
ward10.epi.full[w10.have.therm, "z.mix"] <- ward10.epi.full[w10.have.therm, "top"]

# format time
ward10.epi.full[,"datetime"] <- as.POSIXct(ward10.epi.full[,"datetime"], tz="GMT")


# ==========================================
# = Calculate Light data for Ward Epi 2010 =
# ==========================================
# Calculate the average zmix for each day
w10.light000 <- LakeMetabolizer:::addNAs(ward10.epi.full)[,c("datetime","doy","z.mix", "irr")]
# w12.light000 <- cbind("lake"="Ward", LakeMetabolizer:::addNAs(ward12.epi.full)[,c("datetime","doy","z.mix", "irr")])
# l10.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul10.epi.full)[,c("datetime","doy","z.mix", "irr")])
# l12.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul12.epi.full)[,c("datetime","doy","z.mix", "irr")])

w10.light000[,"doy"] <- trunc(w10.light000[,"doy"])
w10.light000[,"year"] <- format.Date(w10.light000[,"datetime"], format="%Y")
w10.light00 <- w10.light000[complete.cases(w10.light000),]
light0 <- merge(w10.light00, w10.light.prof, all=TRUE)
w10.light0 <- doLight(light0)
w10.light <- w10.light0[,c("datetime", "irr.sonde", "irr.bot", "dz.sonde", "dz.bot", "kd")]

ward10.epi.full <- merge(ward10.epi.full, w10.light, all.x=TRUE)

# ==========================================
# = Calculate other params for Ward Epi 10 =
# ==========================================
# scale wind, calculate K, convert to K O2, scale K to sampling frequency
ward10.epi.full[,"wnd"] <- LakeMetabolizer:::wind.scale.base(ward10.epi.full[,"wnd"], 2) # scale wind speed to 10 m
ward10.epi.full[,"k600.cole"] <- LakeMetabolizer:::k.cole.base(ward10.epi.full[,"wnd"]) # calculate k600 using Cole & Caraco method
ward10.epi.full[,"kgas.cole"] <- LakeMetabolizer:::k600.2.kGAS.base(ward10.epi.full[,"k600.cole"], ward10.epi.full[,"wtr"], gas="O2")
w10.baro <- ward10.epi.full[,"baro"]
w10.baro[is.na(ward10.epi.full[,"baro"])] <- 960
ward10.epi.full[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(ward10.epi.full[,"wtr"], baro=ward10.epi.full[,"baro"])

# calculate k.read
ward10.sw <- par.to.sw.base(ward10.epi.full[,"irr"])
ward10.lw <- calc.lw.net.base(
	dateTime=ward10.epi.full[,"datetime"], 
	sw=ward10.sw, 
	Ts=ward10.epi.full[,"wtr"], 
	lat=46.28, 
	atm.press=ward10.epi.full[,"baro"], 
	airT=ward10.epi.full[,"airTemp"],
	RH=ward10.epi.full[,"RH"]
	)

# miss.read <- is.na(Ts=ward10.epi.full[,"wtr"]) | is.na(ward10.epi.full[,"kd"])
miss.read <- is.na(ward10.epi.full[,"wtr"]) | is.na(ward10.epi.full[,"kd"]) | is.na(ward10.epi.full[,"baro"]) | is.na(ward10.epi.full[,"airTemp"]) | is.na(ward10.epi.full[,"RH"]) | is.na(ward10.epi.full[,"z.mix"])

ward10.epi.full[!miss.read,"k600.read"] <- k.read.base(
	wnd.z=10,
	Kd=ward10.epi.full[!miss.read,"kd"],
	lat=46.28,
	lake.area=19000,
	atm.pres=ward10.epi.full[!miss.read,"baro"],
	dateTime=ward10.epi.full[!miss.read,"datetime"],
	Ts=ward10.epi.full[!miss.read,"wtr"],
	z.aml=ward10.epi.full[!miss.read,"z.mix"],
	airT=ward10.epi.full[!miss.read,"airTemp"],
	wnd=ward10.epi.full[!miss.read,"wnd"],
	RH=ward10.epi.full[!miss.read,"RH"],
	sw=ward10.sw[!miss.read],
	lwnet=ward10.lw[!miss.read]
	)

miss.read2 <- is.na(ward10.epi.full[,"k600.read"])

# dev.new(width=3.5, height=5)
# par(mfrow=c(2,1), mar=c(2.25, 2.25, 0.5, 0.5), mgp=c(1.25, 0.35, 0), tcl=-0.25, ps=10, family="Times")
# plot(ward10.epi.full[!miss.read,"datetime"], ward10.epi.full[!miss.read,"k600.read"], type="l", xlim=range(ward10.epi.full[,"datetime"]), xlab="", ylab=bquote(k600~~(m~d^-1)))
# lines(ward10.epi.full[,"datetime"], ward10.epi.full[,"k600.cole"], type="l", col="red")
# legend("topright", legend=c("k.cole", "k.read"), lty=1, lwd=1.5, col=c("red","black"))
# plot(log(ward10.epi.full[!miss.read,"k600.cole"]), ward10.epi.full[!miss.read,"k600.read"], xlab=bquote(log[e]~k600.cole~~(m~d^-1)), ylab=bquote(k600.read~~(m~d^-1)))

ward10.log.cole <- log(ward10.epi.full[,"k600.cole"])
ward10.lm.read.cole <- lm(ward10.epi.full[,"k600.read"]~ward10.log.cole)
ward10.read.from.cole <- predict(ward10.lm.read.cole, newdata=data.frame(ward10.log.cole=ward10.log.cole[miss.read2]))
ward10.epi.full[miss.read2,"k600.read"] <- as.numeric(ward10.read.from.cole)
ward10.epi.full[,"kgas.read"] <- LakeMetabolizer:::k600.2.kGAS.base(ward10.epi.full[,"k600.read"], ward10.epi.full[,"wtr"], gas="O2")

w10.noK <- ward10.epi.full[,"z.mix"] < 0.7 & !is.na(ward10.epi.full[,"z.mix"])
ward10.epi.full[w10.noK,"kgas.cole"] <- 0L
ward10.epi.full[w10.noK,"kgas.read"] <- 0L

ward10.epi <- ward10.epi.full[,c("datetime", "do.obs", "do.sat", "kgas.cole", "kgas.read", "z.mix",  "irr", "wtr", "wnd")]




# ======================
# = Organize 2010 Meta =
# ======================
ward10.meta00 <- rbind(sondes[[3]], sondes[[8]], sondes[[10]])
ward10.meta0 <- merge(ward10.meta00, weather.full.2010, all.x=TRUE)

# fix names
names(ward10.meta0) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "z1perc", "irr", "wnd", "airTemp", "RH", "baro")

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
save(ward10.therm, ward10.epi.full, ward10.epi, ward10.meta.full, ward10.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")

write.table(ward10.therm, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.therm.txt", row.names=FALSE)
write.table(ward10.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.epi.full.txt", row.names=FALSE)
write.table(ward10.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.epi.txt", row.names=FALSE)
write.table(ward10.meta.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.meta.full.txt", row.names=FALSE)
write.table(ward10.meta, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward10.meta.txt", row.names=FALSE)





