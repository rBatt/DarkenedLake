
library(plyr)
# detach(package:LakeMetabolizer, unload=TRUE)
install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lib/LakeMetabolizer", type="source", repos=NULL)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
# update.packages("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lib/LakeMetabolizer", type="source", repos=NULL)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")
library("rLakeAnalyzer")
library("zoo")

ciClasses <- c("character","numeric", "character", rep("numeric",3), rep("integer",5), rep("numeric",2))
calInfo <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/WardSondes/WardCalData.csv", sep=",", header=TRUE, colClasses=ciClasses)

sites <- unique(calInfo[,"Site"])

# =============
# = Functions =
# =============
Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}
a.nc <- function(x) as.numeric(as.character(x)) # function to convert factor to numeric, but safe if starting as character, numeric, integer

# =======================
# = round.time Function =
# =======================
round.time <- function(x, units, input.format=NULL, output.format="%Y-%m-%d %H:%M:%S"){
	# x = head(t.sonde0.na2[,"date"], 20) + 120
	# units = "df.345 min"
	# Check for invalid input classes
	stopifnot(is.character(units)&is.character(output.format)&(is.null(input.format)|is.character(input.format)))
	
	# Determine time unit
	unit.choices <- c("sec", "min", "hour", "day")
	choices.or <- paste(unit.choices, collapse="|")
	unit.pattern <- paste(".*(", choices.or, ").*", sep="")
	unit <- gsub(unit.pattern, "\\1", units)
	if(is.na(unit)){stop("not a valid unit, use sec, min, hour, or day")}
	which.choice <- which(unit==unit.choices)
	
	# Determine time interval
	u.time.pattern <- "(?:[0-9]+\\.[0-9]+)|(?:[0-9]+\\.)|(?:\\.[0-9]+)|(?:[0-9]+)"
	u.time.char <- regmatches(units, regexpr(u.time.pattern, units, perl=TRUE))
	u.time <- as.numeric(u.time.char)
	u.time <- ifelse(is.na(u.time), 1, u.time)
	
	unit.cutoff <- switch(unit, sec=60, min=60, hour=24, day=1)
	
	# =========================================================================
	# = Check for invalid input (before slow [attempted] conversion to POSIX) =
	# =========================================================================
	if(sign(u.time)==-1L){
		stop("time interval must be positive")
	}
	# Deal with case where units are 1 second (or less)
	if(unit=="sec" & u.time<=1L){
		return(format.Date(x, format=output.format))
	} else
	
	# Fractional time intervals – convert to smaller unit
	if((trunc(u.time)-u.time)!=0){
		if(sign(u.time)==1L){
			while((trunc(u.time)-u.time)!=0){
				if(unit=="sec"){stop("time interval must be an integer when converted to units of seconds")}
				unit <- unit.choices[which.choice-1]
				which.choice <- which(unit==unit.choices)
				unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
				u.time <- unit.cutoff*u.time
			}
		}else{
			stop("time interval must be positive")
		}
	} else 
	
	# Deal with case where units are days
	if(unit=="day"){
		if(u.time==1){
			return(format.Date(trunc.POSIXt(x + 43200, units = units), format=output.format))
		}else{
			stop("units must be <= 1 day")
		}
	} else 
	
	# Deal w/ cases where time interval is 1 unit
	if(u.time==1){
			unit <- unit.choices[which.choice-1]
			which.choice <- which(unit==unit.choices)
			unit.cutoff <- switch(unit, sec=60, min=60, hour=24)
			u.time <- unit.cutoff
	} 
	
	# Deal with cases where time interval is > 1 of a larger unit
	# Note that this follows up on case where u.time is > 1 and is not an integer
	if(u.time>unit.cutoff){
		u.time <- u.time%%unit.cutoff
		mod.mess <- paste("Rounding to units =", u.time, unit) # may or may not want to make this a warning ...
		warning(mod.mess)
	}
	
	# =============================================
	# = Convert to POSIX, or if can't, give error =
	# =============================================
	if(all(class(x)!="POSIX.ct")){
		if(is.null(input.format)){
			x <- as.POSIXct(x)
		}else{
			x <- as.POSIXct(x, format=input.format)
		}
	}
	
	# ===========================================================
	# = Matching units (e.g., min) and unit multiples (e.g., 5) =
	# ===========================================================
	
	which.choice <- which(unit==unit.choices)
	form.unit <- c("%S", "%M", "%H", "%d")[which.choice]
	mult <- as.integer(format.Date(x, format=form.unit))/u.time
	after <- round(mult, 0)*u.time
	# direction <- sign(after-before)
	
	# trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
	trunc.unit <- unit.choices[min(which.choice+1, length(unit.choices))]
	rounded <- trunc.POSIXt(x, trunc.unit) + switch(unit, sec = 1, min = 60, hour = 3600, day = 86400)*after
	format.Date(rounded, format=output.format)
}

sondes <- list()
for(i in 1:length(sites)){

	t.site <- sites[i]
	t.calInfo <- calInfo[calInfo[,"Site"]==t.site, ]
	t.files <- t.calInfo[,"FileName"]

	for(j in 1:nrow(t.calInfo)){

		# Read in a file
		t.f <- t.files[j]
		t.fn <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/WardSondes/",t.f,sep="")
		t.dat <- read.table(t.fn, header=TRUE, sep=",")[-1,]
		# print(t.dat[1,])
		# Format dates in file
		if(t.f%in%c("Ward09Test.csv")){
			t.date <- as.POSIXct(paste(t.dat[,"Date"], t.dat[,"Time"]), format="%m/%d/%y %H:%M")
		}else{
			t.date <- as.POSIXct(paste(t.dat[,"Date"], t.dat[,"Time"]), format="%m/%d/%y %H:%M:%S")
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
		
		t.sonde0.na2[,"date"] <- round.time(t.sonde0.na2[,"date"], "5 minutes")
		
		t.sonde <- t.sonde0.na2[!duplicated(t.sonde0.na2[,"date"]),]
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
		


# ===================================
# = Read in Peter Weather from 2010 =
# ===================================
PeterWeather00 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/Raw Weather Data/PeterWeather.csv", header=TRUE)[-(96344:96345),]
PeterWeather0 <- PeterWeather00[PeterWeather00[,"Year"]==2010L,]
PeterWeather0 <- PeterWeather0[,c("Date", "Time", "Year", "PAR","WindSpd")]

PeterWeather0[,"date"] <- paste(PeterWeather0[,"Date"], PeterWeather0[,"Time"])
PeterWeather0[,"date"] <- gsub("^(?=[0-9]/)", "0", PeterWeather0[,"date"], perl=TRUE)
PeterWeather0[,"date"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", PeterWeather0[,"date"], perl=TRUE)
PeterWeather0[,"date"] <- as.character(as.POSIXct(PeterWeather0[,"date"], format="%m/%d/%y %I:%M:%S %p"))
PeterWeather2010 <- PeterWeather0[,c("date", "PAR", "WindSpd")]



# ===========================
# = Read in Thermistor data =
# ===========================
# Can't actually use therm (or won't) b/c 1) don't have it for first few months of 2010, 2) need to remove out-of-lake obs
# thermFiles <- list.files("/Users/Battrd/Documents/School&Work/WiscResearch/WardTherm/WardThermComplete_2010/")
# thermDepths <- gsub("WardTherm2010_([0-9]\\.[0-9])m.csv", "\\1", thermFiles)
# 
# for(i in 1:length(thermFiles)){
# 	tFile <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/WardTherm/WardThermComplete_2010/", thermFiles[i], sep="")
# 	tDepth <- as.character(as.numeric(thermDepths[i]))
# 	tTherm0 <- read.table(tFile, sep=",", header=FALSE)[,-1]
# 	names(tTherm0) <- c("DateTime", paste("wtr",tDepth, sep="_"))
# 	tTherm0[,"DateTime"] <- as.character(tTherm0[,"DateTime"])
# 	tTherm0[,"DateTime"] <- gsub("^(?=[0-9]/)", "0", tTherm0[,"DateTime"], perl=TRUE)
# 	tTherm0[,"DateTime"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", tTherm0[,"DateTime"], perl=TRUE)
# 	tTherm0[,"DateTime"] <- as.character(as.POSIXct(tTherm0[,"DateTime"], format="%m/%d/%Y %H:%M"))
# 	tTherm0[,"DateTime"] <- round.time(tTherm0[,"DateTime"], units="5 minutes")
# 	rmDups <- !duplicated(tTherm0[,"DateTime"])
# 	tTherm <- tTherm0[rmDups,]
# 	
# 	# rollapply(tTherm[,2], width=288, by=288, scale)
# 
# 	if(i!=1){
# 		tThermCum <- merge(tThermCum, tTherm, all=TRUE)
# 	}else{
# 		tThermCum <- tTherm
# 	}
# }

# tTherm[15000:15010,]

Mode <- function(x){
		ux <- unique(x)
		ux[which.max(tabulate(match(x, ux)))]
}


# =====================
# = Organize 2010 Epi =
# =====================
# merge weather and sonde
LakeMetabolizer:::pred.merge(sondes[[2]], PeterWeather2010, all=TRUE)
ward10.epi <- merge(sondes[[2]], PeterWeather2010, all.x=TRUE)

# fix names
names(ward10.epi) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "baro", "z1perc", "irr", "wnd")

# format time
ward10.epi[,"datetime"] <- as.POSIXct(ward10.epi[,"datetime"])
ward10.epi[,"doy"] <- LakeMetabolizer:::date2doy(ward10.epi[,"datetime"])

# calculate sampling frequency
freq.ward10.epi <- round(Mode(1/diff(ward10.epi[,"doy"])))

# scale wind, calculate K, convert to K O2, scale K to sampling frequency
ward10.epi.scale10 <- scale.exp.wind.base(ward10.epi[,"wnd"], 2) # scale wind speed to 10 m
ward10.epi.k.cole <- k.cole.base(ward10.epi.scale10) # calculate k600 using Cole & Caraco method
ward10.epi[,"k.gas"] <- k600.2.kGAS.base(ward10.epi.k.cole, ward10.epi[,"wtr"], gas="O2")/freq.ward10.epi
ward10.epi[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(ward10.epi[,"wtr"], baro=ward10.epi[,"baro"])

ward10.epi <- ward10.epi[,c("datetime","year", "doy", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr", "wnd")]

ward10.epi2 <- ward10.epi[, c("year","doy","datetime","do.obs", "do.sat", "k.gas", "z.mix", "irr", "wtr")]



# ======================
# = Organize 2010 Meta =
# ======================
ward10.meta00 <- rbind(sondes[[3]], sondes[[8]], sondes[[10]])
LakeMetabolizer:::pred.merge(ward10.meta00, PeterWeather2010, all=TRUE)
ward10.meta0 <- merge(ward10.meta00, PeterWeather2010, all.x=TRUE)

# fix names
names(ward10.meta0) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "top", "bot", "sensor_depth", "baro", "z1perc", "irr", "wnd")

# format time
ward10.meta0[,"datetime"] <- as.POSIXct(ward10.meta0[,"datetime"])
ward10.meta0[,"doy"] <- LakeMetabolizer:::date2doy(ward10.meta0[,"datetime"])

# scale wind, calculate K, convert to K O2, scale K to sampling frequency
ward10.meta0[,"k.gas"] <- 0
ward10.meta0[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(ward10.meta0[,"wtr"], baro=ward10.meta0[,"baro"])

# pull in values needed to do a daily smoother of temperature
ward10.meta0[,"watts"] <- watts.in(ward10.meta0[,"top"], ward10.meta0[,"bot"], ward10.meta0[,"irr"], ward10.meta0[,"z1perc"])
ward10.meta0[,"roundDoy"] <- trunc(ward10.meta0[,"doy"])

ward10.meta.wtrSmooth <- ddply(ward10.meta0, .variables="roundDoy", function(x)cbind(x, "smooth.wtr"=temp.kalman(x[,"wtr"], x[,"watts"], ampH=50)))

ward10.meta0[,"wtr"] <- ward10.meta.wtrSmooth[,"smooth.wtr"]
ward10.meta <- ward10.meta0[,c("datetime", "year", "doy", "do.obs", "do.sat", "k.gas", "z.mix", "irr", "wtr")]




# ==================
# = Epi Metabolism =
# ==================
ward10.epi2 <- ward10.epi[, c("year","doy","datetime","do.obs", "do.sat", "k.gas", "z.mix", "irr", "wtr")]
# MLE
ward10.epi.mle <- LakeMetabolizer:::metab(ward10.epi, "mle")
ward10.epi.mle.res <- ward10.epi.mle[[2]][,c("doy","GPP","R", "NEP")]
ward10.epi.mle.res <- merge(data.frame("doy"=138:239), ward10.epi.mle.res, all=TRUE)

# Kalman
ward10.epi.kal <- LakeMetabolizer:::metab(ward10.epi, "kalman")
ward10.epi.kal.res <- ward10.epi.kal[[3]][,c("doy","GPP","R", "NEP")]
ward10.epi.kal.res <- merge(data.frame("doy"=138:239), ward10.epi.kal.res, all=TRUE)

# Bayesian
ward10.epi.bay <- LakeMetabolizer:::metab(ward10.epi2, "bayesian")
ward10.epi.bay.res <- ward10.epi.bay[[4]][,c("doy","GPP","R", "NEP")]
ward10.epi.bay.res <- merge(data.frame("doy"=138:239), ward10.epi.bay.res, all=TRUE)

# Reformat data for bookkeeping, do BK

ward10.epi3 <- ward10.epi2
bk.irr <- as.integer(LakeMetabolizer:::is.day(46.28, ward10.epi3[,"datetime"]))
ward10.epi3[,"irr"] <- bk.irr
ward10.epi.bk <- LakeMetabolizer:::metab(ward10.epi3, "bookkeep")
ward10.epi.bk.res <- ward10.epi.bk[,c("doy","GPP","R", "NEP")]
ward10.epi.bk.res <- merge(data.frame("doy"=138:239), ward10.epi.bk.res, all=TRUE)

# OLS
ward10.epi.ols <- LakeMetabolizer:::metab(ward10.epi2, "ols")
ward10.epi.ols.res <- ward10.epi.ols[,c("doy","GPP","R", "NEP")]
ward10.epi.ols.res <- merge(data.frame("doy"=138:239), ward10.epi.ols.res, all=TRUE)

# merge all results
w10e.r1 <- merge(cbind("method"="mle",ward10.epi.mle.res), cbind("method"="kalman",ward10.epi.kal.res), all=TRUE)
w10e.r2 <- merge(w10e.r1, cbind("method"="bayes", ward10.epi.bay.res), all=TRUE)
w10e.r3 <- merge(w10e.r2, cbind("method"="bookkeep", ward10.epi.bk.res), all=TRUE)
w10e.r4 <- merge(w10e.r3, cbind("method"="ols", ward10.epi.ols.res), all=TRUE)


# ================================================
# = Plot Ward 2010 Epi Metabolism across Methods =
# ================================================
# dev.new(width=3.5, height=7)
pdf("~/Desktop/ward2010_epilimnion.pdf", width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="ols","R"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="mle","R"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="mle","R"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","R"], w10e.r4[w10e.r4[,"method"]=="kalman","R"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","R"], w10e.r4[w10e.r4[,"method"]=="bayes","R"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nRespiration", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="ols","GPP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="mle","GPP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="mle","GPP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","GPP"], w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","GPP"], w10e.r4[w10e.r4[,"method"]=="bayes","GPP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nGross Primary Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="ols","NEP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="mle","NEP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="bookkeep","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="mle","NEP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="ols","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="mle","NEP"], w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10e.r4[w10e.r4[,"method"]=="mle","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10e.r4[w10e.r4[,"method"]=="kalman","NEP"], w10e.r4[w10e.r4[,"method"]=="bayes","NEP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Epi\nNet Ecosystem Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)

dev.off()




# ==================
# = Meta Metabolism =
# ==================
# MLE
ward10.meta.mle <- LakeMetabolizer:::metab(ward10.meta, "mle")
ward10.meta.mle.res <- ward10.meta.mle[[2]][,c("doy","GPP","R", "NEP")]
ward10.meta.mle.res <- merge(data.frame("doy"=138:239), ward10.meta.mle.res, all=TRUE)

# Kalman
ward10.meta.kal <- LakeMetabolizer:::metab(ward10.meta, "kalman")
ward10.meta.kal.res <- ward10.meta.kal[[3]][,c("doy","GPP","R", "NEP")]
ward10.meta.kal.res <- merge(data.frame("doy"=138:239), ward10.meta.kal.res, all=TRUE)

# Bayesian
ward10.meta2 <- ward10.meta[, c("year","doy","datetime","do.obs", "do.sat", "k.gas", "z.mix", "irr", "wtr")]
ward10.meta.bay <- LakeMetabolizer:::metab(ward10.meta2, "bayesian")
ward10.meta.bay.res <- ward10.meta.bay[[4]][,c("doy","GPP","R", "NEP")]
ward10.meta.bay.res <- merge(data.frame("doy"=138:239), ward10.meta.bay.res, all=TRUE)

# Reformat data for bookkemetang, do BK
ward10.meta3 <- ward10.meta2
bk.irr <- as.integer(LakeMetabolizer:::is.day(46.28, ward10.meta3[,"datetime"]))
ward10.meta3[,"irr"] <- bk.irr
ward10.meta.bk <- LakeMetabolizer:::metab(ward10.meta3, "bookkeep")
ward10.meta.bk.res <- ward10.meta.bk[,c("doy","GPP","R", "NEP")]
ward10.meta.bk.res <- merge(data.frame("doy"=138:239), ward10.meta.bk.res, all=TRUE)

# OLS
ward10.meta.ols <- LakeMetabolizer:::metab(ward10.meta2, "ols")
ward10.meta.ols.res <- ward10.meta.ols[,c("doy","GPP","R", "NEP")]
ward10.meta.ols.res <- merge(data.frame("doy"=138:239), ward10.meta.ols.res, all=TRUE)

# merge all results
w10m.r1 <- merge(cbind("method"="mle",ward10.meta.mle.res), cbind("method"="kalman",ward10.meta.kal.res), all=TRUE)
w10m.r2 <- merge(w10m.r1, cbind("method"="bayes", ward10.meta.bay.res), all=TRUE)
w10m.r3 <- merge(w10m.r2, cbind("method"="bookkeep", ward10.meta.bk.res), all=TRUE)
w10m.r4 <- merge(w10m.r3, cbind("method"="ols", ward10.meta.ols.res), all=TRUE)



# ================================================
# = Plot Ward 2010 Epi Metabolism across Methods =
# ================================================
pdf("~/Desktop/ward2010_metalimnion.pdf", width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="ols","R"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="mle","R"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="mle","R"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","R"], w10m.r4[w10m.r4[,"method"]=="kalman","R"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","R"], w10m.r4[w10m.r4[,"method"]=="bayes","R"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\n Respiration", outer=TRUE, line=0, side=3, font=2, cex=1.2)



par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="ols","GPP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="mle","GPP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="mle","GPP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","GPP"], w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","GPP"], w10m.r4[w10m.r4[,"method"]=="bayes","GPP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\nGross Primary Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)


# dev.new(width=3.5, height=7)
par(mfcol=c(5,2), mar=c(2, 2, 0.5, 0.5), oma=c(0.1, 0.1, 2, 0.1), ps=9, mgp=c(1, 0.25, 0), tcl=-0.35, family="Times", cex=1)

plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="ols","NEP"], xlab="BK", ylab="OLS"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="mle","NEP"], xlab="BK", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="BK", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="bookkeep","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="BK", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="mle","NEP"], xlab="OLS", ylab="MLE"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="OLS", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="ols","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="OLS", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="mle","NEP"], w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], xlab="MLE", ylab="Kalman"); abline(a=0, b=1)
plot(w10m.r4[w10m.r4[,"method"]=="mle","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="MLE", ylab="Bayes"); abline(a=0, b=1)

plot(w10m.r4[w10m.r4[,"method"]=="kalman","NEP"], w10m.r4[w10m.r4[,"method"]=="bayes","NEP"], xlab="Kalman", ylab="Bayes"); abline(a=0, b=1)

mtext("Ward 2010 Meta\nNet Ecosystem Production", outer=TRUE, line=0, side=3, font=2, cex=1.2)

dev.off()


