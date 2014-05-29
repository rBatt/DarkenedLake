
library(plyr)
install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lib/LakeMetabolizer", type="source", repos=NULL)
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
		print(t.dat[1,])
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
		# t.sonde0.na1 <- LakeMetabolizer:::addNAs(t.sonde0[t.lake.start:t.lake.stop, ]) # KEEP! BELOW IS JUST TEST
		t.sonde0.na1 <- LakeMetabolizer:::addNAs(t.sonde0[t.lake.start:t.lake.stop, ][-c(2240:2250),]) # TEST ONLY!
		t.elapsed.doy <- c(0, cumsum(diff(t.sonde0.na1[,"doy"])))
		t.cumul.drift.sat <- t.sat.per.time*t.elapsed.doy
		t.cumul.drift.conc <- Sat2Conc(SDO2sat=t.cumul.drift.sat, SDtemp=t.sonde0.na1[,"wtr"])
		
		# Create a data frame that has drift-corrected oxygen
		t.sonde0.na2 <- t.sonde0.na1
		t.sonde0.na2[,"dosat"] <- t.sonde0.na1[,"dosat"] - t.cumul.drift.sat
		t.sonde0.na2[,"doobs"] <- t.sonde0.na1[,"doobs"] - t.cumul.drift.conc
		
		t.sonde0.na2[,"date"] <- round.time(t.sonde0.na2[,"date"], "5 minutes")
		
		t.sonde <- t.sonde0.na2[!duplicated(t.sonde0.na2[,"date"]),]
		if(t.site=="Epi" | t.site=="Litt"){
			t.sonde[,"zmix"] <- t.calInfo[j,"LayerBot"]
		}
		
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
PeterWeather00 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Metabolism/Raw Weather Data/PeterWeather.csv", header=TRUE)
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
thermFiles <- list.files("/Users/Battrd/Documents/School&Work/WiscResearch/WardTherm/WardThermComplete_2010/")
thermDepths <- gsub("WardTherm2010_([0-9]\\.[0-9])m.csv", "\\1", thermFiles)

for(i in 1:length(thermFiles)){
	tFile <- paste("/Users/Battrd/Documents/School&Work/WiscResearch/WardTherm/WardThermComplete_2010/", thermFiles[i], sep="")
	tDepth <- as.character(as.numeric(thermDepths[i]))
	tTherm0 <- read.table(tFile, sep=",", header=FALSE)[,-1]
	names(tTherm0) <- c("DateTime", paste("wtr",tDepth, sep="_"))
	tTherm0[,"DateTime"] <- as.character(tTherm0[,"DateTime"])
	tTherm0[,"DateTime"] <- gsub("^(?=[0-9]/)", "0", tTherm0[,"DateTime"], perl=TRUE)
	tTherm0[,"DateTime"] <- gsub("(?<=[0-9]{2}/)([0-9])(?=/)", "0\\1", tTherm0[,"DateTime"], perl=TRUE)
	tTherm0[,"DateTime"] <- as.character(as.POSIXct(tTherm0[,"DateTime"], format="%m/%d/%Y %H:%M"))
	tTherm0[,"DateTime"] <- round.time(tTherm0[,"DateTime"], units="5 minutes")
	rmDups <- !duplicated(tTherm0[,"DateTime"])
	tTherm <- tTherm0[rmDups,]
	
	# rollapply(tTherm[,2], width=288, by=288, scale)

	if(i!=1){
		tThermCum <- merge(tThermCum, tTherm, all=TRUE)
	}else{
		tThermCum <- tTherm
	}
}

tTherm[15000:15010,]


LakeMetabolizer:::pred.merge(sondes[[2]], PeterWeather2010, all=TRUE)
test <- merge(sondes[[2]], PeterWeather2010, all.x=TRUE)
names(test) <- c("datetime","year", "doy", "do.obs", "do.sat", "wtr", "z.mix", "irr", "wnd")
test[,"k.gas"] <- k.cole(test[,c("datetime","wnd")])[,2] # NEED TO USE k600.2.kGAS
test[,"datetime"] <- as.POSIXct(test[,"datetime"])
test[,"doy"] <- LakeMetabolizer:::date2doy(test[,"datetime"])
test[,"do.sat"] <- LakeMetabolizer:::o2.at.sat(test[,"wtr"]) # NEED TO FINISH THIS USE OF o2.at.sat

LakeMetabolizer:::metab(test, "mle")






