

# =================
# = Load packages =
# =================
library(rLakeAnalyzer)
library(LakeMetabolizer)

# ===================
# = Write Functions =
# ===================
Sat2Conc <- function(SDO2sat, SDtemp){ # convert DO % sat to mg/L
	SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5+0.000002672016*SDtemp^4+(-0.0001884085)*SDtemp^3+0.009778012*SDtemp^2+(-0.4147241)*SDtemp+14.621)
	return(SDO2mgL)
}

# =====================
# = Load Weather Data =
# =====================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PeterWeather_2010&2012/irr_wnd_2012.RData")



# ================================
# = Read in Paul 2012 Thermistor =
# ================================
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/PaulSondes_2012/"
paul12.therm00 <- read.table(paste(datDir, "Paul_tChain_2012.csv", sep=""), sep=",", header=TRUE)
paul12.therm00[,"datetime"] <- as.POSIXct(paul12.therm00[,"datetime"], tz="GMT")

# Deal with breaking thermistors (rLakeAnalyzer does not accept NAs)
zPatt <- function(x){
	x <- as.numeric(x[-1])
	paste(cumsum(!is.na(x) & x>0 & x<40), collapse="")
}
depthPresent00 <- factor(apply(paul12.therm00, 1, zPatt))

paul12.therm0 <- paul12.therm00[depthPresent00!=levels(depthPresent00)[1],]
depthPresent0 <- factor(depthPresent00[depthPresent00!=levels(depthPresent00)[1]])

paul12.therm <- paul12.therm0
for(i in 1:length(levels(depthPresent0))){ # this loop is only 2 iterations; slowness due to rLakeAnalyzer
	t.lvl <- levels(depthPresent0)[i]
	t.ind <- depthPresent0==t.lvl
	
	t.good.col <- apply(paul12.therm[t.ind,], 2, function(x){all(!is.na(x))})	
	
	paul12.therm[t.ind, "z.mix"] <- ts.meta.depths(paul12.therm[t.ind, t.good.col])[,"top"]
}


# Read in manual zmix
man.zimx00 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv") # contains all weekly
manSub.ind <- man.zimx00[,"Lake"]=="Paul" & man.zimx00[,"Layer"]=="PML" & man.zimx00[,"Year"]==2012L # subset manual zmix data to needed

man.zimx0 <- man.zimx00[manSub.ind, c("Date","Zmix")] # add in manual zmix where tChain zmix is missing
man.zimx0[,"datetime"] <- as.POSIXct(man.zimx0[,"Date"], tz="GMT")

mz.interp.date <- seq(range(man.zimx0[,"datetime"])[1], to=range(man.zimx0[,"datetime"])[2]+(287+288)*5*60, by=60*5)
mz.interp.zmix <- approx(x=man.zimx0[,"datetime"], y=man.zimx0[,"Zmix"], xout=mz.interp.date, method="constant", rule=2, f=0)$y
man.z.mix <- data.frame(datetime=mz.interp.date, manZ=mz.interp.zmix)

paul12.zmix0 <- merge(paul12.therm[,c("datetime","z.mix")], man.z.mix, all=TRUE)
miss.tchain.z <- is.na(paul12.zmix0[,"z.mix"])
paul12.zmix0[miss.tchain.z,"z.mix"] <- paul12.zmix0[miss.tchain.z,"manZ"]

paul12.zmix <- paul12.zmix0[,c("datetime", "z.mix")]







# ================================
# = Read in Paul 2012 Sonde Data =
# ================================
FileNames <- c("50312L1D.csv", "52312L1D.csv", "61512L1D.csv", "70612L1D.csv", "72712L1D.csv", "81712L1D.csv")
for(i in 1:length(FileNames)){
	tName <- FileNames[i]
	ReadClassTxt <- c(NA, NA, rep("numeric", 10))
	tDat <- read.table(paste(datDir, tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt) # read in data, skipping gibberish
	names(tDat) <- c("Date", "Time", "wtr", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "do.obs", "Chl_conc", "Chl_RFU", "Battery")
	tDat[,"datetime"] <- as.POSIXct(paste(tDat[,"Date"], tDat[,"Time"]), tz="GMT")
	tDat[,"datetime"] <- LakeMetabolizer:::round.time(tDat[,"datetime"], "5 minutes") # round to the nearest 5 minute interval
	
	tDat <- tDat[,c("datetime","wtr","SpCond","pH", "DOsat", "do.obs", "Chl_conc")] # subset to useful column

	if(i == 1){
		paul12.epi000 <- tDat
	}else{
		paul12.epi000 <- rbind(paul12.epi000, tDat) # accumulate
	}
	
}


# merge sonde data w/ tchain + manual zmix
paul12.epi00 <- merge(paul12.epi000, paul12.zmix, all.x=TRUE)

# a few more hours of z.mix data needed to be added to beginning via constant interpolation rule=2
fill.zmix.ind <- is.na(paul12.epi00[,"z.mix"])
fill.zmix <- approx(x=paul12.epi00[,"datetime"], y=paul12.epi00[,"z.mix"], xout=paul12.epi00[,"datetime"], rule=2, f=0, method="constant")$y

# fill in interp'd zmix at beginning (and, actually, one or two other places too)
paul12.epi0 <- paul12.epi00
paul12.epi0[fill.zmix.ind,"z.mix"] <- fill.zmix[fill.zmix.ind] 

paul12.epi.full0 <- merge(paul12.epi0, irr_wnd_2012, all.x=TRUE) # merge weather and the sonde+tChain data

# add various components needed for metabolism
paul12.epi.full0[,"do.sat"] <- LakeMetabolizer:::o2.at.sat.base(paul12.epi.full0[,"wtr"], 960) # calculate saturated o2
# paul12.epi.full0[,"wnd"] <- LakeMetabolizer:::scale.exp.wind.base(paul12.epi.full0[,"wnd"], 2) # scale wind up to 10m NOTE already done for 2012
paul12.epi.full0[,"k600.cole"] <- LakeMetabolizer:::k.cole.base(paul12.epi.full0[,"wnd"]) # calculate k600
paul12.epi.full0[,"k.gas"] <- LakeMetabolizer:::k600.2.kGAS.base(paul12.epi.full0[,"k600.cole"], paul12.epi.full0[,"wtr"]) # calculate k.gas

paul12.epi.full <- paul12.epi.full0
paul.12.0k <- paul12.epi.full[,"z.mix"] < 0.7
paul12.epi.full[paul.12.0k,"k.gas"] <- 0


paul12.epi <- paul12.epi.full[,c("datetime", "do.obs", "do.sat", "k.gas", "z.mix",  "irr", "wtr", "wnd")]


# =======================
# = Save organized data =
# =======================
save(paul12.epi.full, paul12.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2012.RData")

write.table(paul12.epi.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul12.epi.full.txt", row.names=FALSE)
write.table(paul12.epi, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/paul12.epi.txt", row.names=FALSE)

