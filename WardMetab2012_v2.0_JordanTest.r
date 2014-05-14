#_v2.0 (14-Feb-2013) I need to reorganize everything so badly.  For right now, I want to read in the metalimnion data from 2012.  If I do this, then I'll have all the data analyzed... even if it isn't the exact same method or in the exact same place.  Right now there are 2 lakes with 2 years of epi data, and 1 lake with 2 years of meta data.  The metabolism estimates are in about 4 or 5 different places right now.


rm(list=ls())
graphics.off()
library("zoo")
library("plyr")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/SquealMetabolism")
source("MyBookkeepingMetabolism.R")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/PAR_DO_PhaseShift/")
source("Alme.R")
source("ByeShort.R")
source("Chunks.R")
source("Light.R")
source("Metabolism_LM_v3.1.R")

windowsFonts(Times=windowsFont("TT Times New Roman"))

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
Therms <- data.frame("Year"=2012, read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Data/CompiledWardTherms2012_v2.txt", header=TRUE))
ThermDepths <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3)

zmix <- function(temp, depth){
        tempDiff <- temp[-length(temp)] - temp[-1]
        depthDiff <- depth[-1] - depth[-length(depth)]
        ratio <- tempDiff / depthDiff
		zhat <- depth[which(ratio >= 2)][1]
		# Zhat2 <- min(which(zhat>0))
		# Z <- Zhat2
		Z <- ifelse(!(zhat>0), NA, zhat)
       return(Z)
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
PAR_Wind_Therms[c(17004:19330, 12454:12538), 5:16] <- NA
plot(1:length(PAR_Wind_Therms[,"z0.00"]), PAR_Wind_Therms[,"z0.00"])
abline(v=c(17004, 19330), col="red")
abline(v=c(12454, 12538), col="green")
# abline(h=c(27, 33), col="red")

badtherms <- which(PAR_Wind_Therms[,5:15] < 5 | PAR_Wind_Therms[,5:15] > 36, arr.ind=TRUE)
PAR_Wind_Therms[,5:15][PAR_Wind_Therms[,5:15] < 5 | PAR_Wind_Therms[,5:15] > 36]

# Therms2 <- PAR_Wind_Therms[,5:15]
# which(PAR_Wind_Therms[,"DoY"]
# image.plot(x=PAR_Wind_Therms[,"DoY"], y=ThermDepths, z=as.matrix(Therms2), zlim=c(0,40), ylim=c(3,0))

# badtherms <- which(Therms2 < 5 | Therms2 > 36, arr.ind=TRUE)
# Therms2[badtherms] <- 500

ZmixDays <- c(105, 110, 117, 124, 131, 138, 145, 152, 159, 166, 173, 180, 187, 194, 201, 208, 215, 222, 229, 236, 241) + c(0, rep(0.4,19), 0)
ZmixDepths <- c(1.0,2.0,2.0,0.5,0.5, 1.0, 1.0, 1.5, 0.25, 1.0, 1.0, 0.5, 0.5, 1, 1, 1, 1, 1.5, 1.5, 1, 0.25)

PAR_Wind_Therms[, "Zmix"] <- approx(c(PAR_Wind_Therms[, "DoY"], ZmixDays), c(PAR_Wind_Therms[, "Zmix"], ZmixDepths), PAR_Wind_Therms[, "DoY"])$y
# plot(PAR_Wind_Therms[, "DoY"], PAR_Wind_Therms[, "Zmix"], type="l")

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/")
Key2012 <- read.csv("WardSondes2012_FileKey.csv")
ShalKey2012 <- subset(Key2012, Depth=="B")


for(i in 1:nrow(ShalKey2012)){
	tName <- ShalKey2012[i,"Name"]
	tFormat <- ShalKey2012[i,"Format"]
	
	if(tFormat=="CDF"){
		if(tName=="B14APR12.CDF"){
			ReadClassCDF <- c("character", "character", rep("numeric", 8))
			tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU")
			tDoY <- as.numeric(format.Date(as.POSIXct(as.character(tDat[,1]), format="%m/%d/%y"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%m/%d/%y %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%m/%d/%y %H:%M:%S"), units="days"))
		}
		if(all(tName!="B14APR12.CDF" & ShalKey2012[,"Depth"]=="B")){
			ReadClassCDF <- c("character", "character", rep("numeric", 9))
			tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
			names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "Battery")
			tDoY <- as.numeric(format.Date(as.POSIXct(as.character(tDat[,1]), format="%m/%d/%y"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%m/%d/%y %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%m/%d/%y %H:%M:%S"), units="days"))
		}
		# ReadClassCDF <- c("character", "character", rep("numeric", 9))
		# tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=TRUE, skip=2, colClasses=ReadClassCDF)
		# names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "Chl_conc", "Chl_RFU", "DOconc")
		# tDoY <- as.numeric(format.Date(as.POSIXct(as.character(tDat[,1]), format="%m/%d/%y"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%m/%d/%y %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%m/%d/%y %H:%M:%S"), units="days"))
	}
	if(tFormat=="txt"){
		ReadClassTxt <- c(NA, NA, rep("numeric", 10))
		tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/WardSondes2012/", tName, sep=""), sep=",", header=FALSE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
		names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "DOconc", "Chl_conc", "Chl_RFU", "Battery")
		# tDay <- as.numeric(format.Date(as.POSIXct(tDat[,1], format="%Y/%m/%d"), format="%j"))
		# tFrac <- as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%Y/%m/%d %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%Y/%m/%d %H:%M:%S"), units="days"))
		tDoY <- as.numeric(format.Date(as.POSIXct(tDat[,1], format="%Y/%m/%d"), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%Y/%m/%d %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%Y/%m/%d %H:%M:%S"), units="days"))
		tDat <- tDat[,-9]
	}
	tDat <- tDat[,c("Date", "Time", "Temp", "SpCond", "pH", "DOsat", "Chl_conc")]
	
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



DataAll0 <- Alme(X=sDataSonde, Y=PAR_Wind_Therms)
DataAll <- cbind(DataAll0, "DepID"=Chunks(DataAll0[,"DoY"]))




Data <- DataAll[,c("Year", "DoY", "DOsat", "Temp", "PAR", "Zmix", "Wind")]
Est_LM <- Metabolism_LM(Data)
Metabolism_LM(Data[1:5000,])
jData <- Data[1:5000,]
save(jData, file="jData.RData")
write.csv(jData, file="jData.csv", row.names=FALSE)

jData2 <- read.csv("jData.csv")
Metabolism_LM(jData2)



Deps <- DataAll[,"DepID"]
for(i in 1:length(unique(Deps))){
	if(length(which(Deps==i))<288*3){next}
	print(i)
	INDEX <- which(Deps==i)
	Data01 <- Data[INDEX,]
	Sonde <- data.frame("Year"=Data01[,"Year"], "DoY"=trunc(Data01[,"DoY"]), "Fract"=(Data01[,"DoY"]-trunc(Data01[,"DoY"])), "Temp"=Data01[,"Temp"], "DOsat"=Data01[,"DOsat"], DepID=rep(i, length(Data01[,"DOsat"])))
	
	
	if(i==1){
		Est_BK <- Metabolism(Sonde, zmix=Data01[,"Zmix"], Wind=Data01[,"Wind"], DispGraph=FALSE)
	}
	if(i!=1){
		Est_BK <- rbind(Est_BK, Metabolism(Sonde, zmix=Data01[,"Zmix"], Wind=Data01[,"Wind"], DispGraph=FALSE))
	}	
}


Bad_LM <- union(which(Est_LM[,"GPP_raw"]<0), which(Est_LM[,"R_raw"]>0))
Est_Good_LM <- Est_LM[-Bad_LM,]
names(Est_Good_LM) <- c("Year", "DoY", "GPP_LM", "R_LM", "NEP_LM", "sumPAR", "meanTemp", "TotalF", "R2_LM")

Bad_BK <- union(which(Est_BK[,"GPP"]<0), which(Est_BK[,"R"]>0))
Est_Good_BK <- Est_BK[-Bad_BK,]
dimnames(Est_Good_BK) <- list(NULL, c("DoY", "R_BK", "GPP_BK", "NEP_BK"))


Estimates_Good <- merge(Est_Good_BK, Est_Good_LM, by="DoY", all=TRUE)[,c("Year", "DoY", "GPP_BK","R_BK","NEP_BK", "GPP_LM", "R_LM", "NEP_LM")]
Estimates_Good[,"Year"] <- 2012


Estimates_Good_NoNA <- Estimates_Good[complete.cases(Estimates_Good),]


# ==========================================
# = Plot the "good" BK estimates from 2012 =
# ==========================================
dev.new(height=4, width=4)
par(mar=c(4,3,0,0), oma=c(0,0,0,0), family="Times", las=1)
Metab_Ylim <- c(min(Est_Good_BK[,"R_BK"]), max(Est_Good_BK[,"GPP_BK"]))
plot(Est_Good_BK[,"DoY"], Est_Good_BK[,"NEP_BK"], type="o", lty="dashed", col="black", pch=23, bg=gray(0.5), ylim=Metab_Ylim, bty="l", xlab="", ylab="", cex=0.75, cex.axis=0.75)
lines(Est_Good_BK[,"DoY"], Est_Good_BK[,"GPP_BK"], type="o", lty="solid", col="black", pch=24, bg=gray(0), cex=0.75)
lines(Est_Good_BK[,"DoY"], Est_Good_BK[,"R_BK"], type="o", lty="dotted", col="black", pch=25, bg=gray(1), cex=0.75)
legend(x=120, y=150, legend=c(paste("GPP", round(mean(Est_Good_BK[,"GPP_BK"]),0), sep=" = "), paste("NEP", round(mean(Est_Good_BK[,"NEP_BK"]),0), sep=" = "), paste("R", round(mean(Est_Good_BK[,"R_BK"]),0), sep=" = ")), lty=c("solid", "dashed", "dotted"), pch=c(24,23,25), pt.bg=c(gray(0), gray(0.5), gray(1)), bty="n", cex=0.75) #, title=expression(underline(Summer~Means))
mtext(side=2, line=2, text=expression(mu*mol~O[2]~L^-1~day^-1), las=0, cex=0.75)
mtext(side=1, line=1.9, text="Day of Year", cex=0.75)
mtext(side=1, line=2.6, text=paste(as.Date(min(Est_Good_BK[,"DoY"]), origin=paste(2012,"-01-01",sep="")), as.Date(max(Est_Good_BK[,"DoY"]), origin=paste(2012,"-01-01",sep="")), sep=" to "), font=3, cex=0.75)


# ==========================================
# = Plot the "good" LM estimates from 2012 =
# ==========================================
dev.new(height=4, width=4)
par(mar=c(4,3,0,0), oma=c(0,0,0,0), family="Times", las=1)
Metab_Ylim <- c(min(Est_Good_LM[,"R_LM"]), max(Est_Good_LM[,"GPP_LM"]))
plot(Est_Good_LM[,"DoY"], Est_Good_LM[,"NEP_LM"], type="o", lty="dashed", col="black", pch=23, bg=gray(0.5), ylim=Metab_Ylim, bty="l", xlab="", ylab="", cex=0.75, cex.axis=0.75)
lines(Est_Good_LM[,"DoY"], Est_Good_LM[,"GPP_LM"], type="o", lty="solid", col="black", pch=24, bg=gray(0), cex=0.75)
lines(Est_Good_LM[,"DoY"], Est_Good_LM[,"R_LM"], type="o", lty="dotted", col="black", pch=25, bg=gray(1), cex=0.75)
legend(x=120, y=150, legend=c(paste("GPP", round(mean(Est_Good_LM[,"GPP_LM"]),0), sep=" = "), paste("NEP", round(mean(Est_Good_LM[,"NEP_LM"]),0), sep=" = "), paste("R", round(mean(Est_Good_LM[,"R_LM"]),0), sep=" = ")), lty=c("solid", "dashed", "dotted"), pch=c(24,23,25), pt.bg=c(gray(0), gray(0.5), gray(1)), bty="n", cex=0.75) #, title=expression(underline(Summer~Means))
mtext(side=2, line=2, text=expression(mu*mol~O[2]~L^-1~day^-1), las=0, cex=0.75)
mtext(side=1, line=1.9, text="Day of Year", cex=0.75)
mtext(side=1, line=2.6, text=paste(as.Date(min(Est_Good_LM[,"DoY"]), origin=paste(2012,"-01-01",sep="")), as.Date(max(Est_Good_LM[,"DoY"]), origin=paste(2012,"-01-01",sep="")), sep=" to "), font=3, cex=0.75)



# ==============================================
# = Time to look into the 2010 metabolism data =
# ==============================================
#See setwd("/Users/Battrd/Documents/School&Work/WiscResearch/KalmanFilter") for original (KFMetabolismSingle_v3.R)
#See setwd("/Users/Battrd/Documents/School&Work/GradSchool/DissertationProposal/") for what I used in my dissertation proposal (pretty much the same as the original)
load("/Users/Battrd/Documents/School&Work/GradSchool/DissertationProposal/Ward2010_Epi_Metabolism_Data.RData")
Est_BK_2010 <- EpiMetab
Bad_BK_2010 <- union(which(Est_BK_2010[,"GPP"]<0), which(Est_BK_2010[,"R"]>0))
Est_Good_BK_2010 <- Est_BK_2010[-Bad_BK_2010,]
dimnames(Est_Good_BK_2010) <- list(NULL, c("DoY", "R_BK", "GPP_BK", "NEP_BK"))
Est_Good_BK_2010 <- cbind("Year"=2010, Est_Good_BK_2010)
colMeans(Est_Good_BK_2010)

#Working on implementing the LM method for the 2010 data (5-Dec-2012... on my bday!)
Data_2010 <- data.frame("Year"=2010, "DoY"=Sonde1[,"DoY"]+Sonde1[,"Fract"], "DOsat"=Sonde1[,"DOsat"], "Temp"=Sonde1[,"Temp"], "PAR"=PAR1, "Zmix"=Zmix1, "Wind"=Wind1)
Est_LM_2010 <- Metabolism_LM(Data_2010, ExpectedFreq=360)
Bad_LM_2010 <- union(which(Est_LM_2010[,"GPP_raw"]<0), which(Est_LM_2010[,"R_raw"]>0))
Est_Good_LM_2010 <- Est_LM_2010[-Bad_LM_2010,]
names(Est_Good_LM_2010) <- c("Year", "DoY", "GPP_LM", "R_LM", "NEP_LM", "sumPAR", "meanTemp", "TotalF", "R2_LM")

# ==========================================
# = Plot the "good" BK estimates from 2010 =
# ==========================================
dev.new(height=4, width=4)
par(mar=c(4,3,0,0), oma=c(0,0,0,0), family="Times", las=1)
Metab_Ylim <- c(min(Est_Good_BK_2010[,"R_BK"]), max(Est_Good_BK_2010[,"GPP_BK"]))
plot(Est_Good_BK_2010[,"DoY"], Est_Good_BK_2010[,"NEP_BK"], type="o", lty="dashed", col="black", pch=23, bg=gray(0.5), ylim=Metab_Ylim, bty="l", xlab="", ylab="", cex=0.75, cex.axis=0.75)
lines(Est_Good_BK_2010[,"DoY"], Est_Good_BK_2010[,"GPP_BK"], type="o", lty="solid", col="black", pch=24, bg=gray(0), cex=0.75)
lines(Est_Good_BK_2010[,"DoY"], Est_Good_BK_2010[,"R_BK"], type="o", lty="dotted", col="black", pch=25, bg=gray(1), cex=0.75)
legend(x=194, y=89, legend=c(paste("GPP", round(mean(Est_Good_BK_2010[,"GPP_BK"]),0), sep=" = "), paste("NEP", round(mean(Est_Good_BK_2010[,"NEP_BK"]),0), sep=" = "), paste("R", round(mean(Est_Good_BK_2010[,"R_BK"]),0), sep=" = ")), lty=c("solid", "dashed", "dotted"), pch=c(24,23,25), pt.bg=c(gray(0), gray(0.5), gray(1)), bty="n", cex=0.75) #, title=expression(underline(Summer~Means))
mtext(side=2, line=2, text=expression(mu*mol~O[2]~L^-1~day^-1), las=0, cex=0.75)
mtext(side=1, line=1.9, text="Day of Year", cex=0.75)
mtext(side=1, line=2.6, text=paste(as.Date(min(Est_Good_BK_2010[,"DoY"]), origin=paste(2010,"-01-01",sep="")), as.Date(max(Est_Good_BK_2010[,"DoY"]), origin=paste(2010,"-01-01",sep="")), sep=" to "), font=3, cex=0.75)


#Plot the "good" LM estimates from 2010
dev.new(height=4, width=4)
par(mar=c(4,3,0,0), oma=c(0,0,0,0), family="Times", las=1)
Metab_Ylim <- c(min(Est_Good_LM_2010[,"R_LM"]), max(Est_Good_LM_2010[,"GPP_LM"]))
plot(Est_Good_LM_2010[,"DoY"], Est_Good_LM_2010[,"NEP_LM"], type="o", lty="dashed", col="black", pch=23, bg=gray(0.5), ylim=Metab_Ylim, bty="l", xlab="", ylab="", cex=0.75, cex.axis=0.75)
lines(Est_Good_LM_2010[,"DoY"], Est_Good_LM_2010[,"GPP_LM"], type="o", lty="solid", col="black", pch=24, bg=gray(0), cex=0.75)
lines(Est_Good_LM_2010[,"DoY"], Est_Good_LM_2010[,"R_LM"], type="o", lty="dotted", col="black", pch=25, bg=gray(1), cex=0.75)
legend(x=140, y=90, legend=c(paste("GPP", round(mean(Est_Good_LM_2010[,"GPP_LM"]),0), sep=" = "), paste("NEP", round(mean(Est_Good_LM_2010[,"NEP_LM"]),0), sep=" = "), paste("R", round(mean(Est_Good_LM_2010[,"R_LM"]),0), sep=" = ")), lty=c("solid", "dashed", "dotted"), pch=c(24,23,25), pt.bg=c(gray(0), gray(0.5), gray(1)), bty="n", cex=0.75) #, title=expression(underline(Summer~Means))
mtext(side=2, line=2, text=expression(mu*mol~O[2]~L^-1~day^-1), las=0, cex=0.75)
mtext(side=1, line=1.9, text="Day of Year", cex=0.75)
mtext(side=1, line=2.6, text=paste(as.Date(min(Est_Good_LM_2010[,"DoY"]), origin=paste(2010,"-01-01",sep="")), as.Date(max(Est_Good_LM_2010[,"DoY"]), origin=paste(2010,"-01-01",sep="")), sep=" to "), font=3, cex=0.75)


colMeans(Est_Good_BK_2010)
colMeans(subset(Est_Good_BK, "DoY">=138))
colMeans(subset(Est_Good_LM, "DoY">=138))[c("Year", "DoY", "GPP_LM", "R_LM", "NEP_LM")]

Ward_Year <- c(rep(2010, (nrow(Est_Good_BK_2010)+nrow(Est_Good_LM_2010))), rep(2012, (nrow(Est_Good_BK)+nrow(Est_Good_LM))))
Ward_DoY <- c(Est_Good_BK_2010[,"DoY"], Est_Good_LM_2010[,"DoY"], Est_Good_BK[,"DoY"], Est_Good_LM[,"DoY"])
Ward_Week <- as.numeric(format.Date(as.POSIXct(paste(Ward_Year, Ward_DoY, sep="-"), format="%Y-%j"), format="%m"))-4 #This corresponds the the "isotope sampling week"; in 2010, the first week of sampling May, 2nd in June, 3 in July, 4th in August.  In 2012, the first sampling week was in June, but I have been using the convention of "sampling week = month # -4".  So August is the 8th month.  A sample in August corresponds to the 8-4 = 4th sampling week.  Totally weird, I know.  In the figures this will all be turned into the month.
Ward_Method <- c(rep("BK", nrow(Est_Good_BK_2010)), rep("LM", nrow(Est_Good_LM_2010)), rep("BK", nrow(Est_Good_BK)), rep("LM", nrow(Est_Good_LM)))
Ward_GPP <- c(Est_Good_BK_2010[,"GPP_BK"], Est_Good_LM_2010[,"GPP_LM"], Est_Good_BK[,"GPP_BK"], Est_Good_LM[,"GPP_LM"])
Ward_R <- c(Est_Good_BK_2010[,"R_BK"], Est_Good_LM_2010[,"R_LM"], Est_Good_BK[,"R_BK"], Est_Good_LM[,"R_LM"])
Ward_NEP <- c(Est_Good_BK_2010[,"NEP_BK"], Est_Good_LM_2010[,"NEP_LM"], Est_Good_BK[,"NEP_BK"], Est_Good_LM[,"NEP_LM"])
WardMetabolism <- data.frame("Year"=Ward_Year, "DoY"=Ward_DoY, "Week"=Ward_Week, "Method"=Ward_Method, "GPP"=Ward_GPP, "R"=Ward_R, "NEP"=Ward_NEP)
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
# save(WardMetabolism, file="WardMetabolism_WardMetab2012_v1.RData")

Chla0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/Ward_Chla_2010&2012.csv")
Chla_DoY <- as.numeric(format.Date(Chla0[,"Date"], "%j"))
Chla_Week <- as.numeric(format.Date(as.POSIXct(paste(Chla0[,"Year"], Chla_DoY, sep="-"), format="%Y-%j"), format="%m"))-4
Chla <- data.frame("Year"=Chla0[,"Year"], "Week"=Chla_Week, Chla0[,c(2,3,4,5)])
# save(Chla, file="Chla_WardMetab2012_v1.RData")

Photic0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/Ward_Photic_2010&2012.csv")
Photic_DoY <- as.numeric(format.Date(Photic0[,"Date"], "%j"))
Photic_Week <- as.numeric(format.Date(as.POSIXct(paste(Photic0[,"Year"], Photic_DoY, sep="-"), format="%Y-%j"), format="%m"))-4
Photic <- data.frame("Year"=Photic0[,"Year"], "Week"=Photic_Week, Photic0[,c(2,3,4,5)])
# save(Photic, file="Photic_WardMetab2012_v1.RData")


DOM0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/Ward_DOM_2010&2012.csv")
DOM_DoY <- as.numeric(format.Date(DOM0[,"Date"], "%j"))
DOM_Week <- as.numeric(format.Date(as.POSIXct(paste(DOM0[,"Year"], DOM_DoY, sep="-"), format="%Y-%j"), format="%m"))-4
DOM <- data.frame("Year"=DOM0[,"Year"], "Week"=DOM_Week, DOM0[,c(2,4,5)])
# save(DOM, file="DOM_WardMetab2012_v1.RData")

Metabolism_Plots <- expand.grid(c("GPP", "R", "NEP"), c("BK", "LM"))
dev.new(width=9, height=6)
par(mfrow=c(2,3), mar=c(2,3.5,1,1), ps=10, cex=1, oma=c(0,1,2,0))

for(i in 1:nrow(Metabolism_Plots)){
	MetabD <- subset(WardMetabolism, Method==as.character(Metabolism_Plots[i,2]))
	
	MetabMonths <- c("Apr","May","Jun","Jul","Aug")[1+expand.grid(sort(unique(MetabD[,"Week"])), sort(unique(MetabD[,"Year"])))[,1]]
	RepYearCol <- length(MetabMonths)/2
	
	Response <- MetabD[,as.character(Metabolism_Plots[i,1])]
	Pred_Week <- as.numeric(MetabD[,"Week"])
	Pred_Year <- as.factor(MetabD[,"Year"])
	# dev.new()
	boxplot(Response~Pred_Week*Pred_Year, data=MetabD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=MetabMonths)
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	# mtext(as.character(Metabolism_Plots[i,]), side=2, line=2.5)
	if(i==1){mtext("BK", side=2, line=2.5)}
	if(i==1){mtext("GPP", side=3, line=1)}
	if(i==2){mtext("R", side=3, line=1)}
	if(i==3){mtext("NEP", side=3, line=1)}
	if(i==4){mtext("LM", side=2, line=2.5)}
	if(i==2){legend("bottomleft", c("Reference (2010)", "Aquashade (2012)"), text.col=c("red", "blue"), bty="n")}
}

colMeans(WardMetabolism[,-c(1:4)])

summary(lm(GPP~as.factor(Year)*Method + Method*relevel(as.factor(Week),5) +relevel(as.factor(Week),5)*as.factor(Year), data=WardMetabolism))
#The grand mean of gross primary production was 35 µmol O2/ L / d. The volumetric rate of GPP was not significantly different between years (p=0.6853). May GPP was lower in 2012 than in 2010 (14 µmol O2/ L / d lower, p=0.0366).  On the other hand, there were stronger month-to-month differences in GPP in 2012 than in 2010.  For example, July GPP was high in 2012, and May GPP was low in 2012. On average, there was no detectable difference between GPP estimates from the two methods (linear model [LM] and bookkeeping [BK]) (p=0.3692). However, the seasonal pattern in GPP was different between methods. April metabolism was only measured in 2012. Despite the models being generally similar, the April 2012 GPP estimate was lower for the LM method than it was for BK method (18 µmol O2/ L / d lower, p=0.0580). If such month-specific differences between the methods are accounted for, then LM is shown to have slightly higher estimates of GPP in 2012 than BK in 2012 (8 µmol O2/ L / d, p=0.0717).  The adjusted R^2 for this model was 0.2291.

summary(lm(R~as.factor(Year)*Method + Method*relevel(as.factor(Week),5) +relevel(as.factor(Week),5)*as.factor(Year), data=WardMetabolism))
#The grand mean of respiration was 39 µmol O2/ L / d.  The volumetric rate of R was significantly greater in 2012 than it was in 2010 (25 µmol O2/ L / d greater, p=2.63e-06).  Overall, our data did not reveal a difference between the two methods' estimates of R (p=0.313). However, 2012 R estimated by LM was significantly lower than 2012 R estimated by BK (23 µmol O2/ L / d less respiration for LM 2012 than for BK 2012, p=4.56e-06).  This effectively means that both methods estimated greater R in 2012 than in 2010, but the difference between years was much greater for BK (25) than for LM (25-23=2). In general, respiration rates tended to increase over the season (low R in April, high R in August). The LM method estimated much greater rates of April 2012 R than did the BK method (32 µmol O2/ L / d, p=0.001808).  June and July R was much greater in 2012 than in 2010 (24 and 25 µmol O2/ L / d, respectively; p=0.000413 and p=0.000123, respectively).  The adjusted R^2 for this model was 0.4149.


summary(lm(NEP~as.factor(Year)*Method + Method*relevel(as.factor(Week),5) +relevel(as.factor(Week),5)*as.factor(Year), data=WardMetabolism))
#The grand mean of net ecosystem production was -4 µmol O2/ L / d.  In summer, BK estimates that 2012 NEP is much lower than 2010 NEP (27 µmol O2/ L / d lower, p=2.84e-12), while LM adjusts this estimate upwards by 31 µmol O2/ L / d (p<2e-16), such that NEP is ~0 (the intercept was -4 µmol O2/ L / d).  In addition to the overall balance between the years, BK also shows large differences between the monthly rates of NEP, especially in 2010 (in 2012 BK estimates April NEP as positive, while the other weeks are negative).  The adjusted R^2 for this model was 0.381.






Ward_Data_2012 <- sDataSonde
Ward_Data_2010 <- Data_2010
Ward_Tchain_Weather_2012 <- PAR_Wind_Therms


# dev.new()
# plot(sDataSonde[,c("DoY","Temp")], type="l", lwd=2, col="blue")
# lines(Data_2010[,c("DoY", "Temp")], type="l", lwd=2, col="red")

# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data")
# save(Ward_Data_2012, Ward_Data_2010, WardMetabolism, Ward_Tchain_Weather_2012, file="Ward_Metab&Data_2010&2012.RData")
