#_v0 (07-Feb-2013)

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
source("Metabolism_LM_v3.0.R")

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

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulSondes_2012")
Therms <- data.frame("Year"=2012, read.csv("Paul_Tchain_2012.csv", header=TRUE))
ThermNames0 <- names(Therms)
for(i in 3:length(ThermNames0)){
	Therms[,ThermNames0[i]] <- as.numeric(as.character(Therms[,ThermNames0[i]]))
	wrong <- which(Therms[,ThermNames0[i]] <0 | Therms[,ThermNames0[i]] >40)
	Therms[wrong,ThermNames0[i]] <- NA
}

ThermDepths <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.0, 4.5, 5)
zmix <- function(temp, depth){
		miss <- which(is.na(temp))
		temp <- temp[-miss]
		depth <- depth[-miss]
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
ThermsDoY <- as.numeric(format.Date(as.POSIXct(Therms[,2]), format="%j")) + as.numeric(difftime(time1=as.POSIXct(Therms[,2]), time2=as.POSIXct(paste(format.Date(Therms[,2], format="%Y-%m-%d"), "00:00:00"), format="%Y-%m-%d %H:%M:%S"), units="days"))

chainZmix <- data.frame("Year"=2012, "DoY"=ThermsDoY, "Zmix"=Zmixs)
PAR_Wind_Therms <- Alme(X=PAR_Wind, Y=chainZmix, Freq_High_Low="Low", AllX=TRUE)

ManuZmix <- subset(read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv"), Layer=="PML" & Lake=="Paul" & Year==2012, select=c("Year", "Date", "Zmix"))
ManuDoY <- as.numeric(format.Date(as.POSIXct(ManuZmix[,2]), format="%j")) + c(0, rep(0.4,13), 0)
ManuDepths <- ManuZmix[,"Zmix"]
PAR_Wind_Therms[, "Zmix"] <- approx(c(PAR_Wind_Therms[, "DoY"], ManuDoY), c(PAR_Wind_Therms[, "Zmix"], ManuDepths), PAR_Wind_Therms[, "DoY"])$y

plot(ManuDoY, ManuDepths, type="o", lwd=2, ylim=c(4, 0), col="blue", cex=1.5)
points(ThermsDoY, Zmixs, pch=20)
legend("topright", lty=c(1, 0), pch=c(21, 20), col=c("blue", "black"), legend=c("Routine Zmix", "T-chain Zmix"))
# points(PAR_Wind_Therms[,c("DoY","Zmix")], pch=21, col="red")


FileNames <- c("50312L1D.csv", "52312L1D.csv", "61512L1D.csv", "70612L1D.csv", "72712L1D.csv", "81712L1D.csv")
for(i in 1:length(FileNames)){
	tName <- FileNames[i]
	ReadClassTxt <- c(NA, NA, rep("numeric", 10))
	tDat <- read.table(paste("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulSondes_2012/", tName, sep=""), sep=",", header=TRUE, skip=4, colClasses=ReadClassTxt)#[-c(1,2,3,4),]
	names(tDat) <- c("Date", "Time", "Temp", "SpCond", "pH", "BGA_conc", "BGA_RFU", "DOsat", "DOconc", "Chl_conc", "Chl_RFU", "Battery")
	# tDay <- as.numeric(format.Date(as.POSIXct(tDat[,1], format="%Y/%m/%d"), format="%j"))
	# tFrac <- as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%Y/%m/%d %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%Y/%m/%d %H:%M:%S"), units="days"))
	tDoY <- as.numeric(format.Date(as.POSIXct(tDat[,1]), format="%j")) + as.numeric(difftime(time1=as.POSIXct(paste(tDat[,1], tDat[,2]), format="%Y-%m-%d %H:%M:%S"), time2=as.POSIXct(paste(tDat[,1], "00:00:00"), format="%Y-%m-%d %H:%M:%S"), units="days"))
	tDat <- tDat[,-9]
	
	# FirstLastDays <- which(is.element(floor(tDoY), range(floor(tDoY))))
	# tDat <- tDat[-FirstLastDays,]
	# tDoY <- tDoY[-FirstLastDays]

	if(i == 1){
		DataSonde <- data.frame("Year"=2012, "DoY"=tDoY, tDat[,-c(1:2)])
	}else{
		DataSonde <- rbind(DataSonde, data.frame("Year"=2012, "DoY"=tDoY, tDat[,-c(1:2)]))
	}
	
}

sDataSonde <- ByeShort(DataSonde)
DataAll0 <- Alme(X=sDataSonde, Y=PAR_Wind_Therms)
DataAll <- cbind(DataAll0, "DepID"=Chunks(DataAll0[,"DoY"]))



Data <- DataAll[,c("Year", "DoY", "DOsat", "Temp", "PAR", "Zmix", "Wind")]
Est_LM <- Metabolism_LM(Data)



Deps <- DataAll[,"DepID"]
for(i in 1:length(unique(Deps))){
	if(length(which(Deps==i))<288*3){next}
	print(i)
	INDEX <- which(Deps==i)
	Data01 <- Data[INDEX,]
	Sonde <- data.frame("Year"=Data01[,"Year"], "DoY"=trunc(Data01[,"DoY"]), "Fract"=(Data01[,"DoY"]-trunc(Data01[,"DoY"])), "Temp"=Data01[,"Temp"], "DOsat"=Data01[,"DOsat"], DepID=rep(i, length(Data01[,"DOsat"])))
	
	
	if(i==1){
		Est_BK <- Metabolism(Sonde, zmix=Data01[,"Zmix"], Wind=Data01[,"Wind"], DispGraph=FALSE, GraphNew=FALSE)
	}
	if(i!=1){
		Est_BK <- rbind(Est_BK, Metabolism(Sonde, zmix=Data01[,"Zmix"], Wind=Data01[,"Wind"], DispGraph=FALSE, GraphNew=FALSE))
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

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data")
load("Paul_BK&Data_2010.RData")
Bad_BK_2010 <- union(which(L10_BK[,"GPP"]<0), which(L10_BK[,"R"]>0))
Est_Good_BK_2010 <- L10_BK[-Bad_BK_2010,]


PaulData2010 <- L10Sonde
L10Sonde[,"DoY"] <- L10Sonde[,"DoY"]+L10Sonde[,"Fract"]
L10Sonde[,"pH"] <- L10pH
L10Sonde[,"Chl_conc"] <- L10Chl
# PaulID10 <- Chunks(xout10, Sub=288)
PaulWeek10 <- as.numeric(format.Date(as.POSIXct(paste(2010, xout10, sep="-"), format="%Y-%j"), format="%m"))-4
PaulData2010 <- data.frame(L10Sonde[,c("Year", "DoY", "Temp", "pH", "DOsat", "Chl_conc")], "DepID"=1, "Week"=PaulWeek10)
Paul_Tchain_Weather_2010 <- data.frame("Year"=2010, "DoY"=xout10, "PAR"=PAR10, "Wind"=Wind10, "Zmix"=Lzmix10)



PaulID12 <- Chunks(sDataSonde[,"DoY"])
PaulWeek12 <- as.numeric(format.Date(as.POSIXct(paste(2012, sDataSonde[,"DoY"], sep="-"), format="%Y-%j"), format="%m"))-4
PaulData2012 <- data.frame(sDataSonde[,c("Year", "DoY", "Temp", "pH", "DOsat", "Chl_conc")], "DepID"=PaulID12, "Week"=PaulWeek12)
Paul_Tchain_Weather_2012 <- PAR_Wind_Therms





Paul_Year <- c(rep(2010, nrow(Est_Good_BK_2010)), rep(2012, nrow(Est_Good_BK)))
Paul_DoY <- c(Est_Good_BK_2010[,"DoY"], Est_Good_BK[,"DoY"])
Paul_Week <- as.numeric(format.Date(as.POSIXct(paste(Paul_Year, Paul_DoY, sep="-"), format="%Y-%j"), format="%m"))-4 #This corresponds the the "isotope sampling week"; in 2010, the first week of sampling May, 2nd in June, 3 in July, 4th in August.  In 2012, the first sampling week was in June, but I have been using the convention of "sampling week = month # -4".  So August is the 8th month.  A sample in August corresponds to the 8-4 = 4th sampling week.  Totally weird, I know.  In the figures this will all be turned into the month.
Paul_Method <- "BK"
Paul_GPP <- c(Est_Good_BK_2010[,"GPP"], Est_Good_BK[,"GPP_BK"])
Paul_R <- c(Est_Good_BK_2010[,"R"], Est_Good_BK[,"R_BK"])
Paul_NEP <- c(Est_Good_BK_2010[,"NEP"], Est_Good_BK[,"NEP_BK"])
PaulMetabolism <- data.frame("Year"=Paul_Year, "DoY"=Paul_DoY, "Week"=Paul_Week, "Method"=Paul_Method, "GPP"=Paul_GPP, "R"=Paul_R, "NEP"=Paul_NEP)





# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data")
# save(PaulData2010, PaulData2012, Paul_Tchain_Weather_2010, Paul_Tchain_Weather_2012, PaulMetabolism, file="Paul_Metab&Data_2010&2012.RData")