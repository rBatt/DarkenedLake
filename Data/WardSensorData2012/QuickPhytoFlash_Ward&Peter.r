rm(list=ls())
graphics.off()


Light <- function(latitude=46.28){

lat <- latitude/57.297 #lat converted to radians at Paul L. dyke, 57.297 = 180/pi
day <- (1:365) #day of year. could be read in from file
rads <- 2*pi*day/365 

#dec is the declination of the sun f of day  Spencer, J.W. 1971: Fourier
#series representation of the position of the Sun. Search, 2(5), 172.
dec <- 0.006918 - 0.399912 * cos(rads) + 0.070257 * sin(rads) - 0.006758 * cos(2*rads) + 0.000907 * sin(2*rads) - 0.00297 * cos(3*rads) + 0.00148 * sin(3*rads)
x <- (-1*sin(lat)*sin(dec))/(cos(lat)*cos(dec))

#hours in radians before true noon of sunrise.
SR <- (pi/2)-atan(x/(sqrt(1-x^2))) 
SR <- SR*2/.262 #converts SR from radians back to hours of daylight.
#need to adjust clock time which is in Central Daylight Savings time by 1 h.

dates <- (1:365)
sunrise <- (12-(SR*0.5)+1)/24 #Sunrise in hours Central Daylight Savings time
sunset <- (12+(SR*0.5)+1)/24   #Sunset in CDST
RiseSet <- matrix(c(sunrise, sunset), nrow=365, ncol=2)
#print("Data stored in 'RiseSet' matrix; 1st col is sunRise, 2nd col is sunSet -- given in fractions of day.")
return(RiseSet)
}


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/WardSensorData2012/")
Data_PF0 <- read.table("Ward_Phytoflash_23Apr2012_to_10May2012.txt", header=TRUE, comment.char="C") #comment.char="C" so that "Configuration=Ward12" doesn't mess up the read
PF_Inval <- which(Data_PF0[,"Yield"]=="Invalid") #Which rows have an "Invalid" in the "Yield" column?
Data_PF0[PF_Inval, c("Fo","Fm","Fv","Yield")] <- NA #Replace measurements when Yield was Invalid with NA
Data_PF0[, "Yield"] <- as.numeric(as.character(Data_PF0[,"Yield"])) #Turn the "Yield" column into numeric (is left as a character from the Invalids)

PF_DateTime <- as.POSIXct(paste(Data_PF0[, "Date"], Data_PF0[, "Time"], sep=" "), format="%m/%d/%y")
PF_Day <- as.numeric(format(PF_DateTime, format="%j"))
PF_Frac <- as.numeric(as.difftime(format(Data_PF0[, "Time"], format="%H:%M:%S"), units="days"))
PF_DoY <- PF_Day+PF_Frac

Data_PF <- as.matrix(x=cbind("DoY"=PF_DoY, Data_PF0[c("Fo","Fm","Fv","Yield")]), ncol=5, dimnames=list(NA, c("DoY","Fo","Fm","Fv","Yield")))[-which(PF_DoY>125),]

dev.new(height=4.25, width=11)
par(mar=c(2,2,0,0), oma=c(2,2,1,0.5))
plot(Data_PF[,"DoY"], Data_PF[,"Yield"], type="l", xlab="", ylab="")
abline(v=(trunc(min(Data_PF[,"DoY"]))+0.5):125.5, col="lightblue")
abline(v=(trunc(min(Data_PF[,"DoY"]))+0.0):125.0, col="slategray")
DaysUsed <- trunc(Data_PF[,"DoY"])
X_LeftRight <- Light()[unique(DaysUsed),]
rect(unique(DaysUsed)+X_LeftRight[,1],0,unique(DaysUsed)+X_LeftRight[,2],1, col="#3A5FCD25")
mtext(paste("Ward DoY,",2012), side=1, line=0.5, outer=TRUE)
mtext("Yield (Fv/Fm)", side=2, line=0.5, outer=TRUE)

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/SquealMetabolism/")
Peter_PF0 <- read.csv("PeterPhytoflash2008to10.csv", header=TRUE)
# R_PF_Inval <- which(Peter_PF0[,"Yield"]=="Invalid")
Peter_PF_Date <- as.POSIXct(paste(as.character(Peter_PF0[, "Date"]), as.character(Peter_PF0[, "Time"], sep=" ")), format="%d-%b-%Y %H:%M:%S")
Peter_PF_Year <- as.numeric(format(Peter_PF_Date, format="%Y"))
Peter_PF_DoY <- as.numeric(format(Peter_PF_Date, format="%j"))
Peter_PF_Frac <- as.numeric(as.difftime(format(Peter_PF_Date, format="%H:%M:%S"), units="days"))

Peter_PF1 <- matrix(c(Peter_PF_Year, (Peter_PF_DoY+Peter_PF_Frac), unlist(Peter_PF0[,c("Fo","Fm","Fv","Yield")])), ncol=6, dimnames=list(NULL,c("Year","DoY", "Fo","Fm","Fv","Yield")))
# Peter_PF <- Peter_PF1[-c(which(diff(Peter_PF1[,"DoY"])<0 & diff(Peter_PF1[, "Year"])==0)),]
Peter_PF <- Peter_PF1[-c(which(duplicated((Peter_PF1[, "Year"]*1000)+Peter_PF1[,"DoY"]))),]
# which(duplicated((Peter_PF1[, "Year"]*1000)+Peter_PF1[,"DoY"]))

SplitN <- 7
PeterYears <- unique(Peter_PF_Year)
Params2Use <- c("Fo", "Fm", "Fv", "Yield")[4]
for(p in 1:length(Params2Use)){
	Param <- Params2Use[p]
	for(ye in 1:length(PeterYears)){
		which(diff(Peter_PF[,"DoY"])<0 & diff(Peter_PF[, "Year"])==0)
		print(which(diff(Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]), "DoY"])<0))
		dev.new(width=11, height=8.5)
		par(mfrow=c(SplitN, 1), mar=c(2,2,0,0), oma=c(2,2,1,0.5))
		Intervals <- cut(1:length(Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]), "DoY"]), SplitN)
		for(i in 1:SplitN){
			plot(Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]),"DoY"][which(Intervals==unique(Intervals)[i])], Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]),Param][which(Intervals==unique(Intervals)[i])], type="l", ylab="", xlab="")
			abline(v=unique(trunc(Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]),"DoY"][which(Intervals==unique(Intervals)[i])])), col="slategray")
			DaysUsed <- trunc(Peter_PF[which(Peter_PF[,"Year"]==PeterYears[ye]),"DoY"][which(Intervals==unique(Intervals)[i])])
			X_LeftRight <- Light()[unique(DaysUsed),]
			rect(unique(DaysUsed)+X_LeftRight[,1],0,unique(DaysUsed)+X_LeftRight[,2],1, col="#3A5FCD25")
		}
		mtext(paste("Peter DoY,",PeterYears[ye]), side=1, line=0.5, outer=TRUE)
		mtext(Param, side=2, line=0.5, outer=TRUE)
	}
}





# graphics.off()


