rm(list=ls())
graphics.off()

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/WardSensorData2012/")
Data_PF0 <- read.table("Ward_Phytoflash_23Apr2012_to_10May2012.txt", header=TRUE, comment.char="C") #comment.char="C" so that "Configuration=Ward12" doesn't mess up the read
PF_Inval <- which(Data_PF0[,"Yield"]=="Invalid") #Which rows have an "Invalid" in the "Yield" column?
Data_PF0[PF_Inval, c("Fo","Fm","Fv","Yield")] <- NA #Replace measurements when Yield was Invalid with NA
Data_PF0[, "Yield"] <- as.numeric(as.character(Data_PF0[,"Yield"])) #Turn the "Yield" column into numeric (is left as a character from the Invalids)

PF_DateTime <- as.POSIXct(paste(Data_PF0[, "Date"], Data_PF0[, "Time"], sep=" "), format="%m/%d/%y")
PF_Day <- as.numeric(format(PF_DateTime, format="%j"))
PF_Frac <- as.numeric(as.difftime(format(Data_PF0[, "Time"], format="%H:%M:%S"), units="days"))
PF_DoY <- PF_Day+PF_Frac

Data_PF <- as.matrix(x=cbind("DoY"=PF_DoY, Data_PF0[c("Fo","Fm","Fv","Yield")]), ncol=5, dimnames=list(NA, c("DoY","Fo","Fm","Fv","Yield")))

plot(Data_PF[,"DoY"], Data_PF[,"Yield"], type="l")
abline(v=(trunc(min(Data_PF[,"DoY"]))+0.5):125.5, col="lightblue")
abline(v=(trunc(min(Data_PF[,"DoY"]))+0.0):125.0, col="slategray")