#Ryan Batt
#13-July-2012
#Assess the concentration of aquashade in ward lake water using spectral absorbance
#Version 3 (19-Aug-2013): just adding the last date of aquashade additions to the AquashadeAdditions object.

rm(list=ls())
graphics.off()

# =======================================================================
# = Set Options for File Name and Range of Lambda to Use as a Reference =
# =======================================================================
FileName <- c(paste("InSeasonWeeklyReadings/","WardFiltrate_",format.Date(Sys.Date(), "%d%b%Y"), ".csv", sep=""))#, "InSeasonWeeklyReadings/1000ppmWardPMLFiltrate22Apr12.csv")[2] #Use this default if you want the file from the current day
# FileName <- "InSeasonWeeklyReadings/WardFiltrate_12Jul2012.csv"
StdLambda <- list(620:630,300:900,400:700, 635)[[1]]

AquashadeAdditions <- data.frame("Date"=as.POSIXct(c("23-Apr-2012", "24-Apr-2012", "26-Apr-2012",  "24-May-2012", "11-Jun-2012", "14-July-2013"), format="%d-%b-%Y"), "Gal"=c(9, 5, 10, 1, 4, 6))

# =============================
# = Set the Working Directory =
# =============================
BeginDirectory <- paste("/Users/",Sys.info()["user"], sep="")
setwd(paste(BeginDirectory,"/TheCascadeProject/WardSpecScans2012", sep=""))

# =======================================================================
# = Read in the Data ('Standards' and Today's Ward Absorbance Spectrum) =
# =======================================================================
MQ2ppmStandard_0 <- read.csv("PreSeasonStandards/2ppmAquashadeMQ23Apr2012.csv", header=TRUE, skip=1)[,1:2]
PreWardH2O_0 <- read.csv("PreSeasonStandards/WardFiltrate23Apr2012.csv", header=TRUE, skip=1)[,1:2]

# ===============================================
# = Clip the Spectra to the Desired Wavelengths =
# ===============================================
StdIndex <- is.element(MQ2ppmStandard_0[,1], StdLambda)
MQ2ppmStandard <- MQ2ppmStandard_0[StdIndex,]
PreWardH2O <- PreWardH2O_0[StdIndex,]
PreLong <- mean(PreWardH2O_0[401:601,2])

# ===============================================================================================
# = TheoreticalWavelengths/TheoreticalConcentration = ObservedWavelengths/ObservedConcentration =
# ===============================================================================================

# list.files("InSeasonWeeklyReadings/")
# file.info(paste("InSeasonWeeklyReadings/",list.files("InSeasonWeeklyReadings/"), sep=""))
ListFileDates_00 <- unlist(strsplit(list.files("InSeasonWeeklyReadings/"), split="WardFiltrate_"))
ListFileDates_0 <- ListFileDates_00[-c(seq(1, (length(ListFileDates_00)-1), by=2))]
ListFileDates <- as.POSIXct(unlist(strsplit(ListFileDates_0, split=".csv")), format="%d%b%Y")

OrderedFileDates <- ListFileDates[order(ListFileDates)]
OrderedFiles <- list.files("InSeasonWeeklyReadings/")[order(ListFileDates)]
EstAquashadeConc <- rep(NA, length(OrderedFiles))
SDAquashadeConc <- rep(NA, length(OrderedFiles))
for(i in 1:length(EstAquashadeConc)){
	ThisWardH2O_0 <- read.csv(paste("InSeasonWeeklyReadings/",OrderedFiles[i],sep=""), header=TRUE, skip=1)[,1:2]
	ThisWardH2O <- ThisWardH2O_0[StdIndex,]
	ThisLong <- mean(ThisWardH2O_0[401:601,2])
	
	This_Adj_For_Base <- (ThisWardH2O[,2]-ThisLong) #The absorbance of recent Ward H2O, adjusted so that baseline absorbance (abs between 700 and 900 nm) is 0
	Pre_Adj_For_Base <- (PreWardH2O[,2]-PreLong) #Adjust pre-aquashade absorbance so that its 700 to 900 nm absorbance is 0
	This_Adj_For_PreAqua <- This_Adj_For_Base - Pre_Adj_For_Base #The absorbance of recent Ward H2O in the Aquashade region (StdIndex), normalized to the absorbance of Ward H2O before the Aquashade addition
	This_Frac_Of_2ppm <- mean(This_Adj_For_PreAqua / MQ2ppmStandard[,2])
	EstAquashadeConc[i] <- This_Frac_Of_2ppm * 2
	# EstAquashadeConc[i] <- mean(((ThisWardH2O[,2]-PreWardH2O[,2] +PreLong-ThisLong)*2)/(MQ2ppmStandard[,2]))
	# SDAquashadeConc[i] <- sd(((ThisWardH2O[,2]-PreWardH2O[,2] +PreLong-ThisLong)*2)/(MQ2ppmStandard[,2]))
}

# Ward2012Light_DoY <- as.POSIXct(c("2012-04-14 CDT", "2012-04-19 CDT", "2012-04-26 CDT", "2012-05-03 CDT", "2012-05-10 CDT", "2012-05-17 CDT", "2012-05-24 CDT", "2012-05-31 CDT", "2012-06-07 CDT", "2012-06-11 CDT", "2012-06-18 CDT")) #c(105, 110, 117, 124, 131, 138, 145, 152, 159, 163, 170)
# Ward2012_1percLight <-  c(2.3, 2.8, 1.85, 1.75, 1.9, 1.9, 1.8, 1.75, 2.0, 1.5, 2.0)
# ZmixDays <- as.POSIXct(c("2012-04-14 CDT", "2012-04-19 CDT", "2012-04-26 CDT", "2012-05-03 CDT", "2012-05-10 CDT", "2012-05-17 CDT", "2012-05-24 CDT", "2012-05-31 CDT", "2012-06-07 CDT", "2012-06-11 CDT", "2012-06-18 CDT"))  #c(105, 110, 117, 124, 131, 138, 145, 152, 159, 163, 170)
# ZmixDepths <- c(1.0,2.0,2.0,0.5,0.5, 1.0, 1.0, 1.5, 0.5, 1.0, 1.0)

LM_z <- c(2,2,2,1,1, 1)
LM_ppm <- c((0.8203199-0), (1.4287381-0.8203199), (1.8709102-1.1050495), (approx(x=c(1,7), y=c(1.5072962, 1.4256601), xout=2)$y-1.5072962), (1.7712213-1.4158847), (1.8219112-1.2064273))
LM_gal <- c(9,5,10,1,4, 6)
ppm_perGal_perZ <- summary(lm(LM_ppm ~ I(LM_gal/LM_z) -1))$coef[1]
Gal2Add <- (2 - EstAquashadeConc[i]) / ppm_perGal_perZ


plot(OrderedFileDates, EstAquashadeConc, type="o", ylab="Aquahsade (ppm)", xlab="Date", bty="l") #, ylim=c(0,max(EstAquashadeConc))
grid()
# abline(v=AquashadeAdditions[,1], col="blue")
legend("bottom", paste(format(OrderedFileDates[i], "%d-%b-%Y")," = ", round(EstAquashadeConc[i],2), " ppm", "\n", "Add (", as.character(round(Gal2Add,2)), " * Zmix) ", "gallons of Aquahshade", sep=""), bty="n", cex=1.25)

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results")
save(OrderedFileDates, EstAquashadeConc, file="Ward2012_AquaConc.RData")


