#_v0.1 (13-Feb-2013) Looking for the effect of zmix on metabolism
rm(list=ls())
graphics.off()
library("plyr")


# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
# load("WardMetabolism_WardMetab2012_v1.RData")
# load("Chla_WardMetab2012_v1.RData")
# load("Photic_WardMetab2012_v1.RData")
# load("DOM_WardMetab2012_v1.RData")
# load("Ward2012_AquaConc.RData")


setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data")
load("Paul_Metab&Data_2010&2012.RData")
load("Ward_Metab&Data_2010&2012.RData")

WardMetabolism2 <- data.frame("Lake"="Ward", subset(WardMetabolism, Method=="BK")[,])
PaulMetabolism2 <- data.frame("Lake"="Paul", subset(PaulMetabolism, Method=="BK")[,])
WardPaul_Metabolism <- subset(rbind(WardMetabolism2, PaulMetabolism2), Week>0 & Week<5)

# AllZ <- rbind(data.frame("Lake"="Paul", Paul_Tchain_Weather_2010[,c("Year","DoY","Zmix")]), data.frame("Lake"="Paul", Paul_Tchain_Weather_2012[,c("Year","DoY","Zmix")]), data.frame("Lake"="Ward", Ward_Data_2010[,c("Year","DoY","Zmix")]), data.frame("Lake"="Ward", Ward_Tchain_Weather_2012[,c("Year","DoY","Zmix")]))
# AllZ[,"DoY"] <- trunc(AllZ[,"DoY"])
# AllZdaily <- aggregate(AllZ[,"Zmix"], by=list(AllZ[,"DoY"], AllZ[,"Lake"], AllZ[,"Year"]), FUN=mean)
# names(AllZdaily) <- c("DoY", "Lake", "Year", "Zmix")
# AllZdaily <- AllZdaily[complete.cases(AllZdaily),]



# AllDO <- rbind(data.frame("Lake"="Paul", PaulData2010[,c("Year","DoY","DOsat")]), data.frame("Lake"="Paul", PaulData2010[,c("Year","DoY","DOsat")]), data.frame("Lake"="Ward", Ward_Data_2010[,c("Year","DoY","DOsat")]), data.frame("Lake"="Ward", Ward_Data_2012[,c("Year","DoY","DOsat")]))
# AllDO[,"DoY"] <- trunc(AllDO[,"DoY"])
# AllDOdaily <- aggregate(AllDO[,"DOsat"], by=list(AllDO[,"DoY"], AllDO[,"Lake"], AllDO[,"Year"]), FUN=mean)
# names(AllDOdaily) <- c("DoY", "Lake", "Year", "DOsat")

# merge(WardPaul_Metabolism, AllDOdaily)
# ====================
# = Metabolism Plots =
# ====================

Metabolism_Plots <- expand.grid(c("Paul", "Ward"), c("GPP", "R", "NEP"))
dev.new(width=6, height=8)
par(mfrow=c(3,2), mar=c(1,3,1,1), ps=10, cex=1, oma=c(1,1,2,0))
# GPPYlim <- range(WardPaul_Metabolism[,"GPP"]
# RYlim <- range(WardPaul_Metabolism[,"R"]
# NEPYlim <- range(WardPaul_Metabolism[,"NEP"]
Ylims <- data.frame("GPP"=range(WardPaul_Metabolism[,"GPP"]), "R"=range(WardPaul_Metabolism[,"R"]), "NEP"=range(WardPaul_Metabolism[,"NEP"]))

for(i in 1:nrow(Metabolism_Plots)){
	MetabD <- subset(WardPaul_Metabolism, Lake==as.character(Metabolism_Plots[i,1]))
	# Ylim <- Ylims[,as.character(Metabolism_Plots[i,2])]
	Xaxt <- ifelse(is.element(i, c(5,6)), "s", "n")
	MetabMonths <- c("Apr","May","Jun","Jul","Aug", "Sept")[1+expand.grid(sort(unique(MetabD[,"Week"])), sort(unique(MetabD[,"Year"])))[,1]]
	RepYearCol <- length(MetabMonths)/2
	
	Response <- MetabD[,as.character(Metabolism_Plots[i,2])]
	Pred_Week <- as.numeric(MetabD[,"Week"])
	Pred_Year <- as.factor(MetabD[,"Year"])
	# dev.new()
	boxplot(Response~Pred_Week*Pred_Year, data=MetabD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=MetabMonths, xaxt=Xaxt)
	if(Xaxt=="n"){
		axis(side=1, labels=FALSE)
	}
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	# mtext(as.character(Metabolism_Plots[i,]), side=2, line=2.5)
	if(i==1){mtext("Paul", side=3, line=1)}
	if(i==1){mtext("GPP", side=2, line=2.5)}
	if(i==2){mtext("Ward", side=3, line=1)}
	if(i==3){mtext("R", side=2, line=2.5)}
	if(i==5){mtext("NEP", side=2, line=2.5)}
	if(i==2){legend("topleft", c("2010", "2012"), text.col=c("red", "blue"), bty="n")}
}





Var2Merge <- c("Year", "Week", "DoY", "Temp", "DOsat")
AllPaulSonde <- data.frame("Lake"="Paul", rbind(PaulData2010[,Var2Merge], PaulData2012[,Var2Merge]))


WardWeek2010 <- as.numeric(format.Date(as.POSIXct(paste(2010, Ward_Data_2010[,"DoY"], sep="-"), format="%Y-%j"), format="%m"))-4
WardWeek2012 <- as.numeric(format.Date(as.POSIXct(paste(2012, Ward_Data_2012[,"DoY"], sep="-"), format="%Y-%j"), format="%m"))-4
Ward_Data_2010[,"Week"] <- WardWeek2010
Ward_Data_2012[,"Week"] <- WardWeek2012
AllWardSonde <- data.frame("Lake"="Ward", rbind(Ward_Data_2010[,Var2Merge], Ward_Data_2012[,Var2Merge]))

AllSonde <- rbind(AllPaulSonde, AllWardSonde)
AllSonde[,"RoundDoY"] <- trunc(AllSonde[,"DoY"])

TempStat <- function(x){
	xtemp <- x[,"Temp"]
	xstat <- data.frame("Lake"=unique(x[,"Lake"]), "Year"=unique(x[,"Year"]), "Week"=unique(x[,"Week"]), "DoY"=unique(x[,"RoundDoY"]), "MeanTemp"=mean(xtemp), "MinTemp"=min(xtemp), "MaxTemp"=max(xtemp), "RangeTemp"=diff(range(xtemp)))
	return(xstat)
}
AllTemp <- ddply(AllSonde, .variables=c("Lake", "Year","RoundDoY"), .fun=TempStat)


AllTemp <- subset(AllTemp, Week>0 & Week<5)[,]
Temp_Plots <- expand.grid(c("Paul", "Ward"), c("MeanTemp", "MinTemp", "MaxTemp", "RangeTemp"))
dev.new(width=6, height=9)
par(mfrow=c(4,2), mar=c(1,1,1,1), ps=10, cex=1, oma=c(1,3,2,0))
# GPPYlim <- range(WardPaul_Metabolism[,"GPP"]
# RYlim <- range(WardPaul_Metabolism[,"R"]
# NEPYlim <- range(WardPaul_Metabolism[,"NEP"]
Ylims <- data.frame("MeanTemp"=range(AllTemp[,"MeanTemp"]), "MinTemp"=range(AllTemp[,"MinTemp"]), "MaxTemp"=range(AllTemp[,"MaxTemp"]), "RangeTemp"=range(AllTemp[,"RangeTemp"]))

for(i in 1:nrow(Temp_Plots)){
	TempD <- subset(AllTemp, Lake==as.character(Temp_Plots[i,1]))
	Ylim <- Ylims[,as.character(Temp_Plots[i,2])]
	Xaxt <- ifelse(is.element(i, c(7,8)), "s", "n")
	Yaxt <- ifelse(is.element(i, c(1,3,5,7)), "s", "n")
	MetabMonths <- c("Apr","May","Jun","Jul","Aug", "Sept")[1+expand.grid(sort(unique(TempD[,"Week"])), sort(unique(TempD[,"Year"])))[,1]]
	RepYearCol <- length(MetabMonths)/2
	
	Response <- TempD[,as.character(Temp_Plots[i,2])]
	Pred_Week <- as.numeric(TempD[,"Week"])
	Pred_Year <- as.factor(TempD[,"Year"])
	# dev.new()
	boxplot(Response~Pred_Week*Pred_Year, data=TempD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=MetabMonths, xaxt=Xaxt, ylim=Ylim, yaxt=Yaxt)
	if(Xaxt=="n"){
		axis(side=1, labels=FALSE)
	}
	if(Yaxt=="n"){
		axis(side=2, labels=FALSE)
	}
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	# mtext(as.character(Metabolism_Plots[i,]), side=2, line=2.5)
	if(i==1){mtext("Paul", side=3, line=1)}
	if(i==1){mtext("MeanTemp", side=2, line=2.5)}
	if(i==2){mtext("Ward", side=3, line=1)}
	if(i==3){mtext("MinTemp", side=2, line=2.5)}
	if(i==5){mtext("MaxTemp", side=2, line=2.5)}
	if(i==7){mtext("RangeTemp", side=2, line=2.5)}
	if(i==7){legend("topleft", c("2010", "2012"), text.col=c("red", "blue"), bty="n")}
}


dim(AllTemp)
dim(WardPaul_Metabolism)

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/PAR_DO_PhaseShift/")
source("Alme.R")
AllMetabTemp <- merge(AllTemp, WardPaul_Metabolism, by=c("Lake", "Year", "DoY"))


dev.new()
plot(AllMetabTemp[,"MinTemp"], log(-AllMetabTemp[,"R"]))

yR <- log(-AllMetabTemp[,"R"])
xTemp <- AllMetabTemp[,"MeanTemp"]
xLake <- AllMetabTemp[,"Lake"]
xYear <- as.factor(AllMetabTemp[,"Year"])
lmsum <- summary(lm(yR~xTemp*xLake+xYear))
# abline(lm(y1~x1*x2))

yGPP <- AllMetabTemp[,"GPP"]
xTemp <- AllMetabTemp[,"MeanTemp"]
xLake <- AllMetabTemp[,"Lake"]
xYear <- as.factor(AllMetabTemp[,"Year"])
lmsum <- summary(lm(yGPP~xLake*xYear))


Resp <- log(-AllMetabTemp[,"R"])
x1 <- AllMetabTemp[,"MinTemp"]
x2 <- AllMetabTemp[,"Lake"]
x3 <- AllMetabTemp[,"Year"]
summary(lm(MeanTemp~Lake*Year +as.factor(Week.x)*Year, data=AllMetabTemp))










DOStat <- function(x){
	xtemp <- x[,"DOsat"]
	xstat <- data.frame("Lake"=unique(x[,"Lake"]), "Year"=unique(x[,"Year"]), "Week"=unique(x[,"Week"]), "DoY"=unique(x[,"RoundDoY"]), "MeanDO"=mean(xtemp), "MinDO"=min(xtemp), "MaxDO"=max(xtemp), "RangeDO"=diff(range(xtemp)))
	return(xstat)
}
AllDO <- ddply(AllSonde, .variables=c("Lake", "Year","RoundDoY"), .fun=DOStat)


AllDO <- subset(AllDO, Week>0 & Week<5)[,]
Temp_Plots <- expand.grid(c("Paul", "Ward"), c("MeanDO"))
dev.new(width=5, height=3)
par(mfrow=c(1,2), mar=c(1,1,1,1), ps=10, cex=1, oma=c(1,3,2,0))
# GPPYlim <- range(WardPaul_Metabolism[,"GPP"]
# RYlim <- range(WardPaul_Metabolism[,"R"]
# NEPYlim <- range(WardPaul_Metabolism[,"NEP"]
Ylims <- data.frame("MeanDO"=range(AllDO[,"MeanDO"]))

for(i in 1:nrow(Temp_Plots)){
	TempD <- subset(AllDO, Lake==as.character(Temp_Plots[i,1]))
	Ylim <- Ylims[,as.character(Temp_Plots[i,2])]
	Xaxt <- ifelse(is.element(i, c(1,2)), "s", "n")
	Yaxt <- ifelse(is.element(i, c(1)), "s", "n")
	MetabMonths <- c("Apr","May","Jun","Jul","Aug", "Sept")[1+expand.grid(sort(unique(TempD[,"Week"])), sort(unique(TempD[,"Year"])))[,1]]
	RepYearCol <- length(MetabMonths)/2
	
	Response <- TempD[,as.character(Temp_Plots[i,2])]
	Pred_Week <- as.numeric(TempD[,"Week"])
	Pred_Year <- as.factor(TempD[,"Year"])
	# dev.new()
	boxplot(Response~Pred_Week*Pred_Year, data=TempD, col=c(rep("#FA807225",RepYearCol), rep("#3A5FCD25",RepYearCol)), border=c(rep("red",RepYearCol), rep("blue",RepYearCol)), names=MetabMonths, xaxt=Xaxt, ylim=Ylim, yaxt=Yaxt)
	if(Xaxt=="n"){
		axis(side=1, labels=FALSE)
	}
	if(Yaxt=="n"){
		axis(side=2, labels=FALSE)
	}
	# abline(v=8.5)
	# text(c(2, 10.5), -240, c("Epilimnion", "Metalimnion"), font=4, cex=1.25)
	# mtext(as.character(Metabolism_Plots[i,]), side=2, line=2.5)
	if(i==1){mtext("Paul", side=3, line=1)}
	if(i==1){mtext("Mean DO %", side=2, line=2.5)}
	if(i==2){mtext("Ward", side=3, line=1)}
	if(i==1){legend("bottomleft", c("2010", "2012"), text.col=c("red", "blue"), bty="n")}
}

# AllSonde <- subset(AllSonde, Week>0 & Week<5)
# aggregate(AllSonde[,"DOsat"], by=list(AllSonde[,"Lake"],AllSonde[,"Year"]), FUN=mean)
# aggregate(AllSonde[,"Temp"], by=list(AllSonde[,"Lake"],AllSonde[,"Year"]), FUN=mean)

aggregate(WardPaul_Metabolism[,c("GPP","R","NEP")], by=list(WardPaul_Metabolism[,"Lake"],WardPaul_Metabolism[,"Year"]), FUN=max)




png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/NEP_PaulWard_2010&2012_Metabolism_v0.2.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
par(mar=c(2.5,3.5,0.5,0.5), ps=10)
pWardPaul_Metabolism <- subset(WardPaul_Metabolism, DoY>=143)
pWardPaul_Metabolism[,"Lake"] <- relevel(pWardPaul_Metabolism[,"Lake"], ref="Paul")
boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
legend("bottomleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext(expression(Net~Production~(mmol~O[2]~m^-3~d^-1)), side=2, line=2.25)
dev.off()

png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/DOpercent_PaulWard_2010&2012_Metabolism_v0.2.png", units="in", res=300, height=3.5, width=3.5, pointsize=10)
par(mar=c(2.5,3.5,0.5,0.5), ps=10)
boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), show.names=FALSE, outline=FALSE, lwd=1.5)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
legend("bottomleft", legend=c(2010, 2012), text.col=c("red","blue"), bty="n")
mtext("Dissolved Oxygen (% saturtation)", side=2, line=2.25)
dev.off()


