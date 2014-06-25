
library(plyr)

# detach(package:LakeMetabolizer, unload=TRUE)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2012.RData")


# ===========================================
# = Make the AllDO variable for back compat =
# ===========================================
Conc2Sat <- function(SDO2conc, SDtemp){  
	SDO2mgL <- SDO2conc
	# SDO2mgL=SDO2conc/1000*(15.999*2)
	SDO2sat <- (SDO2mgL*100)/(-0.00000002057759*SDtemp^5 + 0.000002672016*SDtemp^4 + -0.0001884085*SDtemp^3 + 0.009778012*SDtemp^2 + -0.4147241*SDtemp + 14.621)
	return(SDO2sat)
}

ward10.epi.full[,"DOsat"] <- Conc2Sat(ward10.epi.full[,"do.obs"],  ward10.epi.full[,"wtr"])
ward12.epi.full[,"DOsat"] <- Conc2Sat(ward12.epi.full[,"do.obs"],  ward12.epi.full[,"wtr"])
paul10.epi.full[,"DOsat"] <- Conc2Sat(paul10.epi.full[,"do.obs"],  paul10.epi.full[,"wtr"])
paul12.epi.full[,"DOsat"] <- Conc2Sat(paul12.epi.full[,"do.obs"],  paul12.epi.full[,"wtr"])

allDO0 <- rbind(
	cbind("lake"="ward", ward10.epi.full[,c("datetime","DOsat")]),
	cbind("lake"="ward", ward12.epi.full[,c("datetime","DOsat")]),
	cbind("lake"="paul", paul10.epi.full[,c("datetime","DOsat")]),
	cbind("lake"="paul", paul12.epi.full[,c("datetime","DOsat")])
	)

allDO0[,"doy"] <- as.integer(format.Date(allDO0[,"datetime"], format="%j"))
allDO0[,"year"] <- as.integer(format.Date(allDO0[,"datetime"], format="%Y"))

allDO <- allDO0[allDO0[,"doy"]>=143 & allDO0[,"doy"]<241,]

AllDO <- ddply(allDO, c("lake","year","doy"), function(x)data.frame("MeanDO"=mean(x[,"DOsat"])))
AllDO[,"lake"] <- relevel(AllDO[,"lake"], ref="paul")

save(AllDO, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/AllDO.RData")

# ==============
# = experiment =
# ==============
allDO.time <- LakeMetabolizer:::date2doy(allDO[,"datetime"])
allDO[,"time"] <- factor(allDO.time - trunc(allDO.time))
plot(allDO[allDO[,"lake"]=="ward"&allDO[,"year"]==2010,"time"],allDO[allDO[,"lake"]=="ward"&allDO[,"year"]==2010,"DOsat"])

allDOsum0 <- ddply(allDO, c("lake","year", "time"), function(x)data.frame("muDO"=mean(x[,"DOsat"], na.rm=TRUE), "sdDO"=sd(x[,"DOsat"], na.rm=TRUE)))


XY_ErrorBars <- function(X, Y, sdX, sdY, Plot=FALSE){#This is a function that creates error bars given means and sd's
	#Vertical Error Bars
	XnaughtV <- X
	XoneV <- X
	YnaughtV <- Y-sdY
	YoneV <- Y+sdY
	if(Plot==TRUE){polygon(x=c(XnaughtV,rev(XoneV)), y=c(YnaughtV,rev(YoneV)), col="gray")}

}#End function

dev.new()
par(mfrow=c(2,2))
do.lims <- c(min(allDOsum0[,"muDO"]-allDOsum0[,"sdDO"]), max(allDOsum0[,"muDO"]+allDOsum0[,"sdDO"]))
for(l in 1:2){
	for(y in 1:2){
		t.dat <- allDOsum0[allDOsum0[,"lake"]==c("paul","ward")[l]&allDOsum0[,"year"]==c(2010,2012)[y],]
		
		t.lims <- XY_ErrorBars(X=as.numeric(as.character(t.dat[,"time"])), Y=t.dat[,"muDO"], sdX=0, sdY=t.dat[,"sdDO"])
		
		plot(as.numeric(as.character(t.dat[,"time"])), t.dat[,"muDO"], type="n", ylim=do.lims)
		
		XY_ErrorBars(X=as.numeric(as.character(t.dat[,"time"])), Y=t.dat[,"muDO"], sdX=0, sdY=t.dat[,"sdDO"], Plot=TRUE)
		
		points(as.numeric(as.character(t.dat[,"time"])), t.dat[,"muDO"], pch=20)
		
	}
}



