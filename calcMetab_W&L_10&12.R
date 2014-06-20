
library(plyr)

# detach(package:LakeMetabolizer, unload=TRUE)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2012.RData")

source("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/lightFunctions.R")

# =============
# = Ward 2010 =
# =============
# Ward 2010 Epilimnion
# Kalman
ward10.epi.kal <- metab(ward10.epi, "kalman")


# Ward 2010 Metalimnion
# Kalman
ward10.meta.kal <- metab(ward10.meta, "kalman")




# =============
# = Ward 2012 =
# =============
# Ward 2012 Epi Metabolism
# Kalman
ward12.epi.kal <- metab(ward12.epi, "kalman")


# Ward 2012 Meta Metabolism
# Kalman
ward12.meta.kal <- metab(ward12.meta, "kalman")




# =============
# = Paul 2010 =
# =============
# Paul 2010 Epilimnion
# Kalman
paul10.epi.kal <- metab(paul10.epi, "kalman")



# =============
# = Paul 2012 =
# =============
# Paul 2012 Epilimnion
# Kalman
paul12.epi.kal <- metab(paul12.epi, "kalman")



# ===========
# = Summary =
# ===========
colMeans(ward10.epi.kal[ward10.epi.kal[,"GPP"]>0&ward10.epi.kal[,"R"]<0,])
colMeans(ward12.epi.kal[ward12.epi.kal[,"GPP"]>0&ward12.epi.kal[,"R"]<0&ward12.epi.kal[,"doy"]>137,])

colMeans(ward10.meta.kal[ward10.meta.kal[,"GPP"]>0&ward10.meta.kal[,"R"]<0,])
colMeans(ward12.meta.kal[ward12.meta.kal[,"GPP"]>0&ward12.meta.kal[,"R"]<0&ward12.meta.kal[,"doy"]>137,])

colMeans(paul10.epi.kal[paul10.epi.kal[,"GPP"]>0&paul10.epi.kal[,"R"]<0&paul10.epi.kal[,"doy"]<241,])
colMeans(paul12.epi.kal[paul12.epi.kal[,"GPP"]>0&paul12.epi.kal[,"R"]<0&paul12.epi.kal[,"doy"]>135,])

kf.epi.good <- rbind(
	cbind("lake"="Ward", ward10.epi.kal[ward10.epi.kal[,"GPP"]>0&ward10.epi.kal[,"R"]<0,]),
	cbind("lake"="Ward", ward12.epi.kal[ward12.epi.kal[,"GPP"]>0&ward12.epi.kal[,"R"]<0&ward12.epi.kal[,"doy"]>137,]),
	cbind("lake"="Paul", paul10.epi.kal[paul10.epi.kal[,"GPP"]>0&paul10.epi.kal[,"R"]<0&paul10.epi.kal[,"doy"]<241,]),
	cbind("lake"="Paul", paul12.epi.kal[paul12.epi.kal[,"GPP"]>0&paul12.epi.kal[,"R"]<0&paul12.epi.kal[,"doy"]>135,])
	)
kf.epi.good[,"datetime"] <- as.POSIXct(do.call(paste, as.list(kf.epi.good[,c("year","doy")])), format="%Y %j", tz="GMT")
rownames(kf.epi.good) <- NULL

# ===================
# = Calculate light =
# ===================

# Read in light profiles
datDir <- "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/"
light.prof00 <- read.table(paste(datDir,"PaulWard_Weekly_2010&2012/PaulWard_Light_2010&2012.csv", sep=""), sep=",", header=TRUE)
light.prof00[,"datetime"] <- as.POSIXct(light.prof00[,"datetime"], tz="GMT")
light.prof00 <- light.prof00[complete.cases(light.prof00),]
light.prof0  <- reshape(light.prof00, v.names=c("fracLight"), timevar="depth", idvar=c("lake","datetime"), direction="wide")
names(light.prof0) <- gsub("^fracLight\\.", "irr_", names(light.prof0))


# Calculate the average zmix for each day
w10.light000 <- cbind("lake"="Ward", LakeMetabolizer:::addNAs(ward10.epi.full)[,c("datetime","doy","z.mix", "irr")])
w12.light000 <- cbind("lake"="Ward", LakeMetabolizer:::addNAs(ward12.epi.full)[,c("datetime","doy","z.mix", "irr")])
l10.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul10.epi.full)[,c("datetime","doy","z.mix", "irr")])
l12.light000 <- cbind("lake"="Paul", LakeMetabolizer:::addNAs(paul12.epi.full)[,c("datetime","doy","z.mix", "irr")])
light000 <- rbind(w10.light000, w12.light000, l10.light000, l12.light000)


light000[,"doy"] <- trunc(light000[,"doy"])
light000[,"datetime"] <- as.POSIXct(trunc.POSIXt(light000[,"datetime"], "days"), tz="GMT")

light00 <- ddply(.data=light000, .variables=c("lake","datetime","doy"), .fun=function(x)colMeans(x[,c("z.mix","irr")], na.rm=TRUE))
light00 <- light00[complete.cases(light00),]


light0 <- merge(light00, light.prof0, all=TRUE)



light <- doLight(light0)


metab.light <- merge(kf.epi.good, light, all.x=TRUE)

gpp.ll <- function(x, layer){
	# E.g., GPP in bot will be GPP.bot <- irr.bot*(GPP/irr.sonde)*size.bot
	# GPP.top <- irr.top*(GPP/irr.sonde)*size.top
	# GPP.sonde <- irr.sonde*(GPP*irr.sonde)*size.sonde == GPP*size.sonde
	# GPP.total <- GPP.top + GPP.bot + GPP.sonde
	# Note that below, when the sonde is in the epi, I am setting size.top = 0, as GPP.top is redundant with GPP.sonde
	
	irr.n <- paste("irr", layer, sep=".")
	dz.n <- paste("dz", layer, sep=".")
	
	gpp.m3 <- x[,"GPP"]*1000 # convert GPP from mg O2 per L per day to mg O2 per m^3 per day
	
	gpp.l <- x[,irr.n]*gpp.m3/x[,"irr.sonde"]*x[,dz.n] # mg O2 m^-2 d^-1; areal units integrate volumetric over "layer"
	
}

gpp.layer <- data.frame(metab.light[,c("lake","datetime","GPP")], gpp.top=)

# =============================================
# = Notes on calculating GPP variance from KF =
# =============================================
# This was taken from SquealMetabolism_v0.3.1.r
# Starting line 825
# Inv_UprimeU <- solve((t(as.matrix(cbind(x[,2], x[,5]), nol=2))%*%as.matrix(cbind(x[,2], x[,5]), nol=2)))
# ParamCovMat <- Inv_UprimeU*DO_KFnll$par[3] #multiplying by par[3] is multiplying by Queue, which "scales" the covariance matrix by the variance of the residuals.
# GPPcoef_Var <- ParamCovMat[1,1] #solve(DO_KFnll$hessian)[1,1]#*DO_KFnll$par[3] #I'm not sure if solve(hessian) gives me a scaled or unscaled matrix.  If it is unscaled and I need scaled, I think I can get that with solve(hessian)*Q
# Rcoef_Var <- ParamCovMat[2,2] #solve(DO_KFnll$hessian)[2,2]#*DO_KFnll$par[3]
# #Qcoef_Var <- solve(DO_KFnll$hessian)[3,3]#*DO_KFnll$par[3]
# GPP_Var <- sum(GPPcoef_Var*x[,2]^2)
# R_Var <- sum(Rcoef_Var*log(x[,5])^2)
# GPP_R_coefs_Cov <- ParamCovMat[1,2] #solve(DO_KFnll$hessian)[1,2]#*DO_KFnll$par[3]
# #f = aA +/- bB
# #var(f) = a^2*var(A) + b^2*var(B) +/- 2ab*cov(AB)
# #var(NEP) = PAR^2*GPPcoef_Var + log(Temp)^2*Rcoef_Var + 2*PAR*log(Temp)*GPP_R_coefs_Cov
# NEP_Var <- sum(GPPcoef_Var*x[,2]^2 + Rcoef_Var*log(x[,5])^2 + 2*x[,2]*x[,5]*abs(GPP_R_coefs_Cov))



