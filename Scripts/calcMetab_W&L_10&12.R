
library(plyr)

# detach(package:LakeMetabolizer, unload=TRUE)
# install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2010.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_paul2012.RData")

chooseGas <- "kgas.read"

# =============
# = Ward 2010 =
# =============
# Ward 2010 Epilimnion
# Kalman
ward10.epi[,"k.gas"] <- ward10.epi[, chooseGas]
ward10.epi.kal <- metab(ward10.epi, "kalman")
ward10.epi.bk <- metab(ward10.epi, "bookkeep", lake.lat=46.28)


# Ward 2010 Metalimnion
# Kalman
ward10.meta.kal <- metab(ward10.meta, "kalman")




# =============
# = Ward 2012 =
# =============
# Ward 2012 Epi Metabolism
# Kalman
ward12.epi[,"k.gas"] <- ward12.epi[,chooseGas]
ward12.epi.kal <- metab(ward12.epi, "kalman")
ward12.epi.bk <- metab(ward12.epi, "bookkeep", lake.lat=46.28)


# Ward 2012 Meta Metabolism
# Kalman
ward12.meta.kal <- metab(ward12.meta, "kalman")




# =============
# = Paul 2010 =
# =============
# Paul 2010 Epilimnion
# Kalman
paul10.epi[,"k.gas"] <- paul10.epi[,chooseGas]
paul10.epi.kal <- metab(paul10.epi, "kalman")



# =============
# = Paul 2012 =
# =============
# Paul 2012 Epilimnion
# Kalman
paul12.epi[,"k.gas"] <- paul12.epi[,chooseGas]
paul12.epi.kal <- metab(paul12.epi, "kalman")




# ==================
# = Summary Kalman =
# ==================
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

WardPaulMetabolism.kf <- rbind(
	cbind("lake"="Ward", ward10.epi.kal),
	cbind("lake"="Ward", ward12.epi.kal),
	cbind("lake"="Paul", paul10.epi.kal),
	cbind("lake"="Paul", paul12.epi.kal)
	)
	
# =======================
# = Summary Bookkeeping =
# =======================
colMeans(ward10.epi.bk[ward10.epi.bk[,"GPP"]>0&ward10.epi.bk[,"R"]<0,])
colMeans(ward12.epi.bk[ward12.epi.bk[,"GPP"]>0&ward12.epi.bk[,"R"]<0&ward12.epi.bk[,"doy"]>137,])


# ================
# = Save results =
# ================

save(ward10.epi.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward10.epi.kal.RData")
save(ward12.epi.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward12.epi.kal.RData")
save(ward10.meta.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward10.meta.kal.RData")
save(ward12.meta.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/ward12.meta.kal.RData")
save(paul10.epi.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/paul10.epi.kal.RData")
save(paul12.epi.kal, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/paul12.epi.kal.RData")

save(kf.epi.good, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/kf.epi.good.RData")
save(WardPaulMetabolism.kf, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/WardPaulMetabolism.kf.RData")


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



