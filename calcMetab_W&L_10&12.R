
detach(package:LakeMetabolizer, unload=TRUE)
install.packages("/Users/Battrd/Documents/School&Work/WiscResearch/LakeMetabolizer", type="source", repos=NULL)
library("LakeMetabolizer")

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2012.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/sondes_ward2010.RData")

# =============
# = Ward 2010 =
# =============

# Bookkeeping
ward10.epi.bk <- metab(ward10.epi, "bookkeep", lake.lat=46.28)

# OLS
ward10.epi.ols <- metab(ward10.epi, "ols")

# MLE
ward10.epi.mle <- metab(ward10.epi, "mle")

# Kalman
ward10.epi.kal <- metab(ward10.epi, "kalman")

# Bayesian
ward10.epi.bay <- metab(ward10.epi, "bayesian")

# merge all results
w10e.r1 <- merge(cbind("method"="bookkeep",ward10.epi.bk), cbind("method"="ols",ward10.epi.ols), all=TRUE)
w10e.r2 <- merge(w10e.r1, cbind("method"="mle", ward10.epi.mle), all=TRUE)
w10e.r3 <- merge(w10e.r2, cbind("method"="kalman", ward10.epi.kal), all=TRUE)
w10e.r4 <- merge(w10e.r3, cbind("method"="bayes", ward10.epi.bay), all=TRUE)


	# =========================
	# = Ward 2010 Metalimnion =
	# =========================
# Bookkeeping
ward10.meta.bk <- metab(ward10.meta, "bookkeep", lake.lat=46.28)

# OLS
ward10.meta.ols <- metab(ward10.meta, "ols")

# MLE
ward10.meta.mle <- metab(ward10.meta, "mle")

# Kalman
ward10.meta.kal <- metab(ward10.meta, "kalman")

# Bayesian
ward10.meta.bay <- metab(ward10.meta, "bayesian")

# merge Ward 2010 Metalimnion results
w10m.r1 <- merge(cbind("method"="bookkeep",ward10.meta.bk), cbind("method"="ols",ward10.meta.ols), all=TRUE)
w10m.r2 <- merge(w10m.r1, cbind("method"="mle", ward10.meta.mle), all=TRUE)
w10m.r3 <- merge(w10m.r2, cbind("method"="kalman", ward10.meta.kal), all=TRUE)
w10m.r4 <- merge(w10m.r3, cbind("method"="bayes", ward10.meta.bay), all=TRUE)




# =============
# = Ward 2012 =
# =============

	# ============================
	# = Ward 2012 Epi Metabolism =
	# ============================

# Bookkeeping
ward12.epi.bk <- metab(ward12.epi, "bookkeep", lake.lat=46.28)

# OLS
ward12.epi.ols <- metab(ward12.epi, "ols")

# MLE
ward12.epi.mle <- metab(ward12.epi, "mle")

# Kalman
ward12.epi.kal <- metab(ward12.epi, "kalman")

# Bayesian
ward12.epi.bay <- metab(ward12.epi, "bayesian")

# merge Ward 2012 Epilimnion results
w12e.r1 <- merge(cbind("method"="bookkeep",ward12.epi.bk), cbind("method"="ols",ward12.epi.ols), all=TRUE)
w12e.r2 <- merge(w12e.r1, cbind("method"="mle", ward12.epi.mle), all=TRUE)
w12e.r3 <- merge(w12e.r2, cbind("method"="kalman", ward12.epi.kal), all=TRUE)
w12e.r4 <- merge(w12e.r3, cbind("method"="bayes", ward12.epi.bay), all=TRUE)


	# ============================
	# = Ward 2012 Meta Metabolism =
	# ============================

# Bookkeeping
ward12.meta.bk <- metab(ward12.meta, "bookkeep", lake.lat=46.28)

# OLS
ward12.meta.ols <- metab(ward12.meta, "ols")

# MLE
ward12.meta.mle <- metab(ward12.meta, "mle")

# Kalman
ward12.meta.kal <- metab(ward12.meta, "kalman")

# Bayesian
ward12.meta.bay <- metab(ward12.meta, "bayesian")

# merge Ward 2012 Metalimnion results
w12m.r1 <- merge(cbind("method"="bookkeep",ward12.meta.bk), cbind("method"="ols",ward12.meta.ols), all=TRUE)
w12m.r2 <- merge(w12m.r1, cbind("method"="mle", ward12.meta.mle), all=TRUE)
w12m.r3 <- merge(w12m.r2, cbind("method"="kalman", ward12.meta.kal), all=TRUE)
w12m.r4 <- merge(w12m.r3, cbind("method"="bayes", ward12.meta.bay), all=TRUE)








# save.image("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/ward_2010&2012_metabolism_allModels_LakeMetabolizer.RData")



