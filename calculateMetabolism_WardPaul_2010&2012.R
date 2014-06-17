# =============
# = Ward 2010 =
# =============

# Bookkeeping
ward10.epi.bk <- metab(ward10.epi, "bookkeep", lake.lat=46.28)
# ward10.epi.bk.res <- ward10.epi.bk[,c("doy","GPP","R", "NEP")]
# ward10.epi.bk.res <- merge(data.frame("doy"=138:239), ward10.epi.bk.res, all=TRUE)

# OLS
ward10.epi.ols <- metab(ward10.epi, "ols")
# ward10.epi.ols.res <- ward10.epi.ols[,c("doy","GPP","R", "NEP")]
# ward10.epi.ols.res <- merge(data.frame("doy"=138:239), ward10.epi.ols.res, all=TRUE)

# MLE
ward10.epi.mle <- metab(ward10.epi, "mle")
# ward10.epi.mle.res <- ward10.epi.mle[,c("doy","GPP","R", "NEP")]
# ward10.epi.mle.res <- merge(data.frame("doy"=138:239), ward10.epi.mle.res, all=TRUE)

# Kalman
ward10.epi.kal <- metab(ward10.epi, "kalman")
# ward10.epi.kal.res <- ward10.epi.kal[,c("doy","GPP","R", "NEP")]
# ward10.epi.kal.res <- merge(data.frame("doy"=138:239), ward10.epi.kal.res, all=TRUE)

# Bayesian
ward10.epi.bay <- metab(ward10.epi, "bayesian")
# ward10.epi.bay.res <- ward10.epi.bay[,c("doy","GPP","R", "NEP")]
# ward10.epi.bay.res <- merge(data.frame("doy"=138:239), ward10.epi.bay.res, all=TRUE)

# merge all results
w10e.r1 <- merge(cbind("method"="bookkeep",ward12.epi.bk), cbind("method"="ols",ward12.epi.ols), all=TRUE)
w10e.r2 <- merge(w10e.r1, cbind("method"="mle", ward12.epi.mle), all=TRUE)
w10e.r3 <- merge(w10e.r2, cbind("method"="kalman", ward12.epi.kal), all=TRUE)
w10e.r4 <- merge(w10e.r3, cbind("method"="bayes", ward12.epi.bay), all=TRUE)


	# =========================
	# = Ward 2010 Metalimnion =
	# =========================
# Bookkeeping
ward10.meta.bk <- metab(ward10.meta, "bookkeep", lake.lat=46.28)
# ward10.meta.bk.res <- ward10.meta.bk[,c("doy","GPP","R", "NEP")]
# ward10.meta.bk.res <- merge(data.frame("doy"=138:239), ward10.meta.bk.res, all=TRUE)

# OLS
ward10.meta.ols <- metab(ward10.meta, "ols")
# ward10.meta.ols.res <- ward10.meta.ols[,c("doy","GPP","R", "NEP")]
# ward10.meta.ols.res <- merge(data.frame("doy"=138:239), ward10.meta.ols.res, all=TRUE)

# MLE
ward10.meta.mle <- metab(ward10.meta, "mle")
# ward10.meta.mle.res <- ward10.meta.mle[,c("doy","GPP","R", "NEP")]
# ward10.meta.mle.res <- merge(data.frame("doy"=138:239), ward10.meta.mle.res, all=TRUE)

# Kalman
ward10.meta.kal <- metab(ward10.meta, "kalman")
# ward10.meta.kal.res <- ward10.meta.kal[,c("doy","GPP","R", "NEP")]
# ward10.meta.kal.res <- merge(data.frame("doy"=138:239), ward10.meta.kal.res, all=TRUE)

# Bayesian
ward10.meta.bay <- metab(ward10.meta, "bayesian")
# ward10.meta.bay.res <- ward10.meta.bay[,c("doy","GPP","R", "NEP")]
# ward10.meta.bay.res <- merge(data.frame("doy"=138:239), ward10.meta.bay.res, all=TRUE)

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
# ward12.epi.bk <- merge(data.frame("doy"=do.call(":", as.list(range(ward12.epi.bk[,"doy"])))), ward12.epi.bk, all=TRUE)

# OLS
ward12.epi.ols <- metab(ward12.epi, "ols")
# ward12.epi.ols.res <- ward12.epi.ols[,c("doy","GPP","R", "NEP")]
# ward12.epi.ols.res <- merge(data.frame("doy"=138:239), ward12.epi.ols.res, all=TRUE)

# MLE
ward12.epi.mle <- metab(ward12.epi, "mle")
# ward12.epi.mle.res <- ward12.epi.mle[,c("doy","GPP","R", "NEP")]
# ward12.epi.mle.res <- merge(data.frame("doy"=138:239), ward12.epi.mle.res, all=TRUE)

# Kalman
ward12.epi.kal <- metab(ward12.epi, "kalman")
# ward12.epi.kal.res <- ward12.epi.kal[,c("doy","GPP","R", "NEP")]
# ward12.epi.kal.res <- merge(data.frame("doy"=138:239), ward12.epi.kal.res, all=TRUE)

# Bayesian
ward12.epi.bay <- metab(ward12.epi, "bayesian")
# ward12.epi.bay.res <- ward12.epi.bay[,c("doy","GPP","R", "NEP")]
# ward12.epi.bay.res <- merge(data.frame("doy"=138:239), ward12.epi.bay.res, all=TRUE)

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
# ward12.meta.bk.res <- ward12.meta.bk[,c("doy","GPP","R", "NEP")]
# ward12.meta.bk.res <- merge(data.frame("doy"=138:239), ward12.meta.bk.res, all=TRUE)

# OLS
ward12.meta.ols <- metab(ward12.meta, "ols")
# ward12.meta.ols.res <- ward12.meta.ols[,c("doy","GPP","R", "NEP")]
# ward12.meta.ols.res <- merge(data.frame("doy"=138:239), ward12.meta.ols.res, all=TRUE)

# MLE
ward12.meta.mle <- metab(ward12.meta, "mle")
# ward12.meta.mle.res <- ward12.meta.mle[,c("doy","GPP","R", "NEP")]
# ward12.meta.mle.res <- merge(data.frame("doy"=138:239), ward12.meta.mle.res, all=TRUE)

# Kalman
ward12.meta.kal <- metab(ward12.meta, "kalman")
# ward12.meta.kal.res <- ward12.meta.kal[,c("doy","GPP","R", "NEP")]
# ward12.meta.kal.res <- merge(data.frame("doy"=138:239), ward12.meta.kal.res, all=TRUE)

# Bayesian
ward12.meta.bay <- metab(ward12.meta, "bayesian")
# ward12.meta.bay.res <- ward12.meta.bay[,c("doy","GPP","R", "NEP")]
# ward12.meta.bay.res <- merge(data.frame("doy"=138:239), ward12.meta.bay.res, all=TRUE)

# merge Ward 2012 Metalimnion results
w12m.r1 <- merge(cbind("method"="bookkeep",ward12.meta.bk), cbind("method"="ols",ward12.meta.ols), all=TRUE)
w12m.r2 <- merge(w12m.r1, cbind("method"="mle", ward12.meta.mle), all=TRUE)
w12m.r3 <- merge(w12m.r2, cbind("method"="kalman", ward12.meta.kal), all=TRUE)
w12m.r4 <- merge(w12m.r3, cbind("method"="bayes", ward12.meta.bay), all=TRUE)








# ==========
# = lower? =
# ==========
good10 <- ward10.epi.kal[,"GPP"] > 0 & ward10.epi.kal[,"R"] < 0
colMeans(ward10.epi.kal[good10,])

good12 <- ward12.epi.kal[,"GPP"] > 0 & ward12.epi.kal[,"R"] < 0
late12 <-  (ward12.epi.kal[,"doy"] >= min(ward10.epi.kal[,"doy"]))
early12 <-  (ward12.epi.kal[,"doy"] < min(ward10.epi.kal[,"doy"]))
colMeans(ward12.epi.kal[good12&early12,])
colMeans(ward12.epi.kal[good12&late12,])


wa12.kal.g2r <- ward12.epi.kal[good12,"GPP"]/-ward12.epi.kal[good12, "R"]



