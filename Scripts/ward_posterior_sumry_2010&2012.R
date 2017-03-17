load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")
library(plyr)


# ===============
# = Print Means =
# ===============
for(i in 1:length(Cons)){
	NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	myRU <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	myRU <- reshape(myRU, varying=list(c(NamesThisUse)), times=ConsChoices[[Cons[i]]][[GroupChoose]], ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")
	t.means <- ddply(myRU, c("Year", "Consumer", "Grouping", "Source"), function(x)data.frame("Proportion"=mean(x[,"Proportion"])))
	if(i==1){all.means <- t.means}else{all.means <- rbind(all.means, t.means)}
	print(t.means)
}

# ==============
# = Print SD's =
# ==============
for(i in 1:length(Cons)){
	NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	myRU <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	myRU <- reshape(myRU, varying=list(c(NamesThisUse)), times=ConsChoices[[Cons[i]]][[GroupChoose]], ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")
	print(ddply(myRU, c("Year", "Consumer", "Grouping", "Source"), function(x)data.frame("Proportion"=sd(x[,"Proportion"]))))
}


# ================================
# = Posterior of POM Composition =
# ================================
pomPost2010 <- c(apply(bugsOut_pt1_2010$BUGSoutput$sims.matrix, 2, mean)[c("f[3]","f[4]")], apply(bugsOut_pt1_2010$BUGSoutput$sims.matrix, 2, sd)[c("f[3]","f[4]")])
names(pomPost2010) <- c("terr_mu", "phyto_mu", "terr_sd", "phyto_sd")
pomPost2010

pomPost2012 <- c(apply(bugsOut_pt1_2012$BUGSoutput$sims.matrix, 2, mean)[c("f[3]","f[4]")], apply(bugsOut_pt1_2012$BUGSoutput$sims.matrix, 2, sd)[c("f[3]","f[4]")])
names(pomPost2012) <- c("terr_mu", "phyto_mu", "terr_sd", "phyto_sd")
pomPost2012

# ================================
# = Posterior of DOM composition =
# ================================
load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/DOM_Comp_Ward2010&2012.RData")

domPost2010 <- c(apply(domOut2010$BUGSoutput$sims.matrix, 2, mean)[1:4], apply(domOut2010$BUGSoutput$sims.matrix, 2, sd)[1:4])
names(domPost2010) <- c("Terr_mu", "Macroph_mu", "Phyto_mu", "Periphy_mu", "Terr_sd", "Macroph_sd", "Phyto_sd", "Periphy_sd")
domPost2010

domPost2012 <- c(apply(domOut2012$BUGSoutput$sims.matrix, 2, mean)[1:4], apply(domOut2012$BUGSoutput$sims.matrix, 2, sd)[1:4])
names(domPost2012) <- c("Terr_mu", "Macroph_mu", "Phyto_mu", "Periphy_mu", "Terr_sd", "Macroph_sd", "Phyto_sd", "Periphy_sd")
domPost2012