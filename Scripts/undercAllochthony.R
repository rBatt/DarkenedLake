
allo00 <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/undercAllochthony.csv", sep=",", header=TRUE)
allo0 <- allo00[!allo00[,"notes"]%in%c("Nut Addition"),]
# allo <- ddply(allo0, c("consumerGroup","consumer","lake","year", "study"), function(x)data.frame(allochthony=mean(x[,"allochthony"])))
allo <- ddply(allo0, c("consumerGroup","consumer","lake","year"), function(x)data.frame(allochthony=mean(x[,"allochthony"]), study=paste(unique(x[,"study"]), collapse=" ")))

write.table(allo, "/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/undercAllochthony_averaged.csv", sep=",", row.names=FALSE)

summary(lm(allochthony~as.factor(year)*lake, data=allo00[allo00[,"consumerGroup"]=="Chaoborus",]))

lme4(allochthony~lake+)


allo.sd <- ddply(allo, .variables=c("lake", "consumerGroup"), function(x)sd(x[,"allochthony"], na.rm=TRUE))
mean(allo.sd[1:8,3])