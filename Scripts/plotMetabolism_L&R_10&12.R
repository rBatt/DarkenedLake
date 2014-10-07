


load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/kf.epi.good.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/WardPaulMetabolism.kf.RData")

ddply(kf.epi.good, c("lake","year"), function(x)colMeans(x[,c("GPP","R","NEP")])/(15.99*2)*1000)

# kf.epi.good <- WardPaulMetabolism.kf

addBlankMetab <- function(kf.epi.good){
	kf.epi.good <- kf.epi.good[kf.epi.good[,"doy"]>=137 & kf.epi.good[,"doy"]<=240, ]
	new.doy <- seq(137, 240)
	blanks <- !new.doy%in%kf.epi.good[,"doy"]
	r0 <- data.frame(lake=kf.epi.good[1,"lake"], "year"=kf.epi.good[1,"year"], "doy"=new.doy, GPP=NA, R=NA, NEP=NA, datetime=NA)
	r1 <- r0
	r1[!blanks, c("GPP")] <- kf.epi.good[,c("GPP")]
	r1[!blanks, c("R")] <- kf.epi.good[,c("R")]
	r1[!blanks, c("NEP")] <- kf.epi.good[,c("NEP")]
	# r1[!blanks, c("datetime")] <- as.POSIXct(kf.epi.good[,c("datetime")], tz="GMT")
	r1
}
kf.epi.good.blanks <- ddply(kf.epi.good, c("lake","year"), addBlankMetab)

ward10.good <- kf.epi.good.blanks[kf.epi.good.blanks[,"lake"]=="Ward"&kf.epi.good.blanks[,"year"]==2010,]
ward12.good <- kf.epi.good.blanks[kf.epi.good.blanks[,"lake"]=="Ward"&kf.epi.good.blanks[,"year"]==2012,]
paul10.good <- kf.epi.good.blanks[kf.epi.good.blanks[,"lake"]=="Paul"&kf.epi.good.blanks[,"year"]==2010,]
paul12.good <- kf.epi.good.blanks[kf.epi.good.blanks[,"lake"]=="Paul"&kf.epi.good.blanks[,"year"]==2012,]

dev.new()
par(mfrow=c(2,2), mar=c(2,2,0.5,0.5), mgp=c(1.4, 0.35, 0), tcl=-0.25, ps=10, family="Times")
plot(ward10.good[,c("doy","GPP")], xlab="", ylab="GPP", type="o", ylim=range(c(ward10.good[,"GPP"], ward12.good[,"GPP"]),na.rm=TRUE))
plot(ward12.good[,c("doy","GPP")], xlab="", ylab="GPP", type="o", ylim=range(c(ward10.good[,"GPP"], ward12.good[,"GPP"]),na.rm=TRUE))
plot(paul10.good[,c("doy","GPP")], xlab="", ylab="GPP", type="o", ylim=range(c(paul10.good[,"GPP"], paul12.good[,"GPP"]),na.rm=TRUE))
plot(paul12.good[,c("doy","GPP")], xlab="", ylab="GPP", type="o", ylim=range(c(paul10.good[,"GPP"], paul12.good[,"GPP"]),na.rm=TRUE))

