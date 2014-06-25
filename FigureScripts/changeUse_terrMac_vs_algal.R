
# ===============
# = New _v0.4.5 =
# ===============
# setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))



tm25 <- rowSums(cbind(Chosen25th[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)
tmMed <- rowSums(cbind(ChosenMedDiffs[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)
tm75 <- rowSums(cbind(Chosen75th[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)

if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/dAlgae_vs_terrMac.pdf", height=3, width=3.23)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/dAlgae_vs_terrMac.png", units="in", res=300, height=3, width=3.23)}
	if(SaveType==".eps"){setEPS(); postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/dAlgae_vs_terrMac.eps", height=3, width=3.23, pointsize=10)}
}else{
	dev.new(height=4, width=3, pointsize=10, family="Times")
}
par(mfrow=c(1,1), mar=c(2.1,2,0.1,0.1), oma=c(0,0,0,0), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")

aYlim <- c(min(Chosen25th[,"Algae"], na.rm=TRUE), max(Chosen75th[,"Algae"], na.rm=TRUE))*c(1.15,1.05)
tYlim <- c(min(tm25, na.rm=TRUE), max(tm75, na.rm=TRUE))

plot(ChosenMedDiffs[,"Algae"], tmMed, xlim=aYlim, xlab="", xaxt="s", ylab="", ylim=tYlim, pch=NA)
abline(h=0, lty="dashed", col="gray")
abline(v=0, lty="dashed", col="gray")

arrows(x0=Chosen25th[,"Algae"], y0=tmMed, x1=Chosen75th[,"Algae"], y1=tmMed, angle=90, code=0, length=0.05, col="gray")
arrows(x0=ChosenMedDiffs[,"Algae"], y0=tm25, x1=ChosenMedDiffs[,"Algae"], y1=tm75, angle=90, code=0, length=0.05, col="gray")

points(ChosenMedDiffs[,"Algae"], tmMed, pch=20)

ConsNames <- c("U. limi", "S. oregonensis", "Phoxinus spp.", "A. melas", "P. promelas", "Chaoborus spp.", "H. trivolvis")
text(x=ChosenMedDiffs[,"Algae"][-c(1,3)], y=tmMed[-c(1,3)], labels=ConsNames[-c(1,3)], srt=0, xpd=NA,adj=c(1.05,-0.1), cex=1, font=3, ps=9)
text(x=ChosenMedDiffs[,"Algae"][c(1,3)], y=tmMed[c(1,3)], labels=ConsNames[c(1,3)], srt=0, xpd=NA,adj=c(-0.01,-0.2), cex=1, font=3, ps=9)
mtext("Change in terrestrial + macrophyte use", side=2, line=1.1, cex=1)
mtext("Change in algal use", side=1, line=1.1, cex=1)
if(Save){dev.off()}
