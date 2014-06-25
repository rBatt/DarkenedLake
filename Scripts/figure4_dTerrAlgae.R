

load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")

Save <- TRUE
SaveType <- ".png"

if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_dTerrAlgae.pdf", height=3, width=3.23)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_dTerrAlgae.png", units="in", res=600, height=3, width=3.23, type="cairo")}
	if(SaveType==".eps"){setEPS(); postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_dTerrAlgae.eps", height=3, width=3.23, pointsize=10)}
}else{
	dev.new(height=3, width=3, pointsize=10, family="Times")
}
par(mfrow=c(1,1), mar=c(2.1,2,0.1,0.1), oma=c(0,0,0,0), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")

aYlim <- c(min(Chosen25th[,"Algae"], na.rm=TRUE), max(Chosen75th[,"Algae"], na.rm=TRUE))*c(1.15,1.05)
tYlim <- c(min(Chosen25th[,"All.Terrestrial"], na.rm=TRUE), max(Chosen75th[,"All.Terrestrial"], na.rm=TRUE))*c(1.00,1.15)

plot(ChosenMedDiffs[,"Algae"], ChosenMedDiffs[,"All.Terrestrial"], xlim=aYlim, xlab="", xaxt="s", ylab="", ylim=tYlim, pch=NA)
abline(h=0, lty="dashed", col="gray")
abline(v=0, lty="dashed", col="gray")

arrows(x0=Chosen25th[,"Algae"], y0=ChosenMedDiffs[,"All.Terrestrial"], x1=Chosen75th[,"Algae"], y1=ChosenMedDiffs[,"All.Terrestrial"], angle=90, code=0, length=0.05, col="gray")
arrows(x0=ChosenMedDiffs[,"Algae"], y0=Chosen25th[,"All.Terrestrial"], x1=ChosenMedDiffs[,"Algae"], y1=Chosen75th[,"All.Terrestrial"], angle=90, code=0, length=0.05, col="gray")

points(ChosenMedDiffs[,"Algae"], ChosenMedDiffs[,"All.Terrestrial"], pch=20)

arrows(x0=-0.28, y0=0.39, x1=-0.315, y1=.32, length=0.075)
text(x=-0.28, y=0.395, "Complementarity", pos=4, font=2, offset=c(-0.1, 0))
ConsNames <- c("U. limi", "S. oregonensis", "Phoxinus spp.", "A. melas", "P. promelas", "Chaoborus spp.", "H. trivolvis")
text(x=ChosenMedDiffs[,"Algae"][-c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][-c(2,4,6)], labels=ConsNames[-c(2,4,6)], srt=0, xpd=NA,adj=c(1.05,-0.1), cex=1, font=3, ps=9)
text(x=ChosenMedDiffs[,"Algae"][c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][c(2,4,6)], labels=ConsNames[c(2,4,6)], srt=0, xpd=NA,adj=c(-0.01,-0.2), cex=1, font=3, ps=9)
mtext("Change in terrestrial support", side=2, line=1.1, cex=1)
mtext("Change in algal support", side=1, line=1.1, cex=1)
abline(a=0, b=-1)
if(Save){dev.off()}

