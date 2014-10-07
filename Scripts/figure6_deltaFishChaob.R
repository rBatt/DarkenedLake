

load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")
Save <- TRUE
SaveType <- ".png"

# ===============
# = New _v0.4.7 =
# ===============
cI10 <- ResourceUse[,"Consumer"]=="Chaoborus" & ResourceUse[,"Year"]==2010 #index of 2010 chaob posterior
cI12 <- ResourceUse[,"Consumer"]=="Chaoborus" & ResourceUse[,"Year"]==2012 #index of 2012 chaob posterior
fN <- rev(c("FHM", "DAC", "CMM", "BHD1")) #fish names
# fC <- tim.colors(n=18, alpha=1)[4*c(0.75,2.5,3.25,4)]
fC <- rep("black",4)

# dev.new(width=3.5, height=5)
if(Save){
	if(SaveType==".pdf"){pdf("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig6_deltaFishChaob.pdf", width=2.9, height=4.040248)}
	if(SaveType==".png"){png("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig6_deltaFishChaob.png", width=2.9, height=4.040248, units="in", res=200, type="quartz")}
}else{
	dev.new(width=2.9, height=4.040248)
}

par(mfrow=c(1,1), mar=c(2.5,2.3,0,0), oma=c(0,0,0.2,0.2), ps=8, las=1, tcl=-0.25, mgp=c(3,0.35,0), yaxp=c(0,20,10), family="Times", cex=1)
r <- c("All.Phytoplankton", "All.Terrestrial")[1]
pC10 <- ResourceUse[cI10,r] #phytoplankton for chaoborus in 2010
pC12 <- ResourceUse[cI12,r] #phytoplankton for chaoborus in 2012
# plot(density(pC10), xlim=c(0,1), ylim=c(0,7), zero.line=FALSE, main="")
for(f in 1:length(fN)){
	vertOff <- 5*(f-1)
	#2010 distribution of differences
	fI10 <- ResourceUse[,"Consumer"]==fN[f] & ResourceUse[,"Year"]==2010
	pF10 <- ResourceUse[fI10,r]
	d10 <- density(pF10-pC10)
	dx10 <- d10$x
	dy10 <- d10$y + vertOff

	#2012 distribution of differences
	fI12 <- ResourceUse[,"Consumer"]==fN[f] & ResourceUse[,"Year"]==2012
	pF12 <- ResourceUse[fI12,r]
	d12 <- density(pF12-pC12)
	dx12 <- d12$x
	dy12 <- d12$y + vertOff

	h <- max(c(dy10, dy12))+0.25

	#graph 2010
	if(f==1){
		plot(dx10,dy10, col=fC[f],xlim=c(-1,0.5), ylim=c(0,20), pch=NA, yaxt="n")
		# axis(side=2, at=seq(0,20, 0.75), labels=FALSE, tcl=0.1)
		tksAt <- c(0,1,2,3)
		nt <- length(tksAt)
		ats <- rep((5*((1:4)-1)), each=nt) + rep(tksAt, 4)
		labs <- ats - rep((5*((1:4)-1)), each=nt)
		axis(side=2, at=ats, labels=labs, tcl=-0.25, mgp=c(3,0.5,0))
		abline(h=0, col="darkgray", lty="dashed")
		segments(x0=0, x1=0, y0=vertOff, y1=(h), lty="dashed", col="darkgray")
		lines(dx10,dy10, col=fC[f], type="l")
	}else{
		lines(dx10,dy10, col=fC[f])
		abline(h=vertOff, col="darkgray", lty="dashed")
		segments(x0=0, x1=0, y0=vertOff, y1=(h), lty="dashed", col="darkgray")
	}

	#2012 graph
	lines(dx12,dy12, col=fC[f], lwd=3)

	#means of distributions of differences, and arrows showing change between years
	mu10 <- mean(pF10-pC10)
	mu12 <- mean(pF12-pC12)
	arrows(x0=mu10, y0=h, x1=mu12,y1=h, length=0.075, col=fC[f], lwd=2)
	text(x=-1.1, y=h, labels=rev(c("P. promelas", "Phoxinus", "U. limi", "A. melas"))[f], font=3, col=fC[f], pos=4)
}
mtext(quote(Fish~phi1[Phyto]-Chaob.~phi1[Phyto]), side=1, line=1.5, outer=FALSE, cex=1)
mtext("Posterior Density", side=2, line=1.35, cex=1, las=0)
legend("topright", bty="n", legend=c("2010","2012"), lwd=c(1,2), inset=c(0,-0.025))
if(Save){dev.off()}
