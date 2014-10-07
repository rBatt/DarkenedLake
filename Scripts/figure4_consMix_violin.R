

library(DescTools)
library(png)
library(plyr)
source("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/myPlotViolin.default.R")

load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")

Save <- TRUE
SaveType <- ".png"



chaob.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/Chaob.png")
BHD.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/black_bullhead.png")
CMM.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/central_mudminnow.png")
NRD.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/northern_redbelly_dace.png")
calan.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/Soregonensis.png")
htri.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/Helisoma.png")
FHM.png <- readPNG("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/critterPics/FHM.png")

addPNG <- function(ping, xl, xr, yt){
	p.dim <- dim(ping)
	din <- par("din")
	usr <- par("usr")
	fig <- par("fig")
	fin <- par("fin")
	mai <- par("mai")

	region.width <- fin[2] - (mai[2] + mai[4])
	region.height <- fin[1] - (mai[1] + mai[3])

	fig.width <- usr[2] - usr[1]
	fig.height <- usr[4] - usr[3]

	
	# hOw <- (p.dim[1]/p.dim[2]) * (fig.height/fig.width)
	# hOw <- (p.dim[1]/p.dim[2]) * (fig.height/fig.width) * ((fig[2]-fig[1])/(fig[4]-fig[3]))
	# hOw <- (p.dim[1]/p.dim[2]) * (fig.height/fig.width) * (region.height/region.width)
	hOw <- (p.dim[1]/p.dim[2]) * (fig.height/fig.width) * (fin[1]/fin[2])
	
	xw <- xr-xl
	yh <- xw*hOw
	yb <- yt - yh
		
	# (yt-yb) / xw
	
	rasterImage(ping, xleft=xl, ybottom=yb, xright=xr, ytop=yt)
	
}

if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_consMix_Violin.pdf", height=6.204112, width=6.1)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_consMix_Violin.png", units="in", res=300, height=6.204112, width=6.1, type="quartz")}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig4_consMix_Violin.eps", height=6.204112, width=6.1)}
}else{
	dev.new(height=6.26927, width=6.1)
}
layout(matrix(c(1,2,3,4,4,6,5,5,7), ncol=3, byrow=TRUE), widths=c(2.4/7, 2/7, 2.6/7, 2/7, 2/7, 3/7, 2/7, 2/7, 3/7))
# par(mar=c(2.5,0.5,1,0.5), oma=c(0,2,0,0)), cex=1)
par(oma=c(0, 0, 0, 0.2))
for(i in 1:length(Cons)){
	Yaxt <- ifelse(is.element(i, c(1,4,5)), "s", "n")
	# LegPos <- ifelse(is.element(i, c(1,2,4,5)), "topleft", "topright")
	LegPos <- switch(i,
		"topleft",
		"topleft",
		"htri",
		"lfish1",
		"lfish2",
		"rfish",
		"rfish2"
		)
	if(!is.element(i, c(5,7))){
		if(Yaxt=="s"){
			# par(mar=c(2,3,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
			par(mar=c(2,2.5,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}else{
			par(mar=c(2,1,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}
	}else{
		if(Yaxt=="s"){
			# par(mar=c(1.75,3,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
			par(mar=c(1.75,2.5,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}else{
			par(mar=c(1.75,1,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}
	}
	NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	ThisRU0 <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	ThisRU00 <- reshape(ThisRU0, varying=list(c(NamesThisUse)), times=ConsChoices[[Cons[i]]][[GroupChoose]], ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")

	ThisRU <- ThisRU00[,c("Year", "Month", "Consumer", "Grouping", "Source", "Proportion")]
	row.names(ThisRU) <- NULL
	ResourceNames <- ConsChoicesShort[ConsChoices[[Cons[i]]][[GroupChoose]]]
	
	ThisRU[,"Source"] <- factor(ThisRU[,"Source"], levels=ConsChoices[[Cons[i]]][[GroupChoose]], ordered=TRUE)
	
	
	ThisAt_Axis <- c(1.5,3.5,5.5,7.5)[1:length(ResourceNames)]
# PlotViolin(formula=Proportion~Year+Source, data=ThisRU, ylim=c(0,1), names=NA, xaxt="n", col=c(NA, "lightgray"), args.boxplot=list(boxwex=0.05, border=c("black","black"), col=c("black", "black"), show.names=FALSE, outline=FALSE,lwd=1.25))
PlotViolin(formula=Proportion~Year+Source, data=ThisRU, ylim=c(0,1.15), names=NA, xaxt="n", col=c(NA, NA), border=c("red2", "blue2"), lwd=1.15, args.boxplot=list(boxwex=0.125, border=c("red2","blue2"), col=c(NA, NA), show.names=FALSE, outline=FALSE,lwd=1.25))
axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
	if(Yaxt=="n"){axis(side=2, labels=FALSE)}
	X <- c("topleft"=0.15, "topright"=ThisAt_Axis[length(ResourceNames)]+1, "htri"=5.0, "lfish1"=8.5, "lfish2"=7.25, "rfish"=5, "rfish2"=3.75)
	Pos <- c("topleft"=4, "topright"=2, "top"=2, "htri"=2, "lfish1"=2, "lfish2"=2, "rfish"=2, "rfish2"=2)
	text(x=X[LegPos], y=1.14, labels=NeatConsNames[Cons[i]], pos=Pos[LegPos], font=3)
	switch(Cons[i],
		"Calanoid" = addPNG(calan.png, xl=0.3, xr=3, yt=1.17),
		"Chaoborus"=addPNG(chaob.png, xl=0.4, xr=3.15, yt=1.125),
		"Helisoma trivolvis" = addPNG(htri.png, xl=3.25, 4.75, yt=1.12),
		"FHM" = addPNG(FHM.png, xl=4.75, xr=7.75, yt=1.20),
		"DAC" = addPNG(NRD.png, xl=3.5, xr=6.5, yt=1.19),
		"CMM" = addPNG(CMM.png, xl=1.65, xr=5.35, yt=1.17),
		"BHD1" = addPNG(BHD.png, xl=1.95, xr=5.00, yt=1.20)
		)
	if(i==length(Cons)){
		mtext("Proportion of Diet", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}

# ===============
# = Print Means =
# ===============
for(i in 1:length(Cons)){
	NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	myRU <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	myRU <- reshape(myRU, varying=list(c(NamesThisUse)), times=ConsChoices[[Cons[i]]][[GroupChoose]], ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")
	print(ddply(myRU, c("Year", "Consumer", "Grouping", "Source"), function(x)data.frame("Proportion"=mean(x[,"Proportion"]))))
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



