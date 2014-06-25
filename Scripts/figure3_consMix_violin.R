

library(DescTools)
source("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Scripts/myPlotViolin.default.R")


Save <- FALSE
load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")




if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/FiguresconsMix_Violin.pdf", height=6.5, width=6.81)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/FiguresconsMix_Violin.png", units="in", res=600, height=6.5, width=6.81)}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/FiguresconsMix_Violin.eps", height=7, width=6.81)}
}else{
	dev.new(height=7, width=6.811)
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
		"lfish",
		"rfish",
		"lfish",
		"rfish"
		)
	if(!is.element(i, c(5,7))){
		if(Yaxt=="s"){
			par(mar=c(2,3,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}else{
			par(mar=c(2,1,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}
	}else{
		if(Yaxt=="s"){
			par(mar=c(1.75,3,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
		}else{
			par(mar=c(1.75,1,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25, yaxt=Yaxt)
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
PlotViolin(formula=Proportion~Year+Source, data=ThisRU, ylim=c(0,1.15), names=NA, xaxt="n", col=c(NA, "lightblue"), args.boxplot=list(boxwex=0.05, border=c("black","blue"), col=c("black", "blue"), show.names=FALSE, outline=FALSE,lwd=1.25))
axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
	if(Yaxt=="n"){axis(side=2, labels=FALSE)}
	X <- c("topleft"=0.15, "topright"=ThisAt_Axis[length(ResourceNames)]+1)
	Pos <- c("topleft"=4, "topright"=2, "top"=1)
	text(x=X[LegPos], y=0.95, labels=NeatConsNames[Cons[i]], pos=Pos[LegPos], font=3)
	if(i==length(Cons)){
		mtext("Proportion of Diet", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}



