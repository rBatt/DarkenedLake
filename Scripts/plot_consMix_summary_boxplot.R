

Save <- FALSE
load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")


NeatConsNames <- c("Calanoid"="S. oregonensis", "Mesocyclops"="Mesocyclops", "Chaoborus"="Chaoborus spp.", "Helisoma trivolvis"="H. trivolvis", "FHM"="P. promelas", "DAC"="Phoxinus spp.", "BHD1"="A. melas", "BHD2"="A. melas", "CMM"= "U. limi", "PKS"="L. gibbosus", "YWP"="P. flavescens")
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("NeatSummary_", Version, sep=""), ".pdf", sep=""), height=6.5, width=6.81)}
	if(SaveType==".png"){png(file=paste(paste("NeatSummary_", Version, sep=""), ".png", sep=""), units="in", res=600, height=6.5, width=6.81)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("NeatSummary_", Version, sep=""), ".eps", sep=""), height=7, width=6.81)}
}else{
	dev.new(height=7, width=6.811)
}
layout(matrix(c(1,2,3,4,4,6,5,5,7), ncol=3, byrow=TRUE), widths=c(2.4/7, 2/7, 2.6/7, 2/7, 2/7, 3/7, 2/7, 2/7, 3/7))
# par(mar=c(2.5,0.5,1,0.5), oma=c(0,2,0,0)), cex=1)
par(oma=c(0, 0, 0, 0.2))
for(i in 1:length(Cons)){
	Yaxt <- ifelse(is.element(i, c(1,4,5)), "s", "n")
	LegPos <- ifelse(is.element(i, c(1,2,4,5)), "topleft", "topright")
	if(!is.element(i, c(5,7))){
		if(Yaxt=="s"){
			par(mar=c(2,3,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}else{
			par(mar=c(2,1,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}
	}else{
		if(Yaxt=="s"){
			par(mar=c(1.75,3,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}else{
			par(mar=c(1.75,1,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}
	}
	# LegPos <- c("topleft", "topleft", "topright", "topright", "topleft", "topright", "topright")[i]
	NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	ThisRU0 <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	ThisRU00 <- reshape(ThisRU0, varying=list(c(NamesThisUse)), times=ConsChoices[[Cons[i]]][[GroupChoose]], ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")

	ThisRU <- ThisRU00[,c("Year", "Month", "Consumer", "Grouping", "Source", "Proportion")]
	row.names(ThisRU) <- NULL
	ResourceNames <- ConsChoicesShort[ConsChoices[[Cons[i]]][[GroupChoose]]]
	
	ThisRU[,"Source"] <- factor(ThisRU[,"Source"], levels=ConsChoices[[Cons[i]]][[GroupChoose]], ordered=TRUE)
	
	# aggregate(ThisRU[,"Proportion"], by=list("Year"=ThisRU[,1], "Source"=ThisRU[,"Source"]), FUN=median)
	
	ThisAt_Axis <- c(0.5,3.5,6.5,9.5)[1:length(ResourceNames)]
	# boxplot(Proportion~Year+Source, data=ThisRU, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5, yaxt=Yaxt)
	boxplot(Proportion~Year+Source, data=ThisRU, border=c("black","black"), col=c(NA, "lightgray"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.25, yaxt=Yaxt)
	axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
	if(Yaxt=="n"){axis(side=2, labels=FALSE)}
	# mtext(ConsNameMedium[Cons[i]], side=2, line=2)
	# legend(LegPos, NeatConsNames[Cons[i]])
	X <- c("topleft"=-0.9, "topright"=ThisAt_Axis[length(ResourceNames)]+1)
	Pos <- c("topleft"=4, "topright"=2)
	text(x=X[LegPos], y=0.95, labels=NeatConsNames[Cons[i]], pos=Pos[LegPos], font=3)
	if(i==length(Cons)){
		mtext("Proportion of Diet", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}