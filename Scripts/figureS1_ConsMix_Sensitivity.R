
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Orgd_ConsMix_Sensit.RData")
load("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/Cons_Mixture_Ward2010&2012.RData")
library(plyr)

NeatConsNames <- c("Calanoid"="S. oregonensis", "Mesocyclops"="Mesocyclops", "Chaoborus"="Chaoborus spp.", "Helisoma trivolvis"="H. trivolvis", "FHM"="P. promelas", "DAC"="Phoxinus spp.", "BHD1"="A. melas", "BHD2"="A. melas", "CMM"= "U. limi", "PKS"="L. gibbosus", "YWP"="P. flavescens")

Save <- TRUE
SaveType <- c(".png", ".tiff")[2]

# ========================
# = Sensitivity Box Plot =
# ========================
if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS1_Sensitivity_boxplot.pdf", height=6.204112, width=6.1)}
	if(SaveType==".png"){png(file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS1_Sensitivity_boxplot.png", units="in", res=300, height=6.204112, width=6.1)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("Sensitivity", Version, sep=""), ".eps", sep=""), height=7, width=6.81)}
	if(SaveType==".tiff"){tiff("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS1_Sensitivity_boxplot.tiff", width=6.1, height=6.204112, units="in", res=300, compression="lzw", type="quartz")}
}else{
	dev.new(height=6.204112, width=6.1)
}
layout(matrix(c(1,2,3,4,4,6,5,5,7), ncol=3, byrow=TRUE), widths=c(2.4/7, 2/7, 2.6/7, 2/7, 2/7, 3/7, 2/7, 2/7, 3/7))
# par(mar=c(2.5,0.5,1,0.5), oma=c(0,2,0,0)), cex=1)
par(oma=c(0, 0, 0, 0.2))
for(i in 1:length(Cons)){
	Yaxt <- ifelse(is.element(i, c(1,4,5)), "s", "n")
	LegPos <- ifelse(is.element(i, c(1,2,4,5)), "topleft", "topright")
	if(!is.element(i, c(5,7))){
		if(Yaxt=="s"){
			par(mar=c(2,3,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}else{
			par(mar=c(2,1,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}
	}else{
		if(Yaxt=="s"){
			par(mar=c(1.75,3,0.2,0),cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}else{
			par(mar=c(1.75,1,0.2,0),cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
		}
	}
	# LegPos <- c("topleft", "topleft", "topright", "topright", "topleft", "topright", "topright")[i]
	t.names <- ConsChoices[[Cons[i]]][[GroupChoose]]#[!is.na(ConsChoices[[i]][[GroupChoose]])]
	NamesThisUse <- make.names(t.names)
	# NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])
	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	ThisRU0 <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	ThisRU00 <- reshape(ThisRU0, varying=list(c(NamesThisUse)), times=t.names, ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")

	ThisRU <- ThisRU00[,c("Year", "Month", "Consumer", "Grouping", "Source", "Proportion")]
	row.names(ThisRU) <- NULL
	ResourceNames <- ConsChoicesShort[ConsChoices[[Cons[i]]][[GroupChoose]]] #[!is.na(ConsChoices[[i]][[GroupChoose]])]
	tSen <- sensit[[Cons[[i]]]][,NamesThisUse]
	tSen2 <- sensit2[[Cons[[i]]]][,NamesThisUse]
	
	ThisRU[,"Source"] <- factor(ThisRU[,"Source"], levels=t.names, ordered=TRUE)
		
	ThisAt_Axis <- c(0.5,3.5,6.5,9.5)[1:length(ResourceNames)]
	# boxplot(Proportion~Year+Source, data=ThisRU, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5, yaxt=Yaxt)
	boxplot(Proportion~Year+Source, data=ThisRU, border=c("black","black"), col=c(NA, "lightgray"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.25, yaxt=Yaxt)
	xVals <- rep(ThisAt_Axis,each=2)+c(-.5, .5)
	points(xVals, tSen, pch=21, cex=1, lwd=2, bg="transparent", col="black")
	if(is.element(Cons[[i]], c("FHM", "DAC", "CMM", "BHD1"))){
		points(xVals, tSen2, pch=22, cex=1.3, lwd=2, bg="transparent", col="black")
	}

	
	axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
	if(Yaxt=="n"){axis(side=2, labels=FALSE)}
	# mtext(ConsNameMedium[Cons[i]], side=2, line=2)
	# legend(LegPos, NeatConsNames[Cons[i]])
	X <- c("topleft"=-0.9, "topright"=ThisAt_Axis[length(ResourceNames)]+1)
	Pos <- c("topleft"=4, "topright"=2)
	text(x=X[LegPos], y=0.95, labels=NeatConsNames[Cons[i]], pos=Pos[LegPos], font=3)
	if(i==length(Cons)){
		mtext("Proportion of Composition", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}

# ===============
# = Print Means =
# ===============
for(i in 1:length(Cons)){
	t.names <- ConsChoices[[Cons[i]]][[GroupChoose]] #[!is.na(ConsChoices[[i]][[GroupChoose]])]
	NamesThisUse <- make.names(t.names)
	# NamesThisUse <- make.names(ConsChoices[[Cons[i]]][[GroupChoose]])

	AllNamesThisuse <- c("Year", "Month", "Consumer", "Grouping", NamesThisUse)
	myRU <- droplevels(subset(ResourceUse, Consumer==Cons[i] & Grouping==GroupChoose & Month=="Pooled"))[AllNamesThisuse] #try to subset only the sources that were part of the analysis for this consumer
	myRU <- reshape(myRU, varying=list(c(NamesThisUse)), times=t.names, ids=1:nrow(ThisRU0), timevar="Source", v.names="Proportion", direction="long")
	print(ddply(myRU, c("Year", "Consumer", "Grouping", "Source"), function(x)data.frame("Proportion"=mean(x[,"Proportion"]))))
}
