
rm(list=ls())
graphics.off()
setwd("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load("Cons_Mixture_Ward2010&2012_v0.4.7.RData")
FigureFolder <- paste("Figures/Figures_", "v0.4.7", sep="")
SaveType <- c(".pdf", ".png", ".eps")[2]
library("plyr")

# Save <- c(TRUE, FALSE)[1]
# SaveType <- c(".pdf", ".png")[2]

# ===============
# = New _v0.4.5 =
# ===============
setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("dTerrAlgae_", Version, sep=""), ".pdf", sep=""), height=3, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("dTerrAlgae_", Version, sep=""), ".png", sep=""), units="in", res=300, height=3, width=3.23)}
	if(SaveType==".eps"){setEPS(); postscript(file=paste(paste("dTerrAlgae_", Version, sep=""), ".eps", sep=""), height=3, width=3.23, pointsize=10)}
}else{
	dev.new(height=4, width=3, pointsize=10, family="Times")
}
par(mfrow=c(1,1), mar=c(2.1,2,0.1,0.1), oma=c(0,0,0,0), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")

aYlim <- c(min(Chosen25th[,"Algae"], na.rm=TRUE), max(Chosen75th[,"Algae"], na.rm=TRUE))*c(1.15,1.05)
tYlim <- c(min(Chosen25th[,"All.Terrestrial"], na.rm=TRUE), max(Chosen75th[,"All.Terrestrial"], na.rm=TRUE))

plot(ChosenMedDiffs[,"Algae"], ChosenMedDiffs[,"All.Terrestrial"], xlim=aYlim, xlab="", xaxt="s", ylab="", ylim=tYlim, pch=NA)
abline(h=0, lty="dashed", col="gray")
abline(v=0, lty="dashed", col="gray")

arrows(x0=Chosen25th[,"Algae"], y0=ChosenMedDiffs[,"All.Terrestrial"], x1=Chosen75th[,"Algae"], y1=ChosenMedDiffs[,"All.Terrestrial"], angle=90, code=0, length=0.05, col="gray")
arrows(x0=ChosenMedDiffs[,"Algae"], y0=Chosen25th[,"All.Terrestrial"], x1=ChosenMedDiffs[,"Algae"], y1=Chosen75th[,"All.Terrestrial"], angle=90, code=0, length=0.05, col="gray")

points(ChosenMedDiffs[,"Algae"], ChosenMedDiffs[,"All.Terrestrial"], pch=20)

ConsNames <- c("U. limi", "S. oregonensis", "Phoxinus spp.", "A. melas", "P. promelas", "Chaoborus spp.", "H. trivolvis")
text(x=ChosenMedDiffs[,"Algae"][-c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][-c(2,4,6)], labels=ConsNames[-c(2,4,6)], srt=0, xpd=NA,adj=c(1.05,-0.1), cex=1, font=3, ps=9)
text(x=ChosenMedDiffs[,"Algae"][c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][c(2,4,6)], labels=ConsNames[c(2,4,6)], srt=0, xpd=NA,adj=c(-0.01,-0.2), cex=1, font=3, ps=9)
mtext("Change in terrestrial use", side=2, line=1.1, cex=1)
mtext("Change in algal use", side=1, line=1.1, cex=1)
abline(a=0, b=-1)
if(Save){dev.off()}

# ===================
# = END new _v0.4.5 =
# ===================


# ===============
# = New _v0.4.7 =
# ===============
cI10 <- ResourceUse[,"Consumer"]=="Chaoborus" & ResourceUse[,"Year"]==2010 #index of 2010 chaob posterior
cI12 <- ResourceUse[,"Consumer"]=="Chaoborus" & ResourceUse[,"Year"]==2012 #index of 2012 chaob posterior
fN <- rev(c("FHM", "DAC", "CMM", "BHD1")) #fish names
# fC <- tim.colors(n=18, alpha=1)[4*c(0.75,2.5,3.25,4)]
fC <- rep("black",4)

# dev.new(width=3.5, height=5)
setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
if(Save){
	if(SaveType==".pdf"){pdf(paste("deltaFishChaob_", Version, ".pdf", sep=""), width=3.23, height=4.5)}
	if(SaveType==".png"){png(paste("deltaFishChaob_", Version, ".png", sep=""), width=3.23, height=4.5, units="in", res=500)}
}

par(mfrow=c(1,1), mar=c(2.5,2.3,0,0), oma=c(0,0,0.2,0.2), ps=9, las=1, tcl=-0.25, mgp=c(3,0.35,0), yaxp=c(0,20,10), family="Times", cex=1)
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
dev.off()

# ===================
# = END New _v0.4.7 =
# ===================


NeatConsNames <- c("Calanoid"="S. oregonensis", "Mesocyclops"="Mesocyclops", "Chaoborus"="Chaoborus spp.", "Helisoma trivolvis"="H. trivolvis", "FHM"="P. promelas", "DAC"="Phoxinus spp.", "BHD1"="A. melas", "BHD2"="A. melas", "CMM"= "U. limi", "PKS"="L. gibbosus", "YWP"="P. flavescens")
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("NeatSummary_", Version, sep=""), ".pdf", sep=""), height=6.5, width=6.81)}
	if(SaveType==".png"){png(file=paste(paste("NeatSummary_", Version, sep=""), ".png", sep=""), units="in", res=200, height=6.5, width=6.81)}
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


# ===============
# = ZOOPLANKTON =
# ===============
zData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardZoopMass2010&2012.csv")
cData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardChaobMass2010&2012.csv")
zData <- reshape(zData0, varying=list(names(zData0[,4:17])), times=names(zData0[,4:17]), ids=1:nrow(zData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(zData) <- NULL
zData <- subset(zData, DoY>=143)

zData[,"Year"] <- as.factor(zData[,"Year"])

cData <- reshape(cData0, varying=list(names(cData0[,4:7])), times=names(cData0[,4:7]), ids=1:nrow(cData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
row.names(cData) <- NULL
cData <- subset(cData, DoY>=143)
cData[,"Year"] <- as.factor(cData[,"Year"])

sumzData <- aggregate(zData[,"Mass"], by=list(zData[,"Lake"], zData[,"Year"], zData[,"DoY"]), sum)
names(sumzData) <- c("Lake", "Year", "DoY", "Mass")
zYearMean <- aggregate(sumzData[,"Mass"], by=list(sumzData[,"Lake"], sumzData[,"Year"]), mean)
names(zYearMean) <- c("Lake", "Year", "Mass")
sumcData <- aggregate(cData[,"Mass"], by=list(cData[,"Lake"], cData[,"Year"], cData[,"DoY"]), sum)
names(sumcData) <- c("Lake", "Year", "DoY", "Mass")
cYearMean <- aggregate(sumcData[,"Mass"], by=list(sumcData[,"Lake"], sumcData[,"Year"]), mean)
names(cYearMean) <- c("Lake", "Year", "Mass")


if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".pdf", sep=""), height=4.75, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".png", sep=""), units="in", res=200, height=4.75, width=3.23)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".eps", sep=""), height=4.75, width=3.23)}
}else{
	dev.new(height=4.75, width=3.23)
}
par(mfrow=c(2,1), mar=c(1,2.5,0,0), oma=c(0.3, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(Zooplankton~(g/m^2)), side=2, line=1.5, cex=1)

boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(italic(Chaoborus)~spp.~(g/m^2)), side=2, line=1.5, cex=1)
dev.off()



# ============
# = ROUTINES =
# ============
#Taken from L&W_Routine_2010&2012_v3.r and then revised
Data0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
POC0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/POC_PaulWard2010&2012.csv")
ChlaP0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_ChlaProfile_2010&2012.csv")
ChlaP0 <- subset(ChlaP0, Zid!="Meta")
areal <- function(x){
	mean(x[,"Chla"])*max(x[,"Z"])
}
aChl0 <- ddply(ChlaP0, .variables=c("Lake", "Year", "Date"), .fun=areal)
names(aChl0) <- c("Lake", "Year", "Date", "aChla")

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")
PlotNames <- c("L 10", "L 12", "W 10", "W 12")
PML_Data0 <- subset(Data0, Layer=="PML")[,]
PML_Data <- merge(PML_Data0, POC0, all.x=TRUE)
PML_Data[,"C:Chl"] <- PML_Data[,"POC"]/(PML_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "POC","PON","C:Chl", "Temp", "DOC", "DIC", "DO_Conc")
Meta_Data0 <- subset(Data0, Layer=="Meta")[,]
Meta_Data <- merge(Meta_Data0, POC0, all.x=TRUE)
Meta_Data[,"C:Chl"] <- Meta_Data[,"POC"]/(Meta_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "Temp", "DOC", "DIC", "DO_Conc")
Hypo_Data <- subset(Data0, Layer=="Hypo")[,]

VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")

# =============================================
# = 4-panel, areal chlorophyll, pCO2, DO, NEP =
# =============================================

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("LimnoMetab_PaulWard_2010&2012_", Version, sep=""), ".pdf", sep=""), height=7, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("LimnoMetab_PaulWard_2010&2012_", Version, sep=""), ".png", sep=""), units="in", res=200, height=7, width=3.23)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("LimnoMetab_PaulWard_2010&2012_", Version, sep=""), ".eps", sep=""), height=7, width=3.23)}
}else{
	dev.new(height=7, width=3.23)
}
par(mfrow=c(4,1), mar=c(1,2.5,0,0), oma=c(0.2, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)

#Chlorophyll
boxplot(aChla~Year+Lake, data=aChl0, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(Chlorophyll~(mg/m^2)), side=2, line=1.5, cex=1)

#pCO2
ppCO2 <- subset(Data0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(italic(p)*CO[2]~(mu*atm)), side=2, line=1.5, cex=1)

#DO
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_2010&2012_Metabolism_v0.2.RData")
boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(DO~("%"*saturation)), side=2, line=1.5, cex=1)

#NEP
pWardPaul_Metabolism <- subset(WardPaul_Metabolism, DoY>=143)
pWardPaul_Metabolism[,"Lake"] <- relevel(pWardPaul_Metabolism[,"Lake"], ref="Paul")
boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(NEP~(mmol~O[2]~m^-3~d^-1)), side=2, line=1.5, cex=1)
dev.off()

# aggregate(pWardPaul_Metabolism[,"NEP"], by=list(pWardPaul_Metabolism[,"Lake"], pWardPaul_Metabolism[,"Year"]), FUN=mean, na.rm=TRUE)

# DOsummary <- aggregate(AllDO[,"MeanDO"], by=list(AllDO[,"Lake"], AllDO[,"Year"]), FUN=mean, na.rm=TRUE)
# cbind(DOsummary, "localDO"=DOsummary[,"x"]+(100-94.6))

