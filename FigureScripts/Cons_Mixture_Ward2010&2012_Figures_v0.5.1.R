# _v0.5.0 (09-Feb-2014): Cole was concerned about reporting only means of resource contribution. Specifically, he raised the issue of the mean hiding the shape of the posterior. While I think I will stick to reporting means in the text, it is true that it may be difficult (impossible?) to distinguish between bimodality/ uniform/ unimodality in a boxplot. So I am going to try presenting these results in a "violin plot". I have also turned this script into a general purpose figure script. It's still a little messy (loading .RData with same object names but different content, lots of working directory resets, etc), but it includes all of the figures except the photic depth aqua conc figure

# _v0.5.1 (12-Feb-2014): Isotope biplots in Latex in order to get the italicized delta.

rm(list=ls())
graphics.off()
library("plyr")
library("DescTools")
library("cluster")
library("tikzDevice")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load("Cons_Mixture_Ward2010&2012_v0.4.7.RData")
FigureFolder <- paste("Figures/Figures_", "v0.5.1", sep="")
SaveType <- c(".pdf", ".png", ".eps")[2]
Version <- "v0.4.7"
ThisVersion <- "v0.5.1"

Save <- c(TRUE, FALSE)[2]
PlotViolin.default <- function(x, ..., horizontal = FALSE, bw = "SJ", na.rm = FALSE, names = NULL, args.boxplot = NULL){
    vlnplt <- function(x, y, center, horizontal = FALSE, col = NA, 
        border = par("fg"), lty = 1, lwd = 1, density = NULL, 
        angle = 45, fillOddEven = FALSE, ...) {
        x <- c(x, rev(x))
        y <- c(y, -rev(y))
        y <- y + center
        if (horizontal == FALSE) {
            tmp = x
            x = y
            y = tmp
        }
        polygon(x = x, y = y, border = border, col = col, lty = lty, 
            lwd = lwd, density = density, angle = angle, fillOddEven = fillOddEven, 
            ...)
    }
    m <- match.call(expand.dots = FALSE)
    pars <- m$...[names(m$...)[!is.na(match(names(m$...), c("cex", 
        "cex.axis", "cex.lab", "cex.main", "cex.sub", "col.axis", 
        "col.lab", "col.main", "col.sub", "family", "font", "font.axis", 
        "font.lab", "font.main", "font.sub", "las", "tck", "tcl", 
        "xaxt", "xpd", "yaxt")))]]
    oldpar <- par(pars)
    on.exit(par(oldpar))
    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names)) 
        attributes(args)$names != ""
    else rep(FALSE, length = length(args))
    groups <- if (is.list(x)) 
        x
    else args[!namedargs]
    if (0 == (n <- length(groups))) 
        stop("invalid first argument")
    if (length(class(groups))) 
        groups <- unclass(groups)
    if (!missing(names)) 
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names"))) 
            attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }
    xvals <- matrix(0, nrow = 512, ncol = n)
    yvals <- matrix(0, nrow = 512, ncol = n)
    center <- 1:n
    for (i in 1:n) {
        if (na.rm) 
            xi <- na.omit(groups[[i]])
        else xi <- groups[[i]]
        tmp.dens <- density(xi, bw = bw, from=0, to=1)
        xvals[, i] <- tmp.dens$x
        yvals.needtoscale <- tmp.dens$y
        yvals.scaled <- 7/16 * yvals.needtoscale/max(yvals.needtoscale)
        yvals[, i] <- yvals.scaled
    }
    if (horizontal == FALSE) {
        xrange <- c(1/2, n + 1/2)
        yrange <- range(xvals)
    }
    else {
        xrange <- range(xvals)
        yrange <- c(1/2, n + 1/2)
    }
    plot.args <- m$...[names(m$...)[!is.na(match(names(m$...), 
        c("xlim", "ylim", "main", "xlab", "ylab", "panel.first", 
            "panel.last", "frame.plot", "add")))]]
    if (!"xlim" %in% names(plot.args)) 
        plot.args <- c(plot.args, list(xlim = xrange))
    if (!"ylim" %in% names(plot.args)) 
        plot.args <- c(plot.args, list(ylim = yrange))
    if (!"xlab" %in% names(plot.args)) 
        plot.args <- c(plot.args, list(xlab = ""))
    if (!"ylab" %in% names(plot.args)) 
        plot.args <- c(plot.args, list(ylab = ""))
    if (!"frame.plot" %in% names(plot.args)) 
        plot.args <- c(plot.args, list(frame.plot = TRUE))
    if (!"add" %in% names(plot.args)) 
        add <- FALSE
    else add <- plot.args$add
    if (!add) 
        do.call(plot, c(plot.args, list(x = 0, y = 0, type = "n", 
            axes = FALSE)))
    poly.args <- args[names(args)[!is.na(match(names(args), c("border", 
        "col", "lty", "lwd", "density", "angle", "fillOddEven")))]]
    poly.args <- lapply(poly.args, rep, length.out = n)
    for (i in 1:n) do.call(vlnplt, c(lapply(poly.args, "[", i), 
        list(x = xvals[, i]), list(y = yvals[, i]), list(center = center[i]), 
        list(horizontal = horizontal)))
    axes <- Coalesce(unlist(m$...[names(m$...)[!is.na(match(names(m$...), 
        c("axes")))]]), TRUE)
    if (axes) {
        xaxt <- Coalesce(unlist(m$...[names(m$...)[!is.na(match(names(m$...), 
            c("xaxt")))]]), TRUE)
        if (xaxt != "n") 
            if (horizontal == TRUE) 
                axis(1)
            else axis(1, at = 1:n, labels = names)
        yaxt <- Coalesce(unlist(m$...[names(m$...)[!is.na(match(names(m$...), 
            c("yaxt")))]]), TRUE)
        if (yaxt != "n") 
            if (horizontal == TRUE) 
                axis(2, at = 1:n, labels = names)
            else axis(2)
    }
    if (!identical(args.boxplot, NA)) {
        args1.boxplot <- list(col = "black", add = TRUE, boxwex = 0.05, 
            axes = FALSE, outline = FALSE, whisklty = 1, staplelty = 0, 
            medcol = "white")
        args1.boxplot[names(args.boxplot)] <- args.boxplot
        do.call(boxplot, c(list(x, horizontal = horizontal), 
            args1.boxplot))
    }
}


load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Figures_v0.4.0_ModelSelection/Sensitivity_v0.0.5.RData")

# ===============
# = New _v0.4.5 =
# ===============
setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("dTerrAlgae_", Version, sep=""), ".pdf", sep=""), height=3, width=3.23)}
	if(SaveType==".png"){png(file=paste(paste("dTerrAlgae_", Version, sep=""), ".png", sep=""), units="in", res=300, height=3, width=3.23)}
	if(SaveType==".eps"){setEPS(); postscript(file=paste(paste("dTerrAlgae_", Version, sep=""), ".eps", sep=""), height=3, width=3.23, pointsize=10)}
}else{
	dev.new(height=3, width=3, pointsize=10, family="Times")
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

arrows(x0=-0.28, y0=0.39, x1=-0.315, y1=.32, length=0.075)
text(x=-0.28, y=0.395, "Complementarity", pos=4, font=2, offset=c(-0.1, 0))
ConsNames <- c("U. limi", "S. oregonensis", "Phoxinus spp.", "A. melas", "P. promelas", "Chaoborus spp.", "H. trivolvis")
text(x=ChosenMedDiffs[,"Algae"][-c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][-c(2,4,6)], labels=ConsNames[-c(2,4,6)], srt=0, xpd=NA,adj=c(1.05,-0.1), cex=1, font=3, ps=9)
text(x=ChosenMedDiffs[,"Algae"][c(2,4,6)], y=ChosenMedDiffs[,"All.Terrestrial"][c(2,4,6)], labels=ConsNames[c(2,4,6)], srt=0, xpd=NA,adj=c(-0.01,-0.2), cex=1, font=3, ps=9)
mtext("Change in terrestrial support", side=2, line=1.1, cex=1)
mtext("Change in algal support", side=1, line=1.1, cex=1)
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
}else{
	dev.new(width=3.23, height=4.5)
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
if(Save){dev.off()}

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

# ========================
# = NEW _v0.5.0 (violin) =
# ========================

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("Violin_", Version, sep=""), ".pdf", sep=""), height=6.5, width=6.81)}
	if(SaveType==".png"){png(file=paste(paste("Violin_", Version, sep=""), ".png", sep=""), units="in", res=200, height=6.5, width=6.81)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("Violin_", Version, sep=""), ".eps", sep=""), height=7, width=6.81)}
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
	
	ThisAt_Axis <- c(1.5,3.5,5.5,7.5)[1:length(ResourceNames)]
	# boxplot(Proportion~Year+Source, data=ThisRU, col=c("#FA807225","#3A5FCD25"), border=c("red","blue"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.5, yaxt=Yaxt)
	# boxplot(Proportion~Year+Source, data=ThisRU, border=c("black","black"), col=c(NA, "lightgray"), at=rep(ThisAt_Axis,each=2)+c(-.5, .5), show.names=FALSE, outline=FALSE, ylim=c(0,1), lwd=1.25, yaxt=Yaxt)
PlotViolin(formula=Proportion~Year+Source, data=ThisRU, ylim=c(0,1), names=NA, xaxt="n", col=c(NA, "lightgray"), args.boxplot=list(boxwex=0.05, border=c("black","black"), col=c("black", "black"), show.names=FALSE, outline=FALSE,lwd=1.25))
axis(side=1, at=ThisAt_Axis, labels=ResourceNames)
	if(Yaxt=="n"){axis(side=2, labels=FALSE)}
	# mtext(ConsNameMedium[Cons[i]], side=2, line=2)
	# legend(LegPos, NeatConsNames[Cons[i]])
	X <- c("topleft"=0.15, "topright"=ThisAt_Axis[length(ResourceNames)]+1)
	Pos <- c("topleft"=4, "topright"=2)
	text(x=X[LegPos], y=0.95, labels=NeatConsNames[Cons[i]], pos=Pos[LegPos], font=3)
	if(i==length(Cons)){
		mtext("Proportion of Diet", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}
# ============================
# = END new _v0.5.0 (violin) =
# ============================

# ========================
# = Sensitivity Box Plot =
# ========================
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("Sensitivity", Version, sep=""), ".pdf", sep=""), height=6.5, width=6.81)}
	if(SaveType==".png"){png(file=paste(paste("Sensitivity", Version, sep=""), ".png", sep=""), units="in", res=200, height=6.5, width=6.81)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("Sensitivity", Version, sep=""), ".eps", sep=""), height=7, width=6.81)}
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
	tSen <- sensit[[Cons[[i]]]][,NamesThisUse]
	tSen2 <- sensit2[[Cons[[i]]]][,NamesThisUse]
	
	ThisRU[,"Source"] <- factor(ThisRU[,"Source"], levels=ConsChoices[[Cons[i]]][[GroupChoose]], ordered=TRUE)
		
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
		mtext("Proportion of Diet", side=2, line=-1, outer=TRUE, cex=1)
	}
}
if(Save){dev.off()}




# ===============
# = ZOOPLANKTON =
# ===============
# zData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardZoopMass2010&2012.csv")
# cData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/Paul&WardChaobMass2010&2012.csv")
# zData <- reshape(zData0, varying=list(names(zData0[,4:17])), times=names(zData0[,4:17]), ids=1:nrow(zData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
# row.names(zData) <- NULL
# zData <- subset(zData, DoY>=143)
# 
# zData[,"Year"] <- as.factor(zData[,"Year"])
# 
# cData <- reshape(cData0, varying=list(names(cData0[,4:7])), times=names(cData0[,4:7]), ids=1:nrow(cData0), timevar="Taxon", v.names="Mass", direction="long")[,c("Lake","Year","DoY","Taxon","Mass")]
# row.names(cData) <- NULL
# cData <- subset(cData, DoY>=143)
# cData[,"Year"] <- as.factor(cData[,"Year"])
# 
# sumzData <- aggregate(zData[,"Mass"], by=list(zData[,"Lake"], zData[,"Year"], zData[,"DoY"]), sum)
# names(sumzData) <- c("Lake", "Year", "DoY", "Mass")
# zYearMean <- aggregate(sumzData[,"Mass"], by=list(sumzData[,"Lake"], sumzData[,"Year"]), mean)
# names(zYearMean) <- c("Lake", "Year", "Mass")
# sumcData <- aggregate(cData[,"Mass"], by=list(cData[,"Lake"], cData[,"Year"], cData[,"DoY"]), sum)
# names(sumcData) <- c("Lake", "Year", "DoY", "Mass")
# cYearMean <- aggregate(sumcData[,"Mass"], by=list(sumcData[,"Lake"], sumcData[,"Year"]), mean)
# names(cYearMean) <- c("Lake", "Year", "Mass")
# 
# 
# if(Save){
# 	if(SaveType==".pdf"){pdf(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".pdf", sep=""), height=4.75, width=3.23)}
# 	if(SaveType==".png"){png(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".png", sep=""), units="in", res=200, height=4.75, width=3.23)}
# 	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("ZoopChaob_PaulWard_2010&2012_", Version, sep=""), ".eps", sep=""), height=4.75, width=3.23)}
# }else{
# 	dev.new(height=4.75, width=3.23)
# }
# par(mfrow=c(2,1), mar=c(1,2.5,0,0), oma=c(0.3, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
# boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
# axis(side=1, at=c(1,4), labels=FALSE)
# mtext(quote(Zooplankton~(g/m^2)), side=2, line=1.5, cex=1)
# 
# boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
# axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
# mtext(quote(italic(Chaoborus)~spp.~(g/m^2)), side=2, line=1.5, cex=1)
# if(Save){dev.off()}



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

# ==========================================================================
# = New Figure that Combines old Figures 2 & 3 (combined zoop and routine) =
# ==========================================================================
# =======================================
# = Copied from LimnoZoops_CrunchedUp.R =
# =======================================
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

# VarAnalyze <- c("Color","Temp","Zmix","Secchi","Light","DOC","DIC", "pH", "TN", "TP","Chla","POC","PON","C:Chl", "pCO2_water", "DO_Conc")


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


# =============================================
# = 4-panel, areal chlorophyll, pCO2, DO, NEP =
# =============================================

if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("LimnoMetabZoops_PaulWard_2010&2012_", Version, sep=""), ".pdf", sep=""), height=5.5, width=4.33)}
	if(SaveType==".png"){png(file=paste(paste("LimnoMetabZoops_PaulWard_2010&2012_", Version, sep=""), ".png", sep=""), units="in", res=200, height=5.5, width=4.33)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("LimnoMetabZoops_PaulWard_2010&2012_", Version, sep=""), ".eps", sep=""), height=5.5, width=4.33)}
}else{
	dev.new(height=5.5, width=4.33)
}
par(mfrow=c(3,2), mar=c(1,2,0.2,0.25), oma=c(0.2, 0, 0.2, 0.2), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)

#Chlorophyll
boxplot(aChla~Year+Lake, data=aChl0, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(Chlorophyll~(mg/m^2)), side=2, line=1)


#pCO2
par(mar=c(1,2.25,0.2,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
ppCO2 <- subset(Data0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(italic(p)*CO[2]~(mu*atm)), side=2, line=1)

#DO
par(mar=c(1,2,0.2,0.25), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_2010&2012_Metabolism_v0.2.RData")
boxplot(MeanDO~Year+Lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(DO~("%"*saturation)), side=2, line=1)

#NEP
par(mar=c(1,2.25,0.2,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
pWardPaul_Metabolism <- subset(WardPaul_Metabolism, DoY>=143)
pWardPaul_Metabolism[,"Lake"] <- relevel(pWardPaul_Metabolism[,"Lake"], ref="Paul")
boxplot(NEP~Year+Lake, data=pWardPaul_Metabolism, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(NEP~(mmol~O[2]~m^-3~d^-1)), side=2, line=1)

# Zoop biomass
par(mar=c(1,2.25,0,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(Mass~Year+Lake, data=sumzData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(Zooplankton~(g/m^2)), side=2, line=1)

# Chaob biomass
par(mar=c(1,2.25,0,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(Mass~Year+Lake, data=sumcData, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=c("Paul", "Ward"))
mtext(quote(italic(Chaoborus)~spp.~(g/m^2)), side=2, line=1)
if(Save){dev.off()}

# ===============================
# = Taxon-specific zoop biomass =
# ===============================
nZ <- zData[,"Taxon"]=="Nauplii"
zData2 <- zData
zData2[nZ,"Mass"] <- zData2[nZ,"Mass"]/1
if(Save){
	if(SaveType==".pdf"){pdf(file=paste(paste("zoop_byTaxon_2010&2012_", Version, sep=""), ".pdf", sep=""), height=3.22, width=3.22)}
	if(SaveType==".png"){png(file=paste(paste("zoop_byTaxon_2010&2012_", Version, sep=""), ".png", sep=""), units="in", res=200, height=3.22, width=3.22)}
	if(SaveType==".eps"){setEPS();postscript(file=paste(paste("zoop_byTaxon_2010&2012_", Version, sep=""), ".eps", sep=""), height=3.22, width=3.22)}
}else{
	dev.new(height=3.22, width=3.22)
}
par(mfrow=c(2,2), mar=c(1.5,1.5,0.25,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25, oma=c(0,0.5, 0, 0))
tt <- c("Cyclopoid", "Mesocyclops","Calanoid","Nauplii")
for(i in 1:4){
	boxplot(Mass~Year+Lake, data=zData2[zData2[,"Taxon"]==tt[i],], outline=FALSE, col=c(NA, "lightgray"), show.names=FALSE)
	if(i > 2){
		axis(side=1, at=c(1.5,3.5), labels=c("Paul", "Ward"))
	}
	if(i == 1){
		legend("topright",Pos, legend=tt[i], bty="n", inset=c(0, -0.1))
	}else{
		legPos <- ifelse(i==1, "topright", "topleft")
		legend(legPos, legend=tt[i], bty="n", inset=c(-0.2, -0.1))	
	}

}
mtext(bquote(Biomass~(g/m^2)), side=2, line=-0.5, outer=TRUE)
if(Save){dev.off()}

# ===================
# = Isotope biplots =
# ===================
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis")
load(paste("Data+Phyto_", "v0.4.6", ".RData",sep=""))
ThisVersion <- "v0.5.1"
CexScale <- 0.8


CombYears <- which(is.element(Data[,"Type"], c("Macrophyte", "Terrestrial"))) #I don't want to make separate points for samples types that I don't really expect to change between years, such as the macrophytes or the terrestrial samples.  I just want to assume that these were the same in the different years, while the consumers, DOM, POM etc. may have changed.
Data[CombYears,"Year"] <- 2010
ExcludeTaxa <- c("POM","Nija","YWP","PKS","BHD2","Mesocyclops", "Pickerell Weed")
ExcludeSamples <- c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")
EndMemberTaxa <- c("Periphyton","Macrophyte","Terrestrial")
InvertTaxa <- c("Zooplankton","Snail")
BaseData <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, ExcludeSamples) & !is.element(Taxon, ExcludeTaxa) & is.element(Type, c(EndMemberTaxa, InvertTaxa, "Fish"))))#select=c("Year","Trophic","Taxon","d13C","d15N","dD")
TaxID <- as.numeric(BaseData[,"Taxon"])
BaseData <- cbind(BaseData, "TaxID"=TaxID)
XY_ErrorBars <- function(X, Y, sdX, sdY, Plot=FALSE){#This is a function that creates error bars given means and sd's
	#Horizontal Error Bars
	XnaughtH <- X-sdX
	XoneH <- X+sdX
	YnaughtH <- Y
	YoneH <- Y
	if(Plot==TRUE){segments(XnaughtH, YnaughtH, XoneH, YoneH, col="gray")}
	#Vertical Error Bars
	XnaughtV <- X
	XoneV <- X
	YnaughtV <- Y-sdY
	YoneV <- Y+sdY
	if(Plot==TRUE){segments(XnaughtV,YnaughtV,XoneV,YoneV, col="gray")}

	ErrorXlim <- c(min(XnaughtH), max(XoneH))
	ErrorYlim <- c(min(YnaughtV), max(YoneV))
	return(list("ErrorXlim"= ErrorXlim, "ErrorYlim"= ErrorYlim))
}#End function


encircle <- function(ind, cx="d13C", cy=NULL, ...){
	muMat <- as.matrix(mu1[ind,c(cx,cy)])
	sdMat <- as.matrix(sd1[ind,c(cx,cy)])/2
	# f1 <- rbind((muMat-sdMat), (muMat+sdMat))
	f01 <- cbind((muMat[,1]-sdMat[,1]),muMat[,2])
	f02 <- cbind((muMat[,1]+sdMat[,1]),muMat[,2])
	f03 <- cbind(muMat[,1], (muMat[,2]-sdMat[,2]))
	f04 <- cbind(muMat[,1], (muMat[,2]+sdMat[,2]))
	f1 <- rbind(f01, f02, f03, f04)
	# f1 <- muMat
	elo <- ellipsoidhull(f1)
	lines(predict(elo), ...)
	# polygon(predict(elo))
	# f2 <- rbind((muMat*1.1), (muMat-muMat*-0.1))
	# lines(predict(ellipsoidhull(f1)), ...)
	# lines(predict(ellipsoidhull(f2)), col="forestgreen", lty="dashed", lwd=2)
}


MetaPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="MetaPhyto"), "TaxID"=max(unique(TaxID))+2)
EpiPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="EpiPhyto"), "TaxID"=max(unique(TaxID))+1)
Sources_Cons <- rbind(BaseData[,c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")], MetaPhytos, EpiPhytos)
Data <- Sources_Cons
Grouping= c("Taxon", "Year")
IsoNames=c("d13C", "d15N", "dD")
nPlots <- combn(length(IsoNames),2)[,-3]
nPlots <- nPlots[,c(2,1)]
muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), mean)
muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))
sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), sd)
sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
Lims1 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]])
Lims2 <- XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]])


GraphLayers <- c("EndMembers", "Inverts2010&2012", "Fish2010&2012")
for(i in 1:length(GraphLayers)){

	if(GraphLayers[i]=="EndMembers"){
			Data <- droplevels(subset(BaseData, Trophic==0 & is.element(Type, EndMemberTaxa), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	if(GraphLayers[i]=="Inverts2010&2012"){
			Data <- droplevels(subset(BaseData, Trophic>=0 & is.element(Type, InvertTaxa), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	if(GraphLayers[i]=="Fish2010&2012"){
			Data <- droplevels(subset(BaseData, Trophic>=0 & is.element(Type, "Fish"), select=c("Year","Trophic","Taxon","TaxID","d13C","d15N","dD")))
	}
	
	# Data <- droplevels(subset(Data, Trophic>=0 & !is.element(SampleID, c("O-0362", "V-0270", "P-1202", "P-1166", "O-0382", "P-1165", "P-1206", "P-1238", "P-1239", "P-1243", "Z-1110", "Z-1115", "Z-1195", "Z-1170", "O-0405")) & Taxon!="POM" & Taxon!="Nija" & Taxon!="YWP" & Taxon!="PKS" &Taxon!="BHD2" & Taxon!="Mesocyclops", select=c("Year","Trophic","Taxon","d13C","d15N","dD")))



	DataNew <- Data
	UniqueTaxa <- as.character(unique(Data[,"Taxon"]))
	RealTaxa <- as.character(UniqueTaxa) #c("Helisoma trivolvis", "DOM", "Perca flavescens", "Lepomis gibbosus", "Pimaphales promelas", "Periphyton", "Potamogeton pussillus", "Chara sp.", "Nymphaea odorata", "Carex sp.", "Alnus incana subsp. Rugosa", "Skistodiaptomus oregonensis", "Chaoborus spp.", "Nuphar lutea", "Brasenia schreberi", "Trees") #"DOM", 

	if(GraphLayers[i]=="EndMembers"){
		MetaPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="MetaPhyto"), "TaxID"=max(unique(TaxID))+2)
		EpiPhytos <- data.frame("Trophic"=0, subset(EndMembers, Taxon=="EpiPhyto"), "TaxID"=max(unique(TaxID))+1)
		Sources_Cons <- rbind(DataNew, MetaPhytos, EpiPhytos)
	}else{
		Sources_Cons <- DataNew
	}
	
	#BiPlot <- function(Grouping="Taxon", IsoNames=c("d13C", "d15N", "dD"), Data=NULL, SDs=NULL, PlotLabels=TRUE, Legend=TRUE, LegendLoc="topleft", ...){
		#"Grouping" is the name of the column that contains the names of the different varieties of your group
		#"IsoNames" is the name of the columns that contain the isotope values
		#"Data" is the data frame that contains your means (and sd's if you have them)
		#"SD's" is the data frame that contains your standard deviations (see below for more info on making this)
		#"PlotLabels" asks if you want the name of each group variety plotted with it's data point.  If TRUE (the default), names will often overlap
		#"LegendLoc" describes the location of the legend in the event that you don't want the labels to be plotted.  Options include, "topleft"(the default), "bottomleft", "bottomright", "topright", "center", etc. ?legend for more options
		#"..." pass on other arguments to plot(); for example, type 'main="FUN!"' as an argument to the biplot function (i.e., BiPlot(main="FUN!")) to make "FUN!" the title of each plot.

	PubCex <- 1
	Grouping= c("Taxon", "Year")
	IsoNames=c("d13C", "d15N", "dD")
	SDs=NULL
	LegendLoc="topleft"
	Data=Sources_Cons
	PlotLabels=FALSE
	Legend=FALSE

	nPlots <- combn(length(IsoNames),2)[,-3]#Leave this line alone
	# LabelIso <- list(expression(phantom()^13*C), expression(phantom()^15*N), expression(phantom()^2*H))
	GroupNames <- paste("Group",1:length(Grouping),sep=".")#Leave this line alone

	muDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), mean)
	muDataGroup <- subset(muDataGroup, !is.na(get(IsoNames)))

	if(is.null(SDs)){
		sdDataGroup <- aggregate(Data[,c(IsoNames,"Trophic","TaxID")], by=as.list(Data[,Grouping]), sd)
		sdDataGroup <- subset(sdDataGroup, !is.na(get(IsoNames)))
		}else{
			sdDataGroup <- SDs
	}

	nPlots <- nPlots[,c(2,1)]

	LegendTitle <- c("A", "B")
	if(Save){
		
		options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}",
		    "\\usepackage[T1]{fontenc}", "\\usetikzlibrary{calc}",
		    "\\usepackage{amssymb}","\\usepackage{tensor}"))
		tikz(file=paste(FigureFolder,"/BiPlots/Ward2010&2012_IsoBiPlot_", GraphLayers[i], "_", ThisVersion, ".tex", sep=""), width=3.22, height=5.25, standAlone=TRUE, 
		packages = c("\\usepackage{tikz}",
		                 "\\usepackage[active,tightpage,psfixbb]{preview}",
		                 "\\PreviewEnvironment{pgfpicture}",
		                 "\\setlength\\PreviewBorder{0pt}",
		                 "\\usepackage{amssymb}",
						"\\usepackage{tensor}")
		)
	}else{
		dev.new(width=3.5, height=5.25)
	}
	par(family="Times", las=1, mfrow=c(2,1), mar=c(2.1,3.1,0.1,0.1), oma=c(1,0,0,0), cex=PubCex, ps=9, mgp=c(2, 0.5, 0), tcl=-0.4)

	plot(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]], bty="l", xaxt=ifelse(1==1, "n", "s"), cex.axis=PubCex)#
	mtext("$\\delta^2\\mathrm{H}$", side=2, line=2, las=0, cex=PubCex)
	
	# Add Circles
	if(GraphLayers[i]=="EndMembers"){
		mu1 <- muDataGroup #being lazy -- when looping through the groups, I need to use the muDataGroup (and sd) for the end members when drawing circles (called in function)
		sd1 <- sdDataGroup
		floatI <- muDataGroup[,"Taxon"] %in% c("Brasenia schreberi", "Nuphar variegata", "Nymphaea odorata", "Potamogeton nodosus")
		subI <- muDataGroup[,"Taxon"] %in% c("Chara","Najas flexilis", "Potamogeton amplifolius", "Potamogeton pusillus")
		terrI <- muDataGroup[,"Taxon"] %in% c("Alder", "Sedge", "Tamarack", "Tree")
		phytoI <- muDataGroup[,"Taxon"] %in% c("EpiPhyto", "MetaPhyto")
		periI <- muDataGroup[,"Taxon"] %in% c("Periphyton")
	}
	encircle(floatI, cy="dD", col="springgreen4", lty="dashed", lwd=2)
	encircle(subI, cy="dD", col="springgreen4", lty="dashed", lwd=2)
	encircle(terrI, cy="dD", col="burlywood4", lty="solid", lwd=2)
	encircle(phytoI, cy="dD", col="springgreen4", lty="solid", lwd=2)
	encircle(periI, cy="dD", col="springgreen4", lty="solid", lwd=2)
	
	#Add error bars
	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,1]]], muDataGroup[,IsoNames[nPlots[2,1]]], sdDataGroup[,IsoNames[nPlots[1,1]]], sdDataGroup[,IsoNames[nPlots[2,1]]], Plot=T)
	
	# Add points
	# points(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*CexScale, xlab=LabelIso[nPlots[1,1]], ylab=LabelIso[nPlots[2,1]], xlim=Lims1[[1]], ylim=Lims1[[2]])
		points(muDataGroup[,IsoNames[nPlots[1:2,1]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*CexScale, xlab="", ylab="", xlim=Lims1[[1]], ylim=Lims1[[2]])
	TaxID <- as.numeric(muDataGroup[,"TaxID"])
	Move17_Index <- which(muDataGroup[,"TaxID"]==17)
	Move17_d_CND <- c(-0.5,0,8)
	TXTmuDataGroup <- muDataGroup
	if(is.element(15, muDataGroup[,"TaxID"])){
		Move15_Index <- which(TXTmuDataGroup[,"TaxID"]==15)
		Move15_d_CND <- c(0,0,10) #c(0,0,5)
		TXTmuDataGroup[Move15_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move15_Index,c("d13C","d15N","dD")] + Move15_d_CND
		
		Move12_Index <- which(TXTmuDataGroup[,"TaxID"]==12)
		Move12_d_CND <- c(0,0,-0.0)
		TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] <- TXTmuDataGroup[Move12_Index,c("d13C","d15N","dD")] + Move12_d_CND
	}
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,1]]], labels=as.character(TaxID), cex=PubCex*CexScale)
	title(main=LegendTitle[1], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)


	# Create 2nd plot
	plot(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch="", cex=1.2, xlab="", ylab="", xlim=Lims2[[1]], ylim=Lims2[[2]], bty="l", xaxt=ifelse(2==1, "n", "s"), cex.axis=PubCex)#
	mtext("$\\delta^{15}\\mathrm{N}$", side=2, line=2, las=0, cex=PubCex)
	if(2==2){mtext("$\\delta^{13}\\mathrm{C}$", side=1, line=2, outer=FALSE, cex=PubCex)}
	
	# Add circles
	encircle(floatI, cy="d15N", col="springgreen4", lty="dashed", lwd=2)
	encircle(subI, cy="d15N", col="springgreen4", lty="dashed", lwd=2)
	encircle(terrI, cy="d15N", col="burlywood4", lty="solid", lwd=2)
	encircle(phytoI, cy="d15N", col="springgreen4", lty="solid", lwd=2)
	encircle(periI, cy="d15N", col="springgreen4", lty="solid", lwd=2)
	
	# Add error bars
	XY_ErrorBars(muDataGroup[,IsoNames[nPlots[1,2]]], muDataGroup[,IsoNames[nPlots[2,2]]], sdDataGroup[,IsoNames[nPlots[1,2]]], sdDataGroup[,IsoNames[nPlots[2,2]]], Plot=T)

	# Add points
	points(muDataGroup[,IsoNames[nPlots[1:2,2]]], pch=ifelse(muDataGroup[,"Trophic"]==0, 21, 22), bg=ifelse(muDataGroup[,"Year"]==2010, "white", "gray85"), cex=2.7*CexScale, xlab="", ylab="", xlim=Lims2[[1]], ylim=Lims2[[2]])
	text(TXTmuDataGroup[,IsoNames[nPlots[1:2,2]]], labels=as.character(TaxID), cex=PubCex*CexScale)
	#}

	title(main=LegendTitle[2], adj=0.025, line=-1.5, cex.main=PubCex*1.25, font.main=1)
	
	

	# if(GraphLayers[i]=="EndMembers"){		
	# 	encircle(floatI, cy="d15N", col="forestgreen", lty="dashed", lwd=2)
	# 	encircle(subI, cy="d15N", col="forestgreen", lty="dashed", lwd=2)
	# 	encircle(terrI, cy="d15N", col="brown", lty="solid", lwd=2)
	# 	encircle(phytoI, cy="d15N", col="forestgreen", lty="solid", lwd=2)
	# 	encircle(periI, cy="d15N", col="forestgreen", lty="solid", lwd=2)
	# }
	if(Save){dev.off()}
	print(cbind(muDataGroup, "ID"=TaxID)[order(TaxID),])
}
# 
# # ================
# # = Changing Use =
# # ================
# 
# # ===============
# # = New _v0.4.5 =
# # ===============
# setwd(paste("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/",FigureFolder,sep=""))
# 
# tm25 <- rowSums(cbind(Chosen25th[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)
# tmMed <- rowSums(cbind(ChosenMedDiffs[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)
# tm75 <- rowSums(cbind(Chosen75th[,c("All.Terrestrial", "Floating.Macrophytes")]), na.rm=TRUE)
# 
# if(Save){
# 	if(SaveType==".pdf"){pdf(file=paste(paste("dAlgae_vs_terrMac_", Version, sep=""), ".pdf", sep=""), height=3, width=3.23)}
# 	if(SaveType==".png"){png(file=paste(paste("dAlgae_vs_terrMac_", Version, sep=""), ".png", sep=""), units="in", res=300, height=3, width=3.23)}
# 	if(SaveType==".eps"){setEPS(); postscript(file=paste(paste("dAlgae_vs_terrMac_", Version, sep=""), ".eps", sep=""), height=3, width=3.23, pointsize=10)}
# }else{
# 	dev.new(height=4, width=3, pointsize=10, family="Times")
# }
# par(mfrow=c(1,1), mar=c(2.1,2,0.1,0.1), oma=c(0,0,0,0), ps=9, cex=1, mgp=c(3,0.3,0), tcl=-0.25, family="serif")
# 
# aYlim <- c(min(Chosen25th[,"Algae"], na.rm=TRUE), max(Chosen75th[,"Algae"], na.rm=TRUE))*c(1.15,1.05)
# tYlim <- c(min(tm25, na.rm=TRUE), max(tm75, na.rm=TRUE))
# 
# plot(ChosenMedDiffs[,"Algae"], tmMed, xlim=aYlim, xlab="", xaxt="s", ylab="", ylim=tYlim, pch=NA)
# abline(h=0, lty="dashed", col="gray")
# abline(v=0, lty="dashed", col="gray")
# 
# arrows(x0=Chosen25th[,"Algae"], y0=tmMed, x1=Chosen75th[,"Algae"], y1=tmMed, angle=90, code=0, length=0.05, col="gray")
# arrows(x0=ChosenMedDiffs[,"Algae"], y0=tm25, x1=ChosenMedDiffs[,"Algae"], y1=tm75, angle=90, code=0, length=0.05, col="gray")
# 
# points(ChosenMedDiffs[,"Algae"], tmMed, pch=20)
# 
# ConsNames <- c("U. limi", "S. oregonensis", "Phoxinus spp.", "A. melas", "P. promelas", "Chaoborus spp.", "H. trivolvis")
# text(x=ChosenMedDiffs[,"Algae"][-c(1,3)], y=tmMed[-c(1,3)], labels=ConsNames[-c(1,3)], srt=0, xpd=NA,adj=c(1.05,-0.1), cex=1, font=3, ps=9)
# text(x=ChosenMedDiffs[,"Algae"][c(1,3)], y=tmMed[c(1,3)], labels=ConsNames[c(1,3)], srt=0, xpd=NA,adj=c(-0.01,-0.2), cex=1, font=3, ps=9)
# mtext("Change in terrestrial + macrophyte support", side=2, line=1.1, cex=1)
# mtext("Change in algal support", side=1, line=1.1, cex=1)
# if(Save){dev.off()}




