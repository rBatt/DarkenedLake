

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/kf.epi.good.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/AllDO.RData")
# load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/kf.epi.good.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Results/WardPaulMetabolism.kf.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/routine.smry.RData")

Save <- TRUE
SaveType <- ".png"

# =============================================
# = 4-panel, areal chlorophyll, pCO2, DO, NEP =
# =============================================

if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig2_LimnoMetabZoops_PaulWard_2010&2012.pdf", height=5.5, width=4.33)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig2_LimnoMetabZoops_PaulWard_2010&2012.png", units="in", res=600, height=5.5, width=4.33, type="quartz")}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/fig2_LimnoMetabZoops_PaulWard_2010&2012.eps", height=5.5, width=4.33)}
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
ppCO2 <- subset(routineData0, Layer=="PML", select=c("Lake","Year","pCO2_water"))
boxplot(pCO2_water~Year+Lake, data=ppCO2, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(italic(p)*CO[2]~(mu*atm)), side=2, line=1)

#DO
par(mar=c(1,2,0.2,0.25), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
boxplot(MeanDO~year+lake, data=AllDO, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(DO~("%"*saturation)), side=2, line=1)

#NEP
par(mar=c(1,2.25,0.2,0), cex=1, ps=9, family="serif", mgp=c(3,0.3,0), tcl=-0.25)
# pWardPaulMetabolism.kf <- subset(WardPaulMetabolism.kf, doy>=143)
pWardPaulMetabolism.kf <- WardPaulMetabolism.kf[WardPaulMetabolism.kf[,"doy"]>=143 & WardPaulMetabolism.kf[,"GPP"]>0 & WardPaulMetabolism.kf[,"R"]<0,]
pWardPaulMetabolism.kf[,"lake"] <- relevel(pWardPaulMetabolism.kf[,"lake"], ref="Paul")
boxplot(NEP~year+lake, data=pWardPaulMetabolism.kf, at=c(0.5,1.5, 3.5, 4.5), col=c(NA,"lightgray"), show.names=FALSE, outline=FALSE, lwd=1.25)
axis(side=1, at=c(1,4), labels=FALSE)
mtext(quote(NEP~(mg~O[2]~L^-1~d^-1)), side=2, line=1)

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
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/zoop_byTaxon_2010&2012.pdf", height=3.23, width=3.23)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/zoop_byTaxon_2010&2012.png", units="in", res=600, height=3.23, width=3.23)}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/zoop_byTaxon_2010&2012.eps", height=3.23, width=3.23)}
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