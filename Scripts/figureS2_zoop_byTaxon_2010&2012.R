

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/routine.smry.RData")

Save <- TRUE
SaveType <- c(".png", ".tiff")[2]

# ===============================
# = Taxon-specific zoop biomass =
# ===============================
nZ <- zData[,"Taxon"]=="Nauplii"
zData2 <- zData
zData2[nZ,"Mass"] <- zData2[nZ,"Mass"]/1
if(Save){
	if(SaveType==".pdf"){pdf(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS2_zoop_byTaxon_2010&2012.pdf", height=2.9, width=2.9)}
	if(SaveType==".png"){png(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS2_zoop_byTaxon_2010&2012.png", units="in", res=200, height=2.9, width=2.9)}
	if(SaveType==".eps"){setEPS();postscript(file="/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS2_zoop_byTaxon_2010&2012.eps", height=2.9, width=2.9)}
	if(SaveType==".tiff"){tiff("/Users/battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Figures/Supplement/figS2_zoop_byTaxon_2010&2012.tiff", width=2.9, height=2.9, units="in", res=300, compression="lzw", type="quartz")}
}else{
	dev.new(height=2.9, width=2.9)
}
par(mfrow=c(2,2), mar=c(1.5,1.5,0.25,0), cex=1, ps=8, family="serif", mgp=c(3,0.3,0), tcl=-0.25, oma=c(0,0.5, 0, 0))
tt <- c("Cyclopoid", "Mesocyclops","Calanoid","Nauplii")
for(i in 1:4){
	boxplot(Mass~Year+Lake, data=zData2[zData2[,"Taxon"]==tt[i],], outline=FALSE, col=c(NA,NA), border=c("red", "blue"), show.names=FALSE)
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