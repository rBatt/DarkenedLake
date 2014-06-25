
# ============
# = ROUTINES =
# ============
#Taken from L&W_Routine_2010&2012_v3.r and then revised
routineData0 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Data/PaulWard_Weekly_2010&2012/PaulWard_Weekly_2010&2012.csv")
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
PML_Data0 <- subset(routineData0, Layer=="PML")[,]
PML_Data <- merge(PML_Data0, POC0, all.x=TRUE)
PML_Data[,"C:Chl"] <- PML_Data[,"POC"]/(PML_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "POC","PON","C:Chl", "Temp", "DOC", "DIC", "DO_Conc")
Meta_Data0 <- subset(routineData0, Layer=="Meta")[,]
Meta_Data <- merge(Meta_Data0, POC0, all.x=TRUE)
Meta_Data[,"C:Chl"] <- Meta_Data[,"POC"]/(Meta_Data[,"Chla"]/1000)

VarAnalyze <- c("Color", "Chla", "Temp", "DOC", "DIC", "DO_Conc")
Hypo_Data <- subset(routineData0, Layer=="Hypo")[,]

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


save(zData, cData, sumcData, sumzData, aChl0, routineData0, file="/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/routine.smry.RData")
