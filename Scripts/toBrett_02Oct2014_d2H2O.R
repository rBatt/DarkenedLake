
library(plyr)
DataRaw <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/IsotopeData2012/WardIsotopes_2010&2012_17Jan2013.csv", header=TRUE)
ddply(DataRaw, c("Year", "Habitat"), function(x){data.frame(muWater=mean(x[x[,"Type"]=="Water","dD"]))})