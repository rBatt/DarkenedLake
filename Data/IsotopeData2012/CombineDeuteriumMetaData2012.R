rm(list=ls())
graphics.off()

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data/IsotopeData2012")
Meta <- read.csv("IsotopeMetaData.csv")
Data0 <- read.csv("IsotopeNewData2012_dD.csv")

Data <- merge (Data0, Meta, by="SampleID", all.x=TRUE)

write.csv(Data, file="IsotopeData_dD&Meta_2012.csv")

#the written CSV was then manually put in the format to match WardIsotopes2010_CompJan2012, and then I added a few columns (Length, Weight, GutID, FishID...), and saved as WardIsotopes_2010&2012_4Dec2012.csv