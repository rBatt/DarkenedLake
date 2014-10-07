rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/Data/IsotopeData2012")
library("datamerge")
NewCN00 <- read.csv("WardIsotopes_CN_recd07Dec2012.csv")
NewCN <- NewCN00[order(NewCN00[,"SampleID"]),]

Data000 <- read.csv("WardIsotopes_2010&2012_3Jan2013.csv")
Data00 <- Data000[,c("SampleID", "d15N", "d13C")]
Data0 <- Data00[order(Data00[,"SampleID"]),]

MissNames <- as.character(Data0[which(is.na(Data0),  arr.ind=TRUE)[,1], "SampleID"])

FillableNames <- NewCN[,"SampleID"] #as.character(MissNames[is.element(MissNames, NewCN[,"SampleID"])]) #I don't need the long version unless I assume that some of the data in NewCN is already contained in Data0, and that I don't want to replace Data0 data with NewCN data.

Data1 <- Data0
dim(NewCN)
dim(Data1[is.element(Data0[,"SampleID"], FillableNames), 2:3])
which(!is.element(Data1[is.element(Data0[,"SampleID"], FillableNames), 1], NewCN[,"SampleID"]))
which(duplicated(Data1[is.element(Data0[,"SampleID"], FillableNames), 1]))

Data1[is.element(Data0[,"SampleID"], FillableNames), 2:3] <- as.matrix(NewCN[,c("d15N", "d13C")])

Data_NoCN <- Data000[,-which(is.element(names(Data000),c("d15N", "d13C")))]
Data2 <- merge(Data_NoCN, Data1, by="SampleID")
Data3 <- Data2[, names(Data000)]

write.csv(Data3, file="WardIsotopes_2010&2012_09Jan2013.csv")