

DataFull <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataFull.txt", sep="\t", header=TRUE)




ddply(DataFull[DataFull[,"year4"]%in%c(2010L,2012L)&DataFull[,"variable"]=="snow_cm"& DataFull[,"daynum"]<122,], c("location","year4"), function(x)data.frame('snow'=sum(x[,"Data"], na.rm=TRUE)))


Ice_Nrtn000 <- read.csv("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_rawData/Ice_Nrtn.csv")
Ice_Nrtn00 <- Ice_Nrtn000[which(!is.na(Ice_Nrtn000[,"firstopen"]) & !is.na(Ice_Nrtn000[,"lastopen"])),]
Ice_Nrtn00[Ice_Nrtn00[,"year"]%in%c(2010L, 2012L),]