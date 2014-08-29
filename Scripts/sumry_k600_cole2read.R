
# =========================================
# = Regressions between k.cole and k.read =
# =========================================
load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/orgWard10_sensor.RData")
summary(ward10.lm.read.cole) # R2 = 0.9094
# sum(miss.read2_ward2010)

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/orgWard12_sensor.RData")
summary(ward12.lm.read.cole) # R2 = 0.9016
# sum(miss.read2_ward2012)

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/orgPaul10_sensor.RData")
summary(paul10.lm.read.cole) # R2 = 0.904
# sum(miss.read2_paul2010)

load("/Users/Battrd/Documents/School&Work/WiscResearch/Isotopes_2012Analysis/Data/orgPaul12_sensor.RData")
summary(paul12.lm.read.cole) # R2 = 0.888
# sum(miss.read2_paul2012)