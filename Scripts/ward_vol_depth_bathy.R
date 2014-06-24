#Ryan Batt! 8-Oct-09
#Calculating the area and volume between each contour on Ward Lake Bathymetric Map.  
#The volumes are calculated in terms of a)"volume that is <= 2 meters" b)"volume that is between 3 meters and 5 meters"
#Taking these known areas/volumes to transform volumetric metabolism data into absolute metabolism data for the 0-2m layer, and the 3-5m layer.


#***Other Notes****
#0-2m is the epilimnion, and is referred to as "shallow"
#3-5m is a portion of the metalimnion, and is referred to as "deep"
#Metabolism estimates are from simultaneously deployed sondes at 4m and 1m; deployment length was 12 days.
#1m sonde sampled temp/DO every 5 minutes
#4m sonde sampled temp/DO ever 10 minutes (should have been 5)
#There was a strong oxygen and chl-a maximum at 4m when the sonde data was collected

# rm.list=ls()
# graphics.off()

#Convert pixels to hectares
#List of pixels between each set of contours
PixPart <- c(66, 6598, 17408, 24901, 87470, 28071, 24258, 17435, 24731, 44946)



#Total area of Ward Lake in Pixels
PixTot <- sum(PixPart)

#Total area of Ward Lake in hectares
HectTot <- 2.74



#List of areas between each set of contours in hectares
HectPart <- rep(0, 10)
for (i in 1:10) {
	HectPart[i] <- (PixPart[i]*HectTot)/PixTot
	}
#Set sequence to outer-inner
HP <- rev(HectPart)

#Convert area to square meters (SMP = square meter part)
SMP <- HP * 10000 

ContsFt <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 27)
ContsMet <- ContsFt*.305

#Shorten object names (aids in laziness)
HT <- HectTot
CM <- ContsMet

depth.area <- rev(cumsum(PixPart/PixTot*HectTot)*1000) # the area of the lake that is at least this deep (in meters^2)
layer.info <- data.frame(depth1=CM[1:9], depth2=CM[2:10], h=diff(CM[1:10]), R=embed(depth.area,2)[,2], r=embed(depth.area,2)[,1])
layer.vol <- function(R, r, h, ...){
	pi*(R^2+r*R+r^2)*h/3 # formula for volume of the frustrum of a cone
}
layer.info[,"volume"] <- apply(layer.info, 1, function(x){do.call(layer.vol, as.list(x))})

da.hiRes <- approx(x=CM[1:10], y=rev(cumsum(PixPart/PixTot*HectTot)*1000), xout=seq(CM[1], CM[10], length.out=100))
li.hiRes <- data.frame(depth1=head(da.hiRes$x, -1), depth2=tail(da.hiRes$x, -1), h=diff(da.hiRes$x), R=embed(da.hiRes$y,2)[,2], r=embed(da.hiRes$y,2)[,1])
li.hiRes[,"volume"] <- apply(li.hiRes, 1, function(x){do.call(layer.vol, as.list(x))})

mean.depth.ward <- sum(layer.info[,"depth1"]*(layer.info[,"R"]/sum(layer.info[,"R"])))
