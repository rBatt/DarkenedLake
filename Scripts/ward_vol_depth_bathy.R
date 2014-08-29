#Ryan Batt! 08-Oct-2009
#Calculating the area and volume between each contour on Ward Lake Bathymetric Map.  
#The volumes are calculated in terms of a)"volume that is <= 2 meters" b)"volume that is between 3 meters and 5 meters"
#Taking these known areas/volumes to transform volumetric metabolism data into absolute metabolism data for the 0-2m layer, and the 3-5m layer.

# NOTE: this was updated by RDB in Winter 2013/ Spring 2014


# ==============================================================
# = Calculate Ward area (pixels and ha) per-countour and total =
# ==============================================================

# note that the bathy map used for the contours lists 2.74 ha as the area
# But examining Ward from Google satellite view, 
# and comparing the Google Ward area to the Google Paul area (and known Paul ha)
# indicates that Ward is ~1.9 ha.

# Convert pixels to hectares
# List of pixels between each set of contours
# Pixexs were counted by coloring area between contours different colors in Photoshop, then counting the number of pixels of each color
PixPart <- c(66, 6598, 17408, 24901, 87470, 28071, 24258, 17435, 24731, 44946) # the number of pixels, from the inside to the outside

#Total area of Ward Lake in Pixels
PixTot <- sum(PixPart)

#Total area of Ward Lake in hectares

HectTot <- 2.74
# HT <- HectTot # total hectares


#List of areas between each set of contours in hectares
HectPart <- rep(0, 10)
for (i in 1:10) {
	HectPart[i] <- (PixPart[i]*HectTot)/PixTot
	}
#Set sequence to outer-inner
HP <- rev(HectPart)

#Convert area to square meters (SMP = square meter part)
SMP <- HP * 10000 

# =============================
# = Depth contours (ft and m) =
# =============================
ContsFt <- c(0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 27) # depth of the contours on the bathy map
ContsMet <- ContsFt*.305 # convert feet to meters


# =============================
# = Calculate volume and area =
# =============================
layer.vol <- function(R, r, h, ...){
	pi*(R^2+r*R+r^2)*h/3 # formula for volume of the frustrum of a cone
}

depth.area <- rev(cumsum(PixPart/PixTot*HectTot)*1000) # the area of the lake that is at least this deep (in meters^2)

# note that the column names correspond to the formula for the frustrum of a cone
layer.info <- data.frame(depth1=ContsMet[1:9], depth2=ContsMet[2:10], h=diff(ContsMet[1:10]), R=embed(depth.area,2)[,2], r=embed(depth.area,2)[,1])

layer.info[,"volume"] <- apply(layer.info, 1, function(x){do.call(layer.vol, as.list(x))}) # the volume of each depth interval

# =========================================================
# = Calculate volume and area at finer spatial resolution =
# =========================================================
# Calculate depth/ area (da) at a fine spatial resolution (hiRes) via linear interpolation
da.hiRes <- approx(x=ContsMet[1:10], y=rev(cumsum(PixPart/PixTot*HectTot)*1000), xout=seq(ContsMet[1], ContsMet[10], length.out=100))

# Calculate "layer info" (li) at a fine spatial resolution (hiRes) using interpolated depth/ area info
li.hiRes <- data.frame(depth1=head(da.hiRes$x, -1), depth2=tail(da.hiRes$x, -1), h=diff(da.hiRes$x), R=embed(da.hiRes$y,2)[,2], r=embed(da.hiRes$y,2)[,1])
li.hiRes[,"volume"] <- apply(li.hiRes, 1, function(x){do.call(layer.vol, as.list(x))}) # the volume of each depth interval, but now at a finer spatial resolution


# ======================================
# = Calculate Ward mean depth (coarse) =
# ======================================
mean.depth.ward <- sum(layer.info[,"depth1"]*(layer.info[,"R"]/sum(layer.info[,"R"])))





