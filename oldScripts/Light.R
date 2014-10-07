Light <- function(latitude=46.28){
	lat <- latitude/57.297 #lat converted to radians at Paul L. dyke, 57.297 = 180/pi
	day <- (1:365) #day of year. could be read in from file
	rads <- 2*pi*day/365 

	#dec is the declination of the sun f of day  Spencer, J.W. 1971: Fourier
	#series representation of the position of the Sun. Search, 2(5), 172.
	dec <- 0.006918 - 0.399912 * cos(rads) + 0.070257 * sin(rads) - 0.006758 * cos(2*rads) + 0.000907 * sin(2*rads) - 0.00297 * cos(3*rads) + 0.00148 * sin(3*rads)
	x <- (-1*sin(lat)*sin(dec))/(cos(lat)*cos(dec))

	#hours in radians before true noon of sunrise.
	SR <- (pi/2)-atan(x/(sqrt(1-x^2))) 
	SR <- SR*2/.262 #converts SR from radians back to hours of daylight.
	#need to adjust clock time which is in Central Daylight Savings time by 1 h.

	dates <- (1:365)
	sunrise <- (12-(SR*0.5)+1)/24 #Sunrise in hours Central Daylight Savings time
	sunset <- (12+(SR*0.5)+1)/24   #Sunset in CDST
	RiseSet <- matrix(c(sunrise, sunset), nrow=365, ncol=2)
	#print("Data stored in 'RiseSet' matrix; 1st col is sunRise, 2nd col is sunSet -- given in fractions of day.")
	return(RiseSet)
}