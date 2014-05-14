#Ryan Batt
#Created: 07-Jan-10
#Let's make a simple metabolism program in R

#Version 1.0.0, 15-Jan-10
#Version 1.0.1, 16-Apr-2010-- Fixed slow step
#Version 1.2.0, 17-Apr-2010-- Fixed error where NEP estimates where a factor of "TstepMin" lower than they should have been.
#Version 1.3.0, 3-May-2010-- 
#	1) Created compatibility with dataframes instead of csv files (obligatory), 
#	2) correctly implemented "k600=.4" argument (previously the code just used .4, and didn't call k600, thus making it difficult to turn off atmospheric flux), 
#	3) added totals of metabolic parameters (summed days) to the output, 
#	4) renamed and reformatted output//output objects, 
#	5) added appropriate title to figure legend, 
#	6) asked figure to call new graphical device if GraphNew==TRUE (default is TRUE), 
#	7) repositioned legend to "top", 
#	8) changed the y-axis label of the plot to adapt to the value of zmix (zmix can also serve the purpose of unit selection), 
#	9) turned off explanation of how Light() values are organized, 
#	10) added option to graph or not.
#Version 1.3.1, 9-May-2010-- changed subtitles to correctly display the date range (had been 1 day early on both beginning and end)
#Version 1.4.0, 27-Jan-2011-- Calculating k600 from Wind.
#Version 2, 21-Dec-2011-- I'm going to look into why white noise throws off BK estimates.  I'm going to change this function so that it returns and graphs some of the intermediate results (such as the differences in DO, the sums of dDO in the different portions of the day, etc), instead of just the final metabolism estimates.

#Future Additions:
#1. Make use of Deployment ID.  Make it so that metabolism can be calculated for each deployment separately, or as a whole, but realizing that there is a time break.
	#2. Add in calculations for k600.  This should be based on wind, and this can more-or-less just be taken from the Matlab program.
	#3. Add a sensor drift correction.  Linear.
#4. Include gas exchange between water layers.  Bell et al. 2006
#5. Potential to include a built-in noise filtering option.  We'll see how the 535 project goes first ;)
	#6. Option to turn diffusion off..., see #4.

#rm(list=ls())
#graphics.off()

Metabolism <- function(Data,
					FileName="Sonde X",
					GraphNew=TRUE,
					DispGraph=TRUE,
					GraphTitle=NULL,
					zmix=1, #depth of the mixed layer in meters
					latitude=46.28, #the latitude of the lake; default=UNDERC
					atmpress=(.942*760), #atmospheric pressure in mmHg; default=UNDERC
					k600=0.4, #this can alternatively be calculated from wind data. To turn diffusion off, set k600=0.  To specify k600, don't supply wind data.
					Wind=NULL,
					UseLegend=TRUE,
					Default_sub=TRUE,
					Default_xlab=TRUE,
					Default_ylab=TRUE,
					Volumetric=TRUE,
					...){ 


#setwd(Directory)

#***** LET'S LOAD SOME FUNCTIONS!!! *****

#"Light" Function Adapted to R from Matlab (JJColoso) by Ryan Batt, Nov 2009
# ****************************************************************
# Calculate time of sunrise for day of year and #latitude. 
# based on WisconsinLight.m from Jon Cole
# ****************************************************************
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
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Light <- Light(latitude)

# ****************************************************************
#DO Concentration @ 100% Saturation
#Gives you the concentration "µM" of O2 for a given water temperature "celcius" and pressure "mmHg" if the water were 100% saturated with O2.
#Totally stolen from from Coloso/ van de Bogert's Matlab prog.
# ****************************************************************
SatdConc <- function(celcius, mmHg){
# Weiss 1970 Deep Sea Res. 17:721-735
#ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100...
#   + S [B1 + B2 T/100 + B3 (T/100)2]
#where DO is ml/L; convert to mg/L by mult by 1.4276
#where
A1=-173.4292 
A2=249.6339 
A3=143.3483 
A4 = -21.8492
##and 
B1= -0.033096 
B2=0.014259 
B3=-0.001700
#and T = temperature degrees K (C + 273.15) S=salinity (g/kg, o/oo)
T=celcius+273.15
DOt=(exp(((A1 + (A2*100/T) + A3*log(T/100) + A4*(T/100)))))
#pressure correction:
#from USGS memo #81.11 and 81.15 1981, based on empiracal data in Handbook
#of Chemistry and Physics, 1964.
u=10^(8.10765 - (1750.286/(235+celcius)))
DOsat=(DOt*((mmHg-u)/(760-u))) #ml/L
DOsat=DOsat/1000 #L/L
#conversion factor
#convert using standard temperature and pressure.  Similar to calculating
#saturation DO at STP in ml/L, converting to mg?L (at STP), and then doing
#the above temperature and pressure conversions.
R=0.082057 # L atm deg-1 mol-1
O2molwt=15.999*2
convfactor=O2molwt*(1/R)*(1/273.15)*(760/760) #g/L
DOsat=DOsat*convfactor*1000 #mg/L
DOsat=DOsat*1000/(15.999*2) #convert to µM
return(DOsat)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# ****************************************************************
# % Saturation to Concentration (mg/L)
# Converts the % Saturation of oxygen normalized to sea level (this is how the sonde outputs the % sat, which is why after you calibrate %sat is something like 94%... depending on ur elevation) at a given temperature of the water (celcius) into an actual concentration (literally indicating the amount of oxygen present) in units of µM
# ****************************************************************
Sat2Conc <- function(SDO2sat, SDtemp){
#Function originally by MCV
#SDO2sat is "sonde data O2(%saturation)"
#SDtemp is "sonde data temperature", which is in C
#calculate oxygen in mg/L from percent saturation, where percent sat in data is relative to sea level and formula below uses saturation value at sea level to get true O2 concentration in mg/L.  
SDO2mgL=(SDO2sat/100)*(-0.00000002057759*SDtemp^5 + 0.000002672016*SDtemp^4
+ -0.0001884085*SDtemp^3 + 0.009778012*SDtemp^2 + -0.4147241*SDtemp + 14.621)
# change units to uM
SDO2uM=SDO2mgL*1000/(15.999*2)
return(SDO2uM)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Data[,7] <- Sat2Conc(Data[,5], Data[,4]) #Observed DO concentration in µM
Data[,8] <- SatdConc(Data[,4], atmpress) #DO concentration in µM at 100% saturation

#When to start and stop
#I know there is a better way.
Start <- which(subset(Data,DoY==min(Data$DoY),Fract)#Subset 1st day fractions
 				<=Light[min(Data$DoY),1])#that are <= the sunrise fraction for that DoY.  In other words-- set Start equal to the first day in the data if the fraction on the first day is before sunrise--- ha, like that'll ever happen!
 				
if(length(Start)>0){Start<-min(Data$DoY)}else{Start<-min(Data$DoY)+1}#And since it'll never happen (meaning Start before this line would still be "integer(0)"), if the length is still zero, just set the Starting DoY equal to the day after the minimum DoY.  This approach will have to be modified when I am updating this function to hand DepID's.

Stop <- which(subset(Data,DoY==max(Data$DoY),Fract) 				>=Light[max(Data$DoY),2])#Same idea as with start, except now I want to make sure that the Stop day has data past sunset. Again, that'll prolly never happen.
if(length(Stop)>0){Stop<-max(Data$DoY)}else{Stop<-max(Data$DoY)-1}

#How long each time step is (sampling frequency)
Tstep <- Data$Fract[2] - Data$Fract[1]#might be off by 1 100,000th


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# SLOW STEP!!!!! Takes about 5 seconds after it has been run once, and the same for being run the first time.  Pre-allocating memory via pre-determining vector length does not help. 
#This is just subtraction.
#No longer a slow step--- diff function fixed it (16-Apr-2010)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data[,9] <- c(NA,diff(Data[,7]))
#for(m in 2:length(Data[,1])){Data[m,9] <- Data[m,7] - Data[m-1,7]}#obs dDO in units of µM (mmmol/m^3) Stupid Slow Step.  Should try using the diff() function in the future-- it might be faster (I bet it is...lol); O yeah, it is.
	
schmidt <- (1800.6-120.1*Data[,4])+(3.7818*Data[,4]^2)-0.047608*Data[,4]^3

if(k600!=0){
	if(!is.null(Wind)){
		Wind10 <- Wind/(2/10)^0.15 #wind at 10m, from wind at 2m
		k600 <- (2.07+0.215*Wind10^(1.7))/100 *24 #k600 in m/day
		}
	}
TstepMin <- Tstep * 1440 #time step in minutes; 1440 is the number of minutes in a day.
#calculate k for O2 from k600 and Sc from Wanninkhoff 1992
kO2 <-k600*(schmidt/600)^-0.5 #.4 is k600 estimate. kO2 is in units of m/day
kO2 <- kO2/1440 #into units of m/minute
Data[,10] <- kO2*(Data[,8]-Data[,7]) #Change in DO due to gasflux in µM*m/minute (mmol/m^2/minute) #uM*m/min = mmol/m^3*m/min = mmol/m^2/min

Gradient <- Data[,8] - Data[,7]

Data[,11] <- ((Data[,9]/TstepMin)*zmix) - Data[,10] #change in DO due to biology in mmol/m^2/minute #uM/4min*2m = mmol/m^3/min*m = mmol/m^2/min
if(Volumetric==TRUE){
	Data[,11] <- Data[,11] / zmix
	
}

dailyR <- c()
dailyGPP <- c()
dailyNEP <- c()
daytimeNEP <- c()

#*****************************************************************
# Respiration Function
#This function will calculate respiration for a given day or series of days.  Makes use of the 3 previous functions and additions to Data dataframe, so this function is not a completely independent module.
#*****************************************************************
Resp <- function(Day){
	#Respiration is comprised of 4 components: evening of the day before, morning day of, evening day of, morning next day.
	mmRdDO <- c() #Should end up length of 4, containing the average for each component. The "mm" is just stands for "mean of mean".
	Rtoday <- c() #length=length(Day); a vector whose values represent the average respiration for a given day.
	
	#Added 21-Dec-2011
	#Collect the mRdDO's for each day
	df_mRdDO <- data.frame("Day"=Day, "RComp1"=rep(NA, length(Day)), "RComp2"=rep(NA, length(Day)), "RComp3"=rep(NA, length(Day)), "RComp4"=rep(NA, length(Day)))#Added 21-Dec-2011
	mRdDO_duration <- rep(NA,4)#Added 21-Dec-2011
	df_mRdDO_duration <- data.frame("Day"=Day, "RComp1"=rep(NA, length(Day)), "RComp2"=rep(NA, length(Day)), "RComp3"=rep(NA, length(Day)), "RComp4"=rep(NA, length(Day)))#Added 21-Dec-2011
	
	for(l in 1:length(Day)){
		i=Day[l]
		#Respiration is comprised of 4 components: (1) evening of the day before, (2) morning day of, (3) evening day of, (4) morning next day.
		#Determine which components are available for the given day
		#Component1
		if(i==min(Data$DoY)){
			Comp1=FALSE
			#Comp2=FALSE
			}else{
			if(min(subset(Data,DoY==i-1,Fract))<=Light[i-1,2]+(1/24)){
			Comp1=TRUE 
			#Comp2=TRUE
			}else{
			Comp1=FALSE
				#Comp2=FALSE
			}
		} #If Comp1 is true Comp2 is as well. See comment below.
		#Component2
		#I will never put a sonde out past dark, so F this step... the presence of this component can be deduced from the presence of Component 1.
		#if the smallest fraction is less than or equal to 2am, we have the morning-day-of component of night:
		if(min(subset(Data,DoY==i,Fract))<=(2/24)){Comp2=TRUE}else{Comp2=FALSE}
		#Component3
		#Similarly, I'm not picking sondes up at night. F it.
		if(max(subset(Data,DoY==i, Fract))>=22/24){Comp3=TRUE}else{Comp3=FALSE}
		#Component4
		if(i==max(Data$DoY)){
			Comp4=FALSE
			#Comp3=FALSE
			}else{
			if(max(subset(Data,DoY==i+1,Fract))>=Light[i+1,1]-(1/24)){
			Comp4=TRUE 
			#Comp3=TRUE
			}else{Comp4=FALSE
				#Comp3=FALSE
			}
		}
		#So far this function tells me what components are available to use in calculating respiration for this day.  Right now they are a little redundant, and probably unnecessary.  I could've probably just said that all 4 are TRUE unless it's the first day, then only 3&4 are true, or if it's the last day, when only 1&2 are true.

		Comps <- c(Comp1, Comp2, Comp3, Comp4) #The components
		CompsT <- c((i-1), i, i, (i+1)) #Day on which a component is found
		RdDO <- c() #Unaveraged diffrences in DO; vector/object to be reused for each Comp
		mRdDO <- c() #Should end up length of 4, containing the average for each component (unless a day is missing components)
		#Comp1
		if(Comp1==TRUE){ #Is this component available for day i?
			RdDO <- subset(Data,
				 DoY==CompsT[1] & Fract>=(Light[CompsT[1],2]+(1/24)), #When day is equal to what CompsT says it should be for this Comp; when the fraction of the day is indicative of being at or after sunset plus 1 hour...
				 V11) #...give us V11, the ΔO2 due to biology (not diffusion) in µM/minute
				mRdDO[1] <- mean(na.omit(RdDO))
				mRdDO_duration[1] <- length(RdDO[[1]])#Added 21-Dec-2011
		}#The average respiration/minute for this Component
		#Comp2*
		if(Comp2==TRUE){
			RdDO <- subset(Data,DoY==CompsT[2] & Fract<=(Light[CompsT[2],1]-(1/24)), V11)
			mRdDO[2] <- mean(na.omit(RdDO))
			mRdDO_duration[2] <- length(RdDO[[1]])#Added 21-Dec-2011
		}		
		#Comp3*
		if(Comp3==TRUE){
			RdDO <- subset(Data,DoY==CompsT[3] & Fract>=(Light[CompsT[3],2]+(1/24)), V11)
			mRdDO[3] <- mean(na.omit(RdDO))
			mRdDO_duration[3] <- length(RdDO[[1]])#Added 21-Dec-2011
		}	
		#Comp4
		if(Comp4==TRUE){
			RdDO <- subset(Data,DoY==CompsT[4] & Fract<=(Light[CompsT[4],1]-(1/24)), V11)
			mRdDO[4] <- mean(na.omit(RdDO))
			mRdDO_duration[4] <- length(RdDO[[1]])#Added 21-Dec-2011
		}
		#*Note: Comps 2&3 could be combined (only assuming no midnight cutoff) but I did seperate b/c it is more organized/ easier to follow, plus it doesn't really matter!
		df_mRdDO[l,2:5] <- mRdDO#Added 21-Dec-2011
		df_mRdDO_duration[l,2:5] <- mRdDO_duration #Added 21-Dec-2011
		mmRdDO[l] <- mean(na.omit(mRdDO))
		Rtoday[l] <- mmRdDO[l]*1440
	}#End initial for loop (created so the arugment of Resp() could take more than one value)
	return(list(Rtoday, df_mRdDO, df_mRdDO_duration))
}#End Resp function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#*****************************************************************
# NEP Function
#This function calculates NEP over the entire 24 hour period (1st column), as well as just for the daytime (2nd column)
#*****************************************************************
NEP <- function(i){
	NEPs <- matrix(ncol=2, nrow=length(i))
	for(l in 1:length(i)){
		NEPs[l,1] <- TstepMin*sum(na.omit(subset(Data, DoY==i[l], V11)))
		NEPs[l,2] <- TstepMin*sum(subset(Data, DoY==i[l]&Fract>=Light[i[l],1]&Fract<=Light[i[l],2],V11))} #I am SO not using this...
		return(NEPs)}#End NEP function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#*****************************************************************
# GPP Function
#GPP based on extrapolated (to 24hrs) R and 24-hour NEP
#*****************************************************************
GPP <- function(i){
	for(l in 1:length(i)){ #i is the day, or list of days.
		dailyGPP[l] <- NEP(i[l])[,1] - Resp(i[l])[[1]]}#GPP as the balance between 24-hour NEP and 24-hour R.  NEP is a function, whose argument is the day of year.  If the NEP function was indexed as [,2], this would give daytime NEP, as opposed to 24-hr NEP
		return(dailyGPP)}#End GPP function
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Resp, NEP, and GPP are all functions.  Turning functions into numeric-class objects.
tResp <- Resp(Start:Stop)[[1]]#Respiration for all DoY's which have a full day
tNEP <- NEP(Start:Stop)[,1]#Same for NEP
tGPP <- GPP(Start:Stop)#and GPP.
		
#Converting the first day in DoY to first day in yyyy-mm-dd format; used in graph
DateOrigin <- paste(as.character(Data[1,1]),"-01-01",sep="")
Bdate <- as.Date(Start,origin=DateOrigin)
#Bdate <- paste(as.character(Start), Data[1,1], sep="-")#Bdate=beginning date
#Bdate <- as.Date(Bdate, format="%j-%Y")
#Converting the last day in DoY to last day in yyyy-mm-dd format; used in graph
Edate <- as.Date(Stop,origin=DateOrigin)
#Edate <- paste(as.character(Stop), Data[1,1], sep="-")#Edate = ending date
#Edate <- as.Date(Edate, format="%j-%Y")

# xout <- seq(Start, Stop, length.out=length(diffinv(Data[-1,11]*TstepMin)))
# plot(xout, Data[,5], ylab="", type="l")
# mtext("DO %Sat", side=2, line=2, cex=0.66)
# par(new=TRUE)
# plot(xout, Data[,4], col="red", type="l", xaxt="n", yaxt="n", ylab="", xlab="")
# axis(side=4)
# mtext(expression(Temperature~(degree*C)), side=4, line=2, cex=0.66)
# plot(xout, Data[,7], xlab="", ylab="", type="l")
# mtext(expression(mmol~O[2]~m^-3), side=2, line=2, cex=0.66)
# # plot(xout, Data[1,7]+diffinv(Data[-1,11]*TstepMin), xlab="", ylab="", type="l")
# # mtext(expression(mmol~O[2]~m^-3~-~Atm.~Flux), side=2, line=2, cex=0.66)
# plot(xout[-1], (Data[-1,11]*TstepMin), xlab="", ylab="")
# mtext((mmol~O[2]~m^-3~(5~min.)^-1~-~Atm.~Flux), side=2, line=2, cex=0.66)
# legend("topleft", legend=paste("sd=",round(sd((Data[-1,11]*TstepMin)), 2),sep=""), cex=1)
# # barplot(t(as.matrix(Resp(Start:Stop)[[3]])[,-1]), beside=TRUE)
# # barplot(t(as.matrix(Resp(Start:Stop)[[2]][,-1])*5), beside=TRUE, names.arg=Start:Stop)
# barplot(t(as.matrix(Resp(Start:Stop)[[2]][,-1])*1440), beside=TRUE, names.arg=Start:Stop)



if(DispGraph==TRUE){
if(GraphNew==TRUE){dev.new()}
if(Volumetric==1){yLabel=expression(phantom(0) * (mu *M * phantom(0) * O[2] *phantom(0)* day^-1) )}else{yLabel=expression(phantom(0) * (mmol*phantom(0)*O[2] *phantom(0)*m^-2 *phantom(0)* day^-1) )}

if(Default_sub==TRUE){sub2use=paste(Bdate,Edate,sep="  to  ")}else{sub2use= Default_sub}
if(Default_xlab==TRUE){xlab2use="Day of Year"}else{xlab2use= Default_xlab}
if(Default_ylab==TRUE){ylab2use=yLabel}else{ylab2use= Default_ylab}

plot(Start:Stop, tResp, #plot R24
	col="red",
	type= "o",
	pch=19,
	lty=1,
	main=GraphTitle,
	xaxp= c(Start,Stop,length(Start:Stop)-1),#Tick marks; (first,last,number of intervals)
	ylim= c(min(c(tResp, tNEP, tGPP)), max(c(tResp, tNEP, tGPP))),#upper/lower y-axis limits
	sub=if(Default_sub!=FALSE){sub2use}else{NULL},#subtitle (date range... a lot of times it is hard to accurately place DoY in normal date terms, this should be nice:) )
	xlab=if(Default_xlab!=FALSE){xlab2use}else{NULL},
	ylab=if(Default_ylab!=FALSE){ylab2use}else{NULL},
	font.sub=2,
	...
	)
abline(h=0, col="gray")
points(Start:Stop, tNEP, #plot NEP
	col="purple",
	type="o",
	pch=19,
	lty=1
	)
points(Start:Stop, tGPP, #plot GPP
	col="green",
	type="o",
	pch=19,
	lty=1
	)
if(UseLegend==TRUE){
	legend("top", #Legend; labels are metab. param. + the mean numeric result for each
		legend=c(
			paste("GPP",as.character(round(mean(tGPP),2)),sep=" = "), 
			paste("NEP",as.character(round(mean(tNEP),2)),sep=" = "), 
			paste("R",as.character(round(mean(tResp),2)),sep=" = ")),
		text.col=c('green','purple','red'),
		col=c('green','purple','red'),
		ncol=1,
		pch= c(19,19,19),
		title="Deployment Mean",
		title.col="brown",
		bty="n",
		lty=c(1,1,1)
		)
	}	
}#end if loop for graphing (or not)




matrix <- matrix(nrow=length(Start:Stop), ncol=4, data=c(Start:Stop, Resp(Start:Stop)[[1]], GPP(Start:Stop), NEP(Start:Stop)[,1]), dimnames=list(NULL,c("DoY",'R', 'GPP', 'NEP')))

#MeanMatrix <- matrix(nrow=1, ncol=3, data=c(mean(tResp), mean(tGPP), mean(tNEP)), dimnames=list("Means:", c('R', 'GPP', 'NEP')))

#KMat <- matrix(nrow=length(Data[,10]), ncol=2)
#colnames(KMat) <- c("kO2", "Diff")
#KMat[,1] <- TstepMin *kO2
#KMat[,2] <- TstepMin *Data[,10]

#return(list(MeanMatrix,colSums(matrix),matrix, KMat, Gradient))
return(matrix)
}#End Metabolism function ;)