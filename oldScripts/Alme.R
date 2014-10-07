#Alme
#Ryan D. Batt 20-Sept-2012
#Alme = ALign and MErge to data frames containing time series over time unit Al, where Al may be nested withing factors setdiff(Me,Al).
#Write function to use the intersection of a shared column vector of 2 data frames, index the rows of those data frames that belong to the intersections, and then combine the two data frames without duplicating any column names.
Alme <- function(X, Y, Al="DoY", Me=c("Year", "DoY"), CheckMissDups=FALSE, ex=NULL, ey=NULL, Freq_High_Low="High", AllX=FALSE, AllY=FALSE){#This is really only meant for the default for Al at the moment, haven't been thinking outside this box.
	#Freq_High_Low can now be set to "Low" (28-Nov-2012)
	Mode <- function(x) {
		ux <- unique(x)
		ux[which.max(tabulate(match(x, ux)))]
	}
	#3 requirements related to "Me": elements of "Me" need to be the name of columns in both X and Y, be in the same order in both data frames as they are in "Me", and these "Me" columns also need to be the left-most columns of the data frames.
	#first both sets of data should be cleared of any incomplete days
	if(is.null(ex)){
	ex <- round(Mode(1/diff(X[,Al]))) #trunc(1/min(abs(diff(X[,Al])))) #expected frequency of X
	}
	if(is.null(ey)){
	ey <-  round(Mode(1/diff(Y[,Al]))) #trunc(1/min(abs(diff(Y[,Al])))) #expected frequency of Y
	}
	if(CheckMissDups){ #checking for missing data and duplicated data
		X <- as.data.frame(ByeShort(X, Expected=ex)) 
		Y <- as.data.frame(ByeShort(Y, Expected=ey))
	}

	#For each data set to be merged, remove any rows that do not belong to complete days.  Next, identify the highest (or lowest) sampling frequency between the two (if there is a tie, no problem, just use this shared frequency). Next, we need to generate a vector of ideal/ reference sampling times for each data set.
	if(Freq_High_Low=="High") {Freq <- max(ex, ey)}
	if(Freq_High_Low=="Low") {Freq <- min(ex, ey)}
		#For each data set, do the following:
		#1) extract a vector of unique days (w/o fraction of day) ... and make that unique year-day combos, returning the days
	uAlMe <- function(x) unique(x[,Al])
	udaysX <- ddply(trunc(X), Me, uAlMe)[,3] #unique(trunc(X[,Al]))
	udaysY <- ddply(trunc(Y), Me, uAlMe)[,3] #unique(trunc(Y[,Al]))
		#2) repeat each element of the vector Freq times, where Freq is the number of samples per day for the highest sampling frequency between the two data frames
	idaysX <- rep(udaysX, each=Freq)
	idaysY <- rep(udaysY, each=Freq)
		#3) identify the daily interval of this daily sampling frequency (1/Freq)
		#4) repeat this fraction (Freq-1) times, and takes its cumulative sum ... e.g., cumsum(rep((1/288),287)) if Freq=288
	iFracs <- c(0,cumsum(rep((1/Freq),(Freq-1))))
		#5) confirm that the vector in 2) is an even multiple of the vector in 4)
	#do this later
		#6) Add the vector from 4) to the vector in 2)
	idealXd <- idaysX+iFracs
	idealYd <- idaysY+iFracs
		#7) is the time vector from the data identical to the reference/ideal time vecotr?  almost certainly not
	#not going to bother
		#8) if 7) is false, approximate each column of the data given the data time and ideal times ( Data[,i...+2] <- approx(x=Data[,"DoY"], y=Data[,i...+2], xout=ideal)$y ).
	Xnames <- colnames(X)
	Ynames <- colnames(Y)
	# idealx <- matrix(data=NA, nrow=nrow(X), ncol=ncol(X), dimnames=list(NULL, Xnames))
	# idealy <- matrix(data=NA, nrow=nrow(Y), ncol=ncol(Y), dimnames=list(NULL, Ynames))
	idealx <- matrix(data=NA, nrow=length(idealXd), ncol=ncol(X), dimnames=list(NULL, Xnames))
	idealy <- matrix(data=NA, nrow=length(idealYd), ncol=ncol(Y), dimnames=list(NULL, Ynames))
	for(ix in (length(Me)+1):ncol(X)){
		idealx[,ix] <- approx(x=X[,Al], y=X[,ix], xout=idealXd)$y
	}
	for(iy in (length(Me)+1):ncol(Y)){
		idealy[,iy] <- approx(x=Y[,Al], y=Y[,iy], xout=idealYd)$y
	}
		#9) if 7) is false, replace the time column of the data with the ideal time generated in 6) (what i called "ideal" in 8))
	idealx[,Al] <- idealXd
	idealy[,Al] <- idealYd
	PutMeIn <- which(!is.element(Me, Al))
	for(PMI in PutMeIn){
		idealx[,PMI] <- approx(x=X[,Al], y=X[,PMI], xout=idealXd)$y
		idealy[,PMI] <- approx(x=Y[,Al], y=Y[,PMI], xout=idealYd)$y
	}
	# idealx[,PutMeIn] <- X[,PutMeIn]
	# idealy[,PutMeIn] <- Y[,PutMeIn]
	#Use merge() to join the two data frames based on their new time vector, and any other column with which time should intersect in order for data to match (e.g., Year, Lake.. in additino to having the same DoY, the data should only be joined if the DoY's are from the same Lake in the same Year.)
	Ideal <- merge(idealx, idealy, by=Me, all.x=AllX, all.y=AllY)
	#return the joined data frame (or matrix? go df to be safe).
	return(Ideal)
}
