#Chunks
#Ryan D. Batt 20-Sept-2012
#Identify cotinuous chunks of a time series, and identify them with a factor (an ID# for each chunk).  Within each chunk, identify sub-chunks of a certain length (say, every 288*3 observations), and Identify those as separate chunks.  Chunks can be created by continuity (Cont), by number of observations (Sub=), or by both.  Only one factor of ID's will be returned, which will reflect your desired chunking method.
Chunks <- function(x, Cont=TRUE, Sub=NULL){
	Mode <- function(x) {
		ux <- unique(x)
		ux[which.max(tabulate(match(x, ux)))]
	}
	#Calculate deployment ID counter based on changes in the continuity of x (ID goes up 1 every time there is a discontinuity in x)
	dx <- diff(x)
	dx_hat <- round(Mode(dx), 10)
	big_dx <- which(abs(round(dx,10)) > dx_hat*1.50) #give it a --10%-- 50% allowance
	n_chunk <- length(big_dx) + 1
	#3-Dec-2012: The way I'm defining IDs_1 below might be wrong.  I changed "big_dx" to "(big_dx+1)"
	IDs_1 <- approx(x=c(1, (big_dx+1)), y=1:n_chunk, xout=seq_along(x), method="constant", rule=2)$y
	#If changes in the ID counter should be based soley on continuity of x...
	if(Cont & is.null(Sub)){
		return(IDs_1)
	}
	
	#Calculate changes in ID counter simply based on the number of observations (when n obs. > Sub) since the last change in ID
	nout <- length(IDs_1)
	IDs_2 <- rep(1:ceiling(nout/Sub), each=Sub, length.out=nout)
	#If nobs is only criterion for changing the ID counter...
	if(!Cont & !is.null(Sub)){	
		return(IDs_2)
	}
	
	#If the deployment ID counter should change when x is not continuous, OR when there have been "Sub" observations since the last change in ID...
	if(Cont & !is.null(Sub)){
		# big_dx3 <- which(abs(round(dx,10)) > dx_hat | diff(c(IDs_2[1],IDs_2[-length(IDs_2)])) != 0) #everything after the "|" is janky b/c for the first part of the logical statement, an index is returned that indicates where the last of a given ID should be placed in the vector, whereas the second logical statement (the statement after the "|") indicates where the new ID should begin
		d2 <- diff(IDs_2) 
		# w1 <- (d2 != 0)
		# w2 <-  which(abs(round(dx,10)) > dx_hat)+1
		# big_dx3 <- c(w1, w2)
		# which(FALSE) #
		w1 <- which(abs(round(dx,10)) > dx_hat*1.5) +1
		w2 <- which(d2 != 0)+1
		big_dx3 <- sort(union(w1, w2))
		
		# big_dx3 <- which(c1[-length(c1)] > dx_hat | d2 != 0) #everything after the "|" is janky b/c for the first part of the logical statement, an index is returned that indicates where the last of a given ID should be placed in the vector, whereas the second logical statement (the statement after the "|") indicates where the new ID should begin
		n_chunk3 <- length(big_dx3) + 1
		IDs_3 <- approx(x=c(1, big_dx3), y=1:n_chunk3, xout=seq_along(x), method="constant", rule=2)$y
		return(IDs_3)
	}
}
