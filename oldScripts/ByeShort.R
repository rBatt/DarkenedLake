#ByeShort
#Ryan D. Batt 20-Sept-2012
#Remove trunc(ToCount) time units containing missing data, and any rows where ToCount is duplicated within subsets defined by the combinatins of By.
#This function is intended to be used with several time series contained in a data frame, with one of the columns representing the time index in units+fractions of whole units.  E.g., if time is index as day of year, you might see 135.5 for noon on julian day 135.  Given this setup, the function looks for how many observations are expected to be observed per unit time (e.g., measurements per day), and removes any instances when this expectation is not exactly met.  I.e., if there is missing data, there will be too few observations, and if there is duplicated data, there will be too many observations.  On the latter point, the definition of duplicates may be subjected to factors in the By columns of the data frame... e.g., "duplicate" observations may have been made at time (ToCount=DoY) 135.003472222, but it may be that observations were made at this time in multiple years, or in different places, and therefore aren't actually replicates.  As such, "By" indicates how the unit time may be nested in other dimensions, as indicated by column vectors of the matrix (By) whose elements can be interpretted as factors.  Another example would be if parameter P is being measured over time, and my time unit is DoY.  I may be measureing in the same place in the same year.  But maybe I am measuring P simultaneously with several instruments that I would like to compare.  In this case I may have a column whose elements are factors indicating which instrument made the measurement in that row.  I would set ToCount="DoY", and By=c("Instrument", "DoY") or By=c("DoY", "Instrument") (either should work, but only the former has been tested).

ByeShort <- function(X, Expected=288, ToCount="DoY", TruncToCount=TRUE, By=c("Year","DoY")){#there is almost certainly a better way to do this
	#in addition to checking for missing values, also checks for duplicates (checks duplicates first, removes those, then checks for missing)
	dups <- function(x) x[!duplicated(round(x[,ToCount],9)),]
	# X[,"Truncd"] <- trunc(X[,ToCount])
	X <- ddply(X, setdiff(By, ToCount), dups)
	# NotDups <- !duplicated(round([,ToCount],9))
	# X <- X[NotDups,]
	ByInd <- data.frame(X[,By], "IND"=1:nrow(X))
	which_nrow <- function(x){ c("Size"=nrow(x), "Start"=min(x[,"IND"]), "Stop"=max(x[,"IND"]))}
	Sizes <- ddply(trunc(ByInd), By, which_nrow)
	TooShort <- Sizes[which(Sizes[,"Size"] < Expected), c("Start","Stop")]
	Start2Stop <- function(x) x[1]:x[2] #these last two steps could probably be combined into an is.element() approach that would be simpler
	WaveTo <- unlist(apply(TooShort, MARGIN=1, FUN=Start2Stop), use.names=FALSE)
	if(length(WaveTo)!=0){
		Xr <- X[-WaveTo,]
	}else{
		Xr <- X
	}
	return(Xr)
}
