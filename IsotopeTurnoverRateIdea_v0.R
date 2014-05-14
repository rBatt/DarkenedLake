#Ryan Batt
#17-Aug-2012
#How to estimate recent isotopic composition (recent, e.g., after Aquashade was added) using logistic growth to estimate the rate of isotopic turnover

# ===================
# = Logistic Growth =
# ===================
#B=Biomass
#t=Time
#r=Growth Rate Constant (d^-1)
#K= Carrying Capacity Constant (maximum biomass)
#dBdt = Change in biomass per change in time
#dBdtdB = 1/dt = proportional change in biomass per change in time

LogGrowth <- function(B, r, k){
	dBdt <- r*B*(1-B/k)
	dBdtdB <- dBdt/B
	return(dBdtdB)
	}

# ==================================
# = Equilibrium Isotopic Signature =
# ==================================
#Reich et al. 2008. Oecologia. 155:651-663. Equation 1.
#dXEq = Equilibrium isotopic composition (after infinite time incorporating new signature)
#dX0 = Initial isotopic signature (before it began incorporating new signature)
#dXt = The new isotopic signature that it is incorporating
#k= Proportional rate of isotopic incorporation (dBdtdB)
#t= Time
EquiIso <- function(dX0, dXt, k, t){
	dXEq <- (dX0*exp(-k*t) - dXt) / (exp(-k*t) - 1)
	return(dXEq)
}




#How to get "r" from weight of fish of a given species
#Step 1: For each sampling date:
	#Step 1a: Create an empirical frequency distribution of mass
	#Step 1b: Define the number of modes, and their locations (masses)
	#Step 1c: Rank the modes in ascending order (1, 2, 3, ..., N) from mode w/ biggest mass (Rank=1) to smallest mass (Rank=N)
		#Modes appear as rank N
		#After a rank N node appears, another node of rank N will never appear again, even if a node first disappears (i.e., rank identifies a cohort)
		#Modes can go to a "lower" rank, until they arrive at rank 1
		#The number of modes can decrease by 1 when a mode arrives at rank 1
			#If the number of modes decreases under any other condition, the method will not work
			#A mode is said to arrive at rank 1 if the frequency of rank 1 mode increases relative to a standard
				#Standard could be total catch (among species) on that day
				#Say rank 1 is F=100, rank 2 is F=300, and total catch that day is 1500.
				#If next time there is 1 mode, and F1=100...
					#If total catch =750, then F1r0=100/1500, F2r0=300/1500, F1rt=100/750, then the first mode would have a relative increase, so this is allowed
					#If same scenario except the second time the total catch remained 1500 (or became higher), then the total catch would not increase, which means the rank 2 mode died
					#If the realtive (to total interspecies catch) number of individuals missing from the rank 2 mode is much greater than the increase in individuals in the rank 1 mode, this also doesn't work (but I don't know what statistical test to use yet)
		#Modes can maintain rank, and maintain location
		#Modes can maintain rank, but increase in location (mass)
		#Modes cannot decrease in location
		#Modes cannot decrease in rank
#Step 2: Create a matrix for each rank over time 
	#Step 2a: 
	#(time series); start a new vector if a rank disappears, but reappears later
#Step 3: Take the first difference of each rank vector (e.g., Rank2[t,c]-Rank2[t-1,c]..., where "t" is the time index, and "c" is the cohort index (the "c"th time that Rank has 'appeared' in the time series))
		
		
		
		
		