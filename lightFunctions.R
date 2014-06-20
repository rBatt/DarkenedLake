# =================
# = Find 1% light =
# =================
find1 <- function(x){
	which.min(abs(x[x>1E-9]-0.01))
}

# ===========================
# = Calculate average light =
# ===========================
irrInt <- function(I0, z1, z2, kd){
	# I0 = PAR at surface (above water)
	# z1 = shallow depth
	# z2 = deep depth
	# kd = light attenuation cofficient (slope of lm(light%~log(z)) when z ranges between z1 and z2)
	(-exp(-kd*z2)*I0)/kd - (-exp(-kd*z1)*I0)/kd # returns the summation of PAR between z1 and z2
}

# =========================================================
# = Calculate average light in top, sonde, and bot layers =
# =========================================================
doLight <- function(x){
	profCols <- grepl("irr_[0-9]",names(x))
	profDepths <- as.numeric(gsub("irr_", "", names(x)[profCols]))
	profExists <- apply(x, 1, function(x)any(!is.na(x[profCols])))
	
	# Interpolate the light profile
	x[,"irr_0"] <- 1 # all surface is 100%
	prof.irr <- x[profExists,profCols]
	prof.times <- x[profExists,"datetime"]
	for(i in 2:sum(profCols)){
		if(all(is.na(prof.irr[,i]))){
			x[,sum(!profCols)+i] <- 0 # this line could be bad if a depth was skipped
		}else{
			x[,sum(!profCols)+i] <- approx(x=prof.times, y=prof.irr[,i], xout=x[,"datetime"], method="constant", rule=2, f=0)$y
		}	
	}
	prof.irr <- x[,profCols]
	
	# Calculate 1% light
	x[,"z.1perc"] <- profDepths[apply(x[,profCols], 1, find1)]
	
	shalEpi <- x[,"z.mix"] < 0.7
	
	# loop through to calculate kd's and average irr's
	# kd.top <- rep(NA, nrow(x))
	kd.sonde <- rep(NA, nrow(x))
	kd.bot <- rep(NA, nrow(x))
	
	# irr.top <- rep(NA, nrow(x))
	irr.sonde <- rep(NA, nrow(x))
	irr.bot <- rep(NA, nrow(x))
	
	# size.top <- rep(NA, nrow(x))
	size.sonde <- rep(NA, nrow(x))
	size.bot <- rep(NA, nrow(x))
	
	for(i in 1:nrow(x)){
		if(is.na(shalEpi[i])){
			next
		}
		
		t.zmix <- x[i,"z.mix"]
		t.epi.ind <- which.min(abs(profDepths-t.zmix))
		t.irr0 <- as.numeric(prof.irr[i,])
		t.irr.L <- t.irr0>0 & is.finite(t.irr0)
		t.irr <- t.irr0[t.irr.L]
		
		t.bot.ind <- t.irr.L & (cumsum(t.irr.L)>t.epi.ind)
		t.bot.ind.shal <- t.irr.L & (cumsum(t.irr.L)>=5) # index to use for irradiance values that are below 1.0m, but non-0 light
		
		
		kd.all <- log(t.irr)/I(-profDepths[t.irr.L]) # calculate kd as the average attenuation between the surface and each depth
		
		t.irr.surf <- x[i,"irr"]
		
		# kd between 0m and z.mix, or between 0 and 0.5 (if the sonde is below z.mix)
		if(shalEpi[i]){
			# Calculate the thickness of each layer (m), to be used later:
			# E.g., GPP in bot will be GPP.bot <- irr.bot*(GPP/irr.sonde)*size.bot
			# GPP.top <- irr.top*(GPP/irr.sonde)*size.top
			# GPP.sonde <- irr.sonde*(GPP*irr.sonde)*size.sonde == GPP*size.sonde
			# GPP.total <- GPP.top + GPP.bot + GPP.sonde
			# Note that below, when the sonde is in the epi, I am setting size.top = 0, as GPP.top is redundant with GPP.sonde
			# size.top[i] <- 0.5
			size.sonde[i] <- 1.0
			z1.bot <- profDepths[t.bot.ind.shal][1]
			z2.bot <- max(profDepths[t.bot.ind.shal])
			size.bot[i] <- z2.bot-z1.bot
			
			
			# this is if zmix is shallower than the sonde
			
			# calculate kd for each depth in a layer (I0=I@ 0m for each depth), and take average across all depths in layer		
			# kd.top[i] <- sum(kd.all[2:3])/2 # between 0m and 0.5m
			kd.sonde[i] <- sum(kd.all[2:5])/4 # between 0.5 and 1.0
			
			kd.bot[i] <- sum(kd.all[t.bot.ind.shal])/sum(t.bot.ind.shal) # between 1.0 and depth w/ minimum non-zero light
			
			# calculate average irradiance in each layer
			# irr.top[i] <- irrInt(I0=t.irr.surf, z1=0, z2=0.5, kd=kd.top[i])/0.5 # / 0.5-0
			irr.sonde[i] <- irrInt(I0=t.irr.surf, z1=0, z2=1.0, kd=kd.sonde[i]) # /0.5 # / 1.0 - 0.5
			
			if(size.bot[i]>0 & is.finite(size.bot[i])){
				irr.bot[i] <- irrInt(I0=t.irr.surf, z1=z1.bot, z2=z2.bot, kd=kd.bot[i])/(z2.bot-z1.bot)
			}else{
				irr.bot[i] <- 0
			}
			
		
		}else{
			# if zmix is not shallower than the sonde
			
			# Layer thickness
			# size.top[i] <- 0
			size.sonde[i] <- profDepths[t.epi.ind]
			
			if(any(t.bot.ind)){ # if there is any light in the bottom
				z1.bot <- profDepths[t.irr.L][t.bot.ind][1]
				z2.bot <- max(profDepths[t.irr.L][t.bot.ind])
				size.bot[i] <- z2.bot - z1.bot
			}else{ #if there isn't, 
				size.bot[i] <- 0L
			}
			
			
			
			
			# calculate average kd between surface and each depth in a layer
			# kd.top[i] <- sum(kd.all[2:t.epi.ind])/(t.epi.ind-1) # between 0m and zmix
			kd.sonde[i] <- sum(kd.all[2:t.epi.ind])/(t.epi.ind-1) # kd.top[i] # the sonde is also between 0m and zmix
			
			t.bot.ind <- (t.irr.L & (cumsum(t.irr.L)>=t.epi.ind))[t.irr.L] # index to use for irradiance values that are below zmix, but non-0 light
			# Note: subsetting with t.irr.L is necessary, otherwise t.bot.ind doesn't match the length of kd.all, and causes an error when light goes to 0 then back up above 0
			kd.bot[i] <- sum(kd.all[t.bot.ind])/sum(t.bot.ind) # between zmix and depth w/ minimum non-zero light
			
			# calculate average irradiance
			# irr.top[i] <- irrInt(I0=t.irr.surf, z1=0, z2=profDepths[t.epi.ind], kd=kd.top[i])/profDepths[t.epi.ind]
			irr.sonde[i] <- irrInt(I0=t.irr.surf, z1=0, z2=profDepths[t.epi.ind], kd=kd.sonde[i])/profDepths[t.epi.ind] # irr.top[i]

			if(size.bot[i]>0L & is.finite(size.bot[i])){
				irr.bot[i] <- irrInt(I0=t.irr.surf, z1=z1.bot, z2=z2.bot, kd=kd.bot[i])/(z2.bot-z1.bot)
			}else{
				irr.bot[i] <- 0L
			}
			
		} # end if - else
	
	} # end looping through rows
	
	# data.frame(x[,1:5], kd.top=kd.top, kd.sonde=kd.sonde, kd.bot=kd.bot, irr.top=irr.top, irr.sonde=irr.sonde, irr.bot=irr.bot, size.top=size.top, size.sonde=size.sonde, size.bot=size.bot)
	# data.frame(x[,1:5], irr.top=irr.top, irr.sonde=irr.sonde, irr.bot=irr.bot, dz.top=size.top, dz.sonde=size.sonde, dz.bot=size.bot)
	data.frame(x[,1:5], irr.sonde=irr.sonde, irr.bot=irr.bot, dz.sonde=size.sonde, dz.bot=size.bot)
	
}



