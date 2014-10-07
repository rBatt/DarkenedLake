AllCombinations <- t(combn(c("Terrestrial", "Epi. Phytoplankton", "Meta. Phytoplankton", "All Phytoplankton", "Submersed Macrophytes", "Floating Macrophytes", "All Macrophytes", "Periphyton", NA, NA), m=4))

Checks <- function(x){
	
	c1 <- is.element("Terrestrial", x) #YES
	c2 <- TRUE #any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton", "All Phytoplankton"), x)) #YES
	c3 <- any(is.element(c("Epi. Phytoplankton", "Meta. Phytoplankton"), x)) & is.element("All Phytoplankton", x) #NO
	c4 <- any(is.element(c("Floating Macrophytes", "Submersed Macrophytes"), x)) & is.element("All Macrophytes", x) #NO
	c5 <- is.element("Epi. Phytoplankton", x) == is.element("Meta. Phytoplankton", x) #YES
	return(all(c1, c2, !c3, !c4, c5))
}

Checked <- apply(AllCombinations, MARGIN=1, FUN=Checks)

FewerCombns <- AllCombinations[Checked,]
# Combinations <- FewerCombns[!duplicated(FewerCombns),]
Combinations <- unique(FewerCombns)