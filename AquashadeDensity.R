# Aquashade has 1.11 lbs of Acid Blue 9 (Erioglaucine) and 0.09 lbs of Acid Yellow 23 (Tartrazine) per gallon of product.

#there is (1/1000)*1 mL of Aquashade in a 1 L volume of water+Aquashade at 1ppm (the 1ppm mark for the auqashade is a volume fraction)

#therefore, at 1ppm, the mass of aquashade is .......

a0Gal <- 1.2 #pounds of pure aquashade in a gallon container that presumably contains water and other stuff
wGal <- 8.33 # pounds of water in a gall on pure water

#those 1.2 pounds were apparently 13.64% of the volume of the container
# *if* the container was pure water, it would weight 8.33 pounds

#so 86.36% water is 
w0Gal <- 0.8636 * wGal #pounds of water in a container of Aquashade in the 86.36% that isn't dye is pure water

#another way, if that volume fo aquashade was just water, it would weight 1.136212 pounds, which is (1.2 - 1.136216)/1.136216 = 0.05613721 pounds off from the density of water per container of aquashade

cMass <- w0Gal + a0Gal #mass of contents of Aquashade container, given the above assumptions

(cMass-wGal)/wGal #so the density **averaged over all the contents in the container** is 0.7% higher than the density of water ... so instead of 1.5 mg/L, it would 


a0Gal/1E6*1.5 #number of pounds of pure aquashade in a 1.5ppm solution ....

a0Gal/1E6*1.5*544.3104*1000 #number of mg of pure aquashade in a 1.5ppm (volume fraction) solution that has a volume of 1 gallon

a0Gal/1E6*1.5*544.3104*1000/3.78541 #number of mg of pure aquashade in a 1.5ppm (volume fraction) solution that has a volume of 1 gallon
