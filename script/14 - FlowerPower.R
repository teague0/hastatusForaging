#Flower power

#Derive the energy content per ml of balsa nectar from Kays et al 2012
# - Total nectar produced by a flower is *estimated* at 25.5 ml
# - Each flower produced 11.6 kcal per night
# - Large peak is 60 flowers per patch, but there is a normal distribution of the number of flowers open per patch over the season.
# - Flowers open with 4.9 ± 1.3 ml of nectar. There is a sharp decline in nectar production over the night. It starts at 2.6 ml / hour and then declines to 0.4 ml per hour.
# - Nectar sugar concentration decreases over the night from 13.3% at 18h to 7.9 % at 06 h. This an average of 12.4% sugar. This is 0.124 g sucrose per ml nectar
# - Roland estimates 11.75 kcal over 1 night (but we don't care about the full night)
# - 3.87 kcal/g sugar
# - 0.47988 kcal per ml (3.87*0.124)
# - 1 kcal = 4.184 kJ
# - 2.007818 kJ / ml nectar
# - peak flower availability is about 60 per patch or 5 per m2. Mean is more like 20 & 2.5 flowers per m2
# - 602.3454 kJ per crown (2.007818 kJ   / ml * 5 ml / flower * 60 flowers per crown) 

#If we assume that bats drain 5 ml of the most profitable nectar

load("./data/HastatusSegStateBiodatPwr.Rdata") #hastMorph
load("./data/NightSummariesNets.Rdata") #nighVals
load("./data/NightSumValues.Rdata") #nightNets


#How many ml of nectar & # of flowers would a bat need to power RMR + flight power for the day
#This can be added to the nightNets summary table

nightNets$mlRequired <- nightNets$dee.kJ / 2.007818
nightNets$NFlowers <-  nightNets$mlRequired / 5


#Cumulative energy expenditure per bat day.
nas <- which(is.na(hastMorph$PmetTotal))
hastMorph$PmetTotal[nas] <- 0
tlagNas <- which(is.na(hastMorph$tlag))
hastMorph$tlag[tlagNas] <- 0

batDay <- split(hastMorph, f=hastMorph$batIDday)

cumEE <- lapply(batDay, function(x){
  x$cumEE.kJ <- cumsum(x$PmetTotal)/1000
  x$timeElapsed.s <-  cumsum(x$tlag)
  return(x)
})

hastMorph <- do.call("rbind", cumEE)
mlNeeded <- hastMorph$cumEE.kJ / 2.007818
flwrNeeded <- mlNeeded / 5

#Data used for Figure 4, Figure S2
