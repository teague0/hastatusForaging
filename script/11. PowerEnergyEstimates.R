#Energetic costs of flights of Phyllostomus hastatus.
library(tidyverse)

load("./processed data/HastHMMdbClus_flwTrueTlagTrueGSAS.Rdata") #hast_df
source("./script/powercurve_noplot.R")
phbiodat <- read.csv("./data/Hastatus loggers Recovered.csv")

hastMorph <- left_join(hast_df, phbiodat)

#Create power curves for every bat & into a figure to show the curves for everyone
phbio_l <- split(phbiodat, f=phbiodat$batID)

allcurves <- lapply(phbio_l, function(x){
  m <- unique(x$mass.deploy)/1000
  wingspan <- unique(x$winglen.mm)*2/1000
  temp <- 25
  altitude <- 50
  minvel <- 1
  maxvel <- 30
  Vdat <- powercurve(m, wingspan, altitude, temp, minvel, maxvel)
  idVdat <- data.frame(x$batID, Vdat)
  return(idVdat)
})
curveDat <- do.call("rbind", allcurves)
#These data are used for figure S1

#Vmp & Vmr are repeated, so just summarize into a table to get those values

pwrRange <- curveDat %>% group_by(x.batID) %>% 
  summarize(Vmp = mean(Vmp),
            Vmr = mean(Vmr))
spdSums <- pwrRange %>% summarize(meanVmp = mean(MinPwrSpd),
                                  sdVmp = sd(MinPwrSpd),
                                  meanVmr = mean(MaxRgSpd), 
                                  sdVmr = sd(MaxRgSpd))

### Figure S1 ####
# powerAllInds <- ggplot(curveDat, aes(x = V, y  = Pmech, group=x.batID))+
#   geom_line(alpha = 0.3)+
#   labs(x = expression(paste("Airspeed (m ", s^-1, ")", sep="")),
#        y = "Mechanical power (W)")
#pdf("./output/Fig S1 - Powercurves.pdf", height = 6, width = 6)
#print(powerAllInds)
#dev.off()

#Get pennycuick Vmp Vmr power estimates for each bat & add to the morphometric data
phbiodat <- phbiodat %>% left_join(pwrRange, by = c("batID" = "x.batID"))


#Calculate metabolic power based on equation 1 from Ward et al 2001. 
#First, Get pennycuick power estimates
source("./script/pennypower_function.R")
m <- hastMorph$mass.deploy/1000
wingspan <- hastMorph$winglen.mm*2/1000
altitude <- hastMorph$height.above.msl
airspeed <- hastMorph$airspeedRecalc
pwr <- numeric(length = length(m))
for(i in 1:length(m)){
  pwr[i] <- pennypower(m[i], wingspan[i], altitude[i], airspeed[i])
}
hastMorph$pwrPenny <- pwr

Pmet <- numeric(length = length(m))
#Ep is the partial efficiency of the muscle. From Thomas 1975 for P. hastatus this ranges between mean value of 0.13 - 0.34. 0.23 or 0.18 are the bird values (ward 2001). 
Ep <- 0.2466667
Pmet <-  1.1*(((pwr)/Ep)+(23.8 / 360))

hastMorph$Pmet <- Pmet

#Phyllostomus metabolic power curves -- Thomas 1975 between 6-9 m/s
# 0 degrees P (W/kg) = 311.98 - 63.08V + 4.54V^2
#metPwr <- (311.98 - 63.08*airspeed + 4.54*airspeed^2)*m

## Incorporate resting metabolic rate into the state 1 data
#McNab 1969: 84.2 g; 1.19 ccO2/g/hr @ 34.7 degrees
# 1 mL O2 = 20 J
# 23.8 J/g/h - divide by 360 to get W = 0.006611111 W / g

#Any time lag over 90 s (the sleep time for missed fixes -- GPS went to sleep for 300 s, then searched for 90 s) I'm going to assume to be roosting somewhere. 

sleepLag <- which(hastMorph$tlag > 90)
hastMorph$newState[sleepLag] <- "state1"
state1s <- which(hastMorph$newState == "state1")
hastMorph$PmetTotal <- hastMorph$Pmet * hastMorph$tlag #convert W (J/s) to just joules
hastMorph$PmetTotal[state1s] <- (23.8 / 360) * hastMorph$mass.deploy[state1s]* hastMorph$tlag[state1s]

#save(hastMorph, file="./data/HastatusSegStateBiodatPwr.Rdata")

#Just for fun -- exploration of Speakman & Thomas 2003 DEE estimates via body mass scaling.
#Speakman & Thomas 2003
#Rest: log(BMR ml O2 per hour) = 1.08985 + 0.744 * logMb (g)
#Flight: Pf (W) = -1.46 + 0.744 * Mb (g)

#Speakman 2005
#ln DEE (KJ per day) = 2.05 + 0.621 * ln (Mb(g))

logSpeak <- 2.05 + 0.621 * log(phbiodat$mass.deploy)
deeSpeak <- exp(logSpeak)


bat.mass$log.Speak.KJ= 2.05 + 0.621 * log(bat.mass$Mb*1000)
bat.mass$Speak.KJ=exp(bat.mass$log.Speak.KJ)
plot(Speak.KJ~log(Mb*1000), data=bat.mass)