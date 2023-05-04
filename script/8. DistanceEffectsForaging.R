#What do patch use & residency times look like with and without neighbors. 
library(tidyverse)
library(lubridate)
library(parallel)
library(cowplot)
theme_set(theme_cowplot())

## ** NOTE ** this data frame has duplicated locations to accommodate the multiple pairwise distances that were added when there were more than 2 bats tracked together on a night. 

load("./data/7_HastatusSegClusRelevel.Rdata")

forage <- hastP %>% filter(behaviour=="forage")
otherBatPatch <- forage %>% dplyr::select(timestamps, batID, otherBatPatch = patchID)
forage <- forage %>% left_join(otherBatPatch, by=c("timestamps" = "timestamps", "otherBat" = "batID")) 

#What are the distances between bats when they are & aren't in the same patch?
partners <- forage
partners$samePatch <- NA
for(i in 1:length(partners$patchID)){
  partners$samePatch[i] <- ifelse(partners$patchID[i] == partners$otherBatPatch[i], 1, 0)
} 



#Who are relatively close to each other? Calls with peak frequency of 6725 ± 36.3 Hz travel 289.8 m using the echolocation formula, so about 2x that for another receiver. Assume a receiver communication distance of about 500 m to be very conservative (social effects at longer distances)

#close <- which(partners$otherDist < 500)
close <- which(partners$otherDist < 290)
partnersClose <- partners[close,]


#Are multiple groups tracked at the same time? No. There is no overlap among groups of tracking nights. Can only look at within group effects and then the lasting effect of the location. This could let us test if patch use is group specific or if everyone finds the same places. 


save(partnersClose, file="./data/8_Hastatus_PatchPartnersIDd.Rdata")

## These results are used in part for Figure 2.

