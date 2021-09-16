#Social influences on foraging & resting?
#number of foraging / feeding positions when group members are nearby? How much resting? 

#mean foraging area size is 3522.231 ± 6226.811 m^2 (circle diameter of 66 m)

library(tidyverse)
library(lubridate)
library(raster)
library(plyr)

load("./processed data/4_HastatusHMMdbClus_flw.Rdata")

#I want to calculate the pairwise distance between all bats on a given day. Seems like the data should be split by date and then distances calculated.

keepTS <- hastClass_df %>%
  dplyr::group_by(timestamps) %>%
  dplyr::summarize(count = length(utm.x)) %>% 
  dplyr::filter(count > 1)
batsKeep<- hastClass_df %>% dplyr::filter(timestamps %in% keepTS$timestamps)

#Calculate pairwise distances in meters for original data
pair.dist <- ddply(batsKeep, 'timestamps', function(t){
  combi <- t(combn(t$batID, 2)) # all unique pairwise combinations of individuals
  combi <- data.frame(bat1=combi[,1], bat2=combi[,2], stringsAsFactors=F) # set up a new data frame
  combi <- combi[combi$bat1!=combi$bat2,]# just in case that a single bat has two readings for the same timestamp
  
  combi$geoDist <- unlist(lapply(1:nrow(combi), function(i){ # calculate geographic distance between all pairs (in meters!!)
    pos.bat1 <- t[t$batID==combi$bat1[i], c('utm.x', 'utm.y')] # location for first bat
    pos.bat2 <- t[t$batID==combi$bat2[i], c('utm.x', 'utm.y')] # location for second bat
    return(pointDistance(pos.bat1, pos.bat2, lonlat = FALSE))})) # geographic distance between bats (in meters)
  combi$n.bats <- length(unique(t$batID)) # number of bats observed for this timestamp
  return(combi)})

#The distances are given as unique combinations on each row (so the bat1-bat2 combo won't be repeated as bat2-bat1). To add this info to the main data, need to combine by timestamp and bat 1 and then bat2.

dist2 <- pair.dist %>% dplyr::rename(c("bat3"="bat2", "geoDist2"="geoDist")) %>% dplyr::select(-n.bats)
dat <- left_join(hastClass_df, pair.dist, by=c("timestamps" = "timestamps", "batID"="bat1"))
twos <- which(dat$n.bats ==2)

dat2 <- dplyr::left_join(dat, dist2, by=c("timestamps" = "timestamps", "batID"="bat3"))
#pull the pairs together into a single column.
dat2$otherBat <- dplyr::coalesce(dat2$bat2, dat2$bat1)
dat2$otherDist <- dplyr::coalesce(dat2$geoDist, dat2$geoDist2)
#clean out those other columns and some other movebank stuff we really don't need.

hastClass_df <- dat2
#save(hastClass_df, file="./data/6_HastatusClusDist.Rdata")





