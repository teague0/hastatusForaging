#Calculate some basic behavioral descriptions for each bat night.
library(sf)
library(tidyverse)
library(lubridate)
library(geosphere)
library(raster)

load("./data/7_HastatusSegClusRelevel.Rdata")
#The above data have duplicated records that were used to calculate distances among all group members. Remove the duplicates
dups <- which(duplicated(hastP$event.id) ==FALSE) 
hastStart <- hastP[dups,]


#Recalculate timelag by batDay. This eliminates very long timelags that were calculated per animal ID ####
#Timelag goes to the second point. 
hastL <- split(hastStart, f = hastStart$batIDday)

hastL_tl <- lapply(hastL, function(x){
  timelags <- c(NA, difftime(x$timestamps[-1], x$timestamps))
  last <- length(timelags)
  timelags <- timelags[1:(last-1)]
  x$tlag <- as.numeric(timelags)
  coords <- SpatialPoints(cbind(x$utm.x, x$utm.y), proj4string = CRS("+proj=utm +zone=17 +datum=WGS84"))
  x$stepLength <- NA
  for(i in 1:(length(x$utm.x)-1)){
    j = i+1
    x$stepLength[j]  = pointDistance(coords[i], coords[j])
  } 
  return(x)
 })
hast_df <- do.call("rbind", hastL_tl)

#Calculate all GPS location distances to La Gruta. Load coordinates of cave to calculate distance
coords <- cbind(-82.271541, 9.396448)
lagruta <- SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=WGS84"))
lagrutaUTM <- spTransform(lagruta, CRS("+proj=utm +zone=17 +datum=WGS84"))
hast_df$ptDistToCave <- spDists(SpatialPoints(cbind(hast_df$utm.x, hast_df$utm.y), proj4string = CRS("+proj=utm +zone=17 +datum=WGS84")), lagrutaUTM)

hast_df$groundSpeedRecalc <- hast_df$stepLength / hast_df$tlag
hast_df$airspeedRecalc <- sqrt((hast_df$groundSpeedRecalc - hast_df$ws)^2 + (hast_df$cw)^2)
save(hast_df, file="./data/9_HastHMMdbClus_flwTrueTlagTrueGSAS.Rdata")

#Calculate time budgets ####
timeRanges <- hast_df %>% group_by(groupID, batID, batIDday) %>% 
  dplyr::summarise(minTimeUTC = min(timestamps),
            maxTimeUTC = max(timestamps),
            minTimeLocal = min(with_tz(timestamps, tz="America/Panama")),
            maxTimeLocal = max(with_tz(timestamps, tz="America/Panama")),
            timeTrack.min = round(as.numeric(as.duration(maxTimeLocal-minTimeLocal))/60, 2),
            nlocs = n(),
            meanTimeLag.s = mean(tlag, na.rm = TRUE),
            maxTimeLag.s = max(tlag, na.rm = TRUE),
            totalDistance = sum(stepLength, na.rm = TRUE))

#Get sunrise/sunset times
library(maptools)
lagrutaLL <- spTransform(lagrutaUTM, CRS("+proj=longlat +datum=WGS84"))
sunset <- sunriset(coordinates(lagrutaLL), timeRanges$minTimeLocal, direction="sunset", POSIXct.out=TRUE)[2]
sunrise <- sunriset(coordinates(lagrutaLL), timeRanges$maxTimeUTC, direction="sunrise", POSIXct.out=TRUE)[2]
timeRanges$timeSinceSunset <- timeRanges$minTimeLocal - sunset$time
timeRanges$timeToSunrise <- timeRanges$maxTimeUTC - sunrise$time

#Distance of first and last point of the night to LG
startDist <- hast_df %>% group_by(batIDday) %>%
  filter(timestamps == min(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batIDday, firstFixDist = ptDistToCave)
lastDist <- hast_df %>% group_by(batIDday) %>%
  filter(timestamps == max(with_tz(timestamps, tz="America/Panama"))) %>%
  dplyr::select(batIDday, lastFixDist = ptDistToCave)
firstLast <- startDist %>% left_join(lastDist)
timeRanges <- timeRanges %>% left_join(firstLast)


overallSums <- timeRanges %>% ungroup() %>% 
  summarize(meanFirstTime = mean(timeSinceSunset, na.rm=T),
            sdFirstTime = sd(timeSinceSunset, na.rm=T), 
            meanFirsLocDist = mean(firstFixDist[,1], na.rm=T),
            sdFirsLocDist = sd(firstFixDist[,1], na.rm=T),
            meanTimeTrack = mean(timeTrack.min, na.rm=T), 
            sdTimeTrack = sd(timeTrack.min, na.rm=T), 
            meanTotalDist = mean(totalDistance, na.rm=T),
            sdTotalDist =sd(totalDistance, na.rm=T))

#Speeds. Ground speeds over 29 seem to be errors, so let's get rid of those.
outliers <- which(hast_df$groundSpeedRecalc > 29)
hast_df$groundSpeedRecalc[outliers] <- NA


speedSum <- hast_df %>% group_by(behaviour) %>% 
  filter(newState != "state1") %>% 
  summarize(meanGrSpeed = mean(groundSpeedRecalc, na.rm=T),
            sdGrSpeed = sd(groundSpeedRecalc, na.rm=T),
            meanAirspeed = mean(airspeedRecalc, na.rm=T),
            sdAirpeed = sd(airspeedRecalc, na.rm=T),
            ws = mean(ws, na.rm=T))

#Number patches used per bat per night? ####

patchUsed <- hast_df %>% group_by(batID, batIDday, patchID) %>% 
  mutate(Npatch = length(patchID), 
         meanPatchDistCave = mean(patchDistToCave, na.rm = TRUE))

nightVals <- timeRanges %>% left_join(firstLast)
nightVals <- nightVals %>% left_join(patchUsed)
save(nightVals, file="./data/9_NightSummaries.Rdata")
