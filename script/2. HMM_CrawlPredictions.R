#Using crawl to regularize tracks based on a correlated random walk with a Kalman filter. This will help with the HMM segmentation to assign behavioral modes to the GPS locations.
library(momentuHMM)
library(tidyverse)
library(move)
library(lubridate)
library(parallel)

#use object hast_df from script 1 or load in written csv file.
hast_df <- read.csv("./data/Phyllostomus Hastatus Processed GPS Data S1.csv")
hast_df$timestamp <- as.POSIXct(hast_df$timestamp, tz = "UTC")
hast_df$timestamps <- as.POSIXct(hast_df$timestamps, tz = "UTC")

#Re-classify some points for "X74F9F83.4" that are on Colón at the begin of commute
forageDat <- hast_df %>% filter(behaviour=="forage") #use foraging data only
hast_L <- split(hast_df, f=hast_df$batIDday)
lapply(hast_L, function(x){
  time_length(interval(min(x$timestamps), max(x$timestamps)), unit="mins")}) #check the interval of how long each bat was tracked per night


#### WARNING: Go get coffee or a snack because this will take a few minutes. If running on a slower machine, then a leisurely meal is recommended ####

hastLcrw <- mclapply(X=hast_L, mc.cores=detectCores()-2, FUN=function(batHMM) { 
  tryCatch({
    tmp <- batHMM %>% dplyr::select(ID=batIDday, time=timestamps, x=utm.x, y=utm.y)
    id <- unique(tmp$ID)
    crwOut <- crawlWrap(obsData=tmp, 
                        Time.name="time", 
                        timeStep="1 secs",
                        theta=c(7, 0.1), #velocity variation, velocity autocorrelation
                        fixPar=c(NA,NA),
                        retryFits = 5) 
    return(crwOut)
  },error=function(e) finally = print(paste(id, "f'd up")))
})

#save(hastLcrw, file = "./data/Hastatus Crawl Fits.Rdata")

#Now move on to the HMM fits.
