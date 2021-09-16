#Find clusters of GPS locations to then ID behavioral states.

library(tidyverse)
library(adehabitatLT)
library(dbscan)
library(move)
library(fpc)
library(factoextra)
library(parallel)
library(cowplot)
theme_set(theme_cowplot())

#https://www.datanovia.com/en/lessons/dbscan-density-based-clustering-essentials/#why-dbscan

load("./data/3_HastatusGPS_HMMRelevelstates.Rdata")

#One set of values for eps & minpts doesn't really work for all of the tracks. Looks like the break into eps = 5, 10, or 20. Plot diagnostic plots to help sort this out. This is commented out, but run if needed
# mclapply(X=newStates, mc.cores = detectCores()-2, FUN = function(x){
#   tryCatch({
#     tmp <- x %>% dplyr::filter(behaviour=="forage")
#     xy <-  tmp %>% dplyr::select(utm.x, utm.y)
#     ID <- unique(x$batIDday)
#     pdf(file = paste0("./output/",ID, "_dbscanEPSck.pdf"))
#     dbscan::kNNdistplot(xy, k = 10)
#     abline(h=5)
#     abline(h=10, lty=3)
#     abline(h=20, lty=2)
#     dev.off()
#   },error=function(e) finally = print("No Forage"))
# })

#Assign custom EPS based on the diagnostics from above
ids <- names(hast_L)
eps <- rep(7, length(ids)) #turns out it mostly should have been 10
eps[c(5,6,7,8,9,10,11,14,15,16,17,19,21,22,23,26,27,29)] <- 10
eps[25] <- 8
eps[c(1,10)] <- 12
eps[c(1,2,13,24)] <- 20
eps[3] <- 4
eps[5] <- 22
eps[8] <- 14
eps_df <- data.frame(ids, eps)
eps_df$pts <- 15
eps_df$eps[8] <- 14 
eps_df$eps[28] <- 4 
eps_df$pts[8] <- 10

clus <- mclapply(X=newStates, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
  tmp <- x %>% dplyr::filter(behaviour=="forage")
  xy <-  tmp %>% dplyr::select(utm.x, utm.y)
  ID <- unique(x$batIDday)
  epsVal <- eps_df$eps[which(eps_df$ids == ID)]
  ptsVal <- eps_df$pts[which(eps_df$ids == ID)]
  set.seed(123)
  db <- fpc::dbscan(xy, eps = epsVal, MinPts=ptsVal)
  tmp$dbClus <- db$cluster
  return(tmp)
  },error=function(e) finally = print("No Forage"))
})

#To further cluster the data to see if there are any identifiable feeding sites (or just all noise), need to pull out the points that are in the lower speed HMM categories (1 & 2). I don't want to re-run HMMs here yet since there might be some crazy long time lags for animals that leave the patch for awhile and return. 


#Now recluster each larger cluster to ID potential sites of long residency or flowers.
clusFrgFlw <- mclapply(X=clus, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
    w <- x %>% filter(dbClus != 0)
    wlist <- split(w, f=w$dbClus)
    tmp <- lapply(wlist, function(newClus){
      tryCatch({
      set.seed(1234)
      xy2 <- newClus %>% filter(newState %in% c("state2")) %>% dplyr::select(utm.x, utm.y)
      xy3 <- newClus %>% filter(newState %in% c("state3")) %>% dplyr::select(utm.x, utm.y)
      xy4 <- newClus %>% filter(newState %in% c("state4")) %>% dplyr::select(utm.x, utm.y)
      tryCatch({db2 <- fpc::dbscan(xy2, eps = 0.8, MinPts=6)},error=function(e) NULL)
      tryCatch({db3 <- fpc::dbscan(xy3, eps = 0.8, MinPts=6)},error=function(e) NULL)
      tryCatch({db4 <- fpc::dbscan(xy4, eps = 0.8, MinPts=6)},error=function(e) NULL)
      xy2event.id <- newClus %>% filter(newState %in% c("state2")) %>% dplyr::select(event.id)
      xy2$event.id <- xy2event.id$event.id
      xy3event.id <- newClus %>% filter(newState %in% c("state3")) %>% dplyr::select(event.id)
      xy3$event.id <- xy3event.id$event.id
      xy4event.id <- newClus %>% filter(newState %in% c("state4")) %>% dplyr::select(event.id)
      xy4$event.id <- xy4event.id$event.id
      xy2$dbClus_feed2 <- NA
      xy3$dbClus_feed3 <- NA
      xy4$dbClus_feed4 <- NA
      tryCatch({xy2$dbClus_feed2 <- db2$cluster},error=function(e) NULL)
      tryCatch({xy3$dbClus_feed3 <- db3$cluster},error=function(e) NULL)
      tryCatch({xy4$dbClus_feed4 <- db4$cluster},error=function(e) NULL)
      dat2 <- left_join(newClus, xy2)
      dat3 <- left_join(dat2, xy3)
      dat4 <- left_join(dat3, xy4)
      return(dat4)
      },error=function(e) NULL)
    })
    tmpdat <- do.call("rbind", tmp)
    clusterVals <- tmpdat %>% dplyr::select(event.id, dbClus_feed2, dbClus_feed3, dbClus_feed4)
    dat <- left_join(x, clusterVals, by="event.id")
    return(dat)
  },error=function(e) NULL)
  }) 

hastAll_df <- do.call("rbind", newStates)
forage_df <- do.call("rbind", clusFrgFlw)
clusters <- forage_df %>% dplyr::select(event.id, dbClus:dbClus_feed4)
hastClass_df <- left_join(hastAll_df, clusters)
save(hastClass_df, file="./data/4_HastatusHMMdbClus_flw.Rdata")
save(clusFrgFlw, file="./data/4_Clustering_ForagSites.Rdata")
