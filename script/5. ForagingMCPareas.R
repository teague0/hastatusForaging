#Foraging area sizes
library(adehabitatHR)

load(file="./data/4_Clustering_ForagSites.Rdata") #This only includes behaviour=="foraging" data

#Calculate 100% MCP for each cluster 
clusMCPIDs <-  mclapply(X=clusFrgFlw, mc.cores = detectCores()-2, FUN = function(x){
  tryCatch({
    z <- x %>% filter(dbClus != 0)
    ID <- unique(z$batIDday)
    clusterList <- split(z, f=z$dbClus)
    mcpList <- lapply(clusterList, function(w){
      tryCatch({
      coordinates(w) <- ~utm.x+utm.y
      proj4string(w) <- CRS("+proj=utm +zone=17 +datum=WGS84")
      mcp100 <- mcp(w, percent = 100, unin="m", unout = "m2")
      return(mcp100)
    },error=function(e) NULL)
    })
    areas <- unlist(lapply(mcpList, function(z){area <- z@data$area}))
    clusNo <- seq(1:length(areas))
    IDs <- rep(ID, length(areas))
    counts <- z %>% dplyr::select(dbClus, newState) %>% 
      dplyr::count(dbClus, newState)
    pctStatePerClus <- counts %>% group_by(dbClus) %>% 
      mutate(percent = n / sum(n))
    pctsWide <- pctStatePerClus %>% spread(newState, percent)
    pctsWide_t <- pctsWide %>% dplyr::select(pctState2 = state2, 
                                             pctState3 = state3, 
                                             pctState4 = state4)
    pctTabl <- na.omit(pctsWide_t %>% 
                         gather(state, percent,pctState2:pctState4)) %>% spread(state, percent)
    areaDat <- left_join(data.frame(IDs, clusNo, areas), pctTabl, by=c("clusNo"="dbClus"))
    return(areaDat)
  }, error=function(e) NULL)
})

foragingArea <- do.call("rbind", clusMCPIDs)
save(foragingArea, file="./data/5_ForagingMCPAreas_states.Rdata")
