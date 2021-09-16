#Spatial social network.

library(spatsoc)
library(data.table)
library(tidyverse)
load("./data/9_HastHMMdbClus_flwTrueTlagTrueGSAS.Rdata")

groupIDS <- hast_df %>% group_by(groupID, batID) %>% 
  dplyr::summarize(nlocs = length(location.lat)) %>% 
  arrange(batID) %>% 
  as.data.frame()
groupIDS$groupNo <- recode(groupIDS$groupID, blue = 1, brown = 2, yellow = 3)

batNos <- groupIDS %>% arrange(groupID, batID) 
batNos$batNo = seq(1,18, by=1)
batNos <- batNos %>% arrange(batID) 

#Generate a single group spatio-temporal network ####
g3 <- hast_df %>% filter(groupID == "yellow")
g2 <- hast_df %>% filter(groupID == "brown")
g1 <- hast_df %>% filter(groupID == "blue")

hast_dt <- as.data.table(hast_df)
group_times(hast_dt, datetime = "timestamps", threshold = "1 minute")
group_pts(hast_dt, threshold = 30, id = "batID", coords = c("utm.x", "utm.y"), timegroup = "timegroup")
gbi <- get_gbi(hast_dt, group = "group", id = "batID")

library(asnipe)
library(igraph)
net <- get_network(gbi, data_format = "GBI", association_index = "SRI")

## Plot the spatial graph
g <- graph.adjacency(net, 'undirected',
                     diag = FALSE, weighted = TRUE)


#R is dumb and I have to create a triangle shape ####
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/110 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

#Assign shape and color to each bat ID####
coul <- viridisLite::viridis(6)[c(3, 4, 6)]
shape <- c("square", "circle", "triangle")
node.colors <-  coul[groupIDS$groupNo]
node.shape <- shape[groupIDS$groupNo]
set.seed(2)

#Get the layout. Manipulate in tkplot & then save the coordinates. 
#Guide: https://kateto.net/netscix2016.html
g_tk <- tkplot(g, vertex.color=node.colors, vertex.label = batNos$batNo)
lay <- tk_coords(g_tk)

## Metrics for all individuals in all group spatial network ####
observed <- data.table(
  centrality = evcent(g)$vector,
  strength = graph.strength(g),
  degree = degree(g),
  ID = names(degree(g))
)
#save(observed, file="./data/10_ObservedSocNetworkMetrics.Rdata")

#Network of patch use similarities ####

#Generate a 0/1 matrix for the patch occupany
patchSums <- hast_df %>% group_by(batID, patchID) %>% 
  dplyr::summarise(sumTime = ifelse(sum(tlag, na.rm=T)>0, 1, 0)) %>%
  #summarise(sumTime = sum(tlag, na.rm=T)) %>%
  spread(patchID, sumTime) %>%
  mutate_if(is.numeric, funs(replace_na(., 0))) %>%
  as.data.frame()

#asnipe wants a matrix of group by individual instead of individual by group (which is what we have). This needs to be transformed. First the bat IDs are going to be set as rownames, the NA column is taken out, and then transformation happens.

rownames(patchSums) <- patchSums$batID  
patchSums <- patchSums[,2:64]
patchSums.t <- t(patchSums)
adj.m <- get_network(patchSums.t, association_index = "SRI")  
gP <- graph_from_adjacency_matrix(adj.m, "undirected", weighted=T) #the igraph object
coms <- cluster_fast_greedy(gP)

gr1 <- which(coms$membership==1)
gr2 <- which(coms$membership==2)
gr3 <- which(coms$membership==3)

set.seed(2)
g1 <- data.frame("batID" = batNos$batID[gr1], "patchCom" = 1)
g2 <- data.frame("batID" = batNos$batID[gr2], "patchCom" = 2)
g3 <- data.frame("batID" = batNos$batID[gr3], "patchCom" = 3)
patchCommun <- rbind(g1, g2, g3)

#Add in the social network variables to the NightSums from 9. NightSummaries
load("./data/9_NightSummaries.Rdata")
socNet <- as.data.frame(observed)
nightVals <- nightVals %>% left_join(socNet, by = c("batID" = "ID"))
nightVals <- nightVals %>% left_join(patchCommun)
#save(nightVals, file = "./data/10_NightSummariesNets.Rdata")

V(gP)$color <- node.colors
gP_tk <- tkplot(gP, vertex.label = batNos$batNo)
layP <- tk_coords(gP_tk)

#save(lay, layP, file = "./data/NetworkLayouts.Rdata")

plot(gP, layout=layP, 
     vertex.label = batNos$batNo, 
     vertex.shape = node.shape,
     edge.width=E(g)$weight*100, 
     edge.color="grey50", 
     edge.curved=0.3, 
     mark.groups = list(gr1, gr2, gr3)) #plot it
observedPatch <- data.table(
  centrality = evcent(gP)$vector,
  strength = graph.strength(gP),
  degree = degree(gP),
  ID = names(degree(gP))
)
#save(observedPatch, file="./data/10_ObservedPatchSocNetworkMetrics.Rdata")


pdf("./output/SpatialPatchNetworks.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
plot(g, vertex.color=node.colors, vertex.shape = node.shape,
     layout=lay, vertex.label = batNos$batNo,
     edge.width=E(g)$weight*50, 
     edge.color="grey50", 
     edge.curved=0.2)
par(adj=0)
title("A")
plot(gP, layout=layP, 
     vertex.label = batNos$batNo, 
     vertex.shape = node.shape,
     edge.width=E(g)$weight*10, 
     edge.color="grey50", 
     edge.curved=0.1, 
     mark.groups = list(gr1, gr2, gr3))
par(adj=0)
title("B")
dev.off()



