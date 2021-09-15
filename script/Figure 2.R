#Figure 2
library(tidyverse)
library(lme4)
library(cowplot)
theme_set(theme_cowplot())

load("./data/Hastatus_PatchPartnersIDd.Rdata") #partners from 8. DistancEffectsForaging
load("./data/HastatusSegStateBiodatPwr.Rdata") #hast_df from 11. PowerEnergyEstimates
groups <- hast_df %>% dplyr::group_by(batID) %>% 
  dplyr::summarize(groupID = unique(groupID))
partners <- partners %>% left_join(groups)

mycols <- viridisLite::viridis(6)[c(3, 4, 6)]

actBudgIn <- partners %>% filter(!is.na(samePatch)) %>% 
  group_by(groupID, batID, date, samePatch, newState) %>% 
  dplyr::summarize(nlocs = n()) %>% 
  mutate(stateFreq = nlocs / sum(nlocs)) %>% 
  filter(samePatch == 1) %>% 
  ggplot( )+
  geom_half_boxplot(aes(x = newState, y = stateFreq, fill = groupID), nudge = 0.02, outlier.color = NA)+
  geom_half_point(aes(x = newState, y = stateFreq, color = groupID), size = 2, shape = 16)+ 
  scale_x_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("Rest", "Slow\nflight" ,"Move", "Commute"))+
  ylim(0, 0.075)+
  scale_fill_manual(values = mycols[2:3], 
                    breaks = c("brown", "yellow"),
                    labels = c("Group 2", "Group 3"))+
  scale_color_manual(values = mycols[2:3], 
                     breaks = c("brown", "yellow"),
                     labels = c("Group 2", "Group 3"))+
  labs(x = "", 
       y = "frequency", 
       subtitle = "Activity budgets when animals\n occupy the same patches")+
  theme(legend.position = "none")


actBudgOut <- partners %>% filter(!is.na(samePatch)) %>% 
  group_by(groupID, batID, date, samePatch, newState) %>% 
  dplyr::summarize(nlocs = n()) %>% 
  mutate(stateFreq = nlocs / sum(nlocs)) %>% 
  filter(samePatch == 0) %>%
  ggplot( )+
  geom_half_boxplot(aes(x = newState, y = stateFreq, fill = groupID), nudge = 0.02, outlier.color = NA)+
  geom_half_point(aes(x = newState, y = stateFreq, color = groupID), size = 2, shape = 16)+ 
  scale_x_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("Rest", "Slow\nFlight" ,"Move", "Commute"))+
  ylim(0, 0.075)+
  scale_fill_manual(values = rep(mycols, 3), 
                    name = "", 
                    breaks = c("blue", "brown", "yellow"),
                    labels = c("Group 1", "Group 2", "Group 3"))+
  scale_color_manual(values = rep(mycols, 3), 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  labs(x = "", 
       y = "frequency", 
       subtitle = "Activity budgets when animals occupy\n diffferent patches")+
  theme(legend.position = c(0.05, 0.86), 
        legend.background = element_blank(),
        legend.title = element_blank())

# What about synchronization of behavior in a patch? ####
load("./processed data/AllBatsSamePatch 2020-07-06.Rdata")
inPatch <- forage %>% filter(samePatch == 1)
inPatch$batStateRecode <- NA
inPatch$otheBatStateRecode <- NA
inPatch$sameBehav <- NA
for(i in 1:length(inPatch$newState)){
  inPatch$batStateRecode[i] <- ifelse(inPatch$newState[i] == "state2", "forage", 
                                      ifelse(inPatch$newState[i] == "state3", "forage", inPatch$newState[i]))
  inPatch$otheBatStateRecode[i] <- ifelse(inPatch$otherBatState[i] == "state2", "forage", 
                                          ifelse(inPatch$otherBatState[i] == "state3", "forage", inPatch$otherBatState[i]))
}

inPatch$sameBehav <- ifelse(inPatch$batStateRecode == inPatch$otheBatStateRecode,"sameBehav", "diffBehav")
inPatch$sameBehav_orig <- ifelse(inPatch$newState == inPatch$otherBatState,"sameBehav", "diffBehav")


sameSums <- inPatch %>% 
  dplyr::group_by(groupID, batID, date, newState, sameBehav_orig) %>% 
  dplyr::summarize(nObs = n()) %>% 
  ungroup() %>% 
  dplyr::group_by(groupID, batID, date, newState) %>% 
  dplyr::mutate(totalN = sum(nObs), 
                freq = nObs/totalN) %>% 
  dplyr::filter(sameBehav_orig == "sameBehav")

#Lump the synchrony into rest vs other
sameSums <- sameSums %>% mutate(restOther = ifelse(newState == "state1", "rest", "other" ))

mS <- glmer(freq~restOther+(1|batID), family = binomial, data = sameSums)
summary(mS)
Anova(mS)


library(gghalves)
behavSameFreq <- sameSums %>% 
  ggplot()+
  geom_half_boxplot(aes(x = newState, y = freq, fill = groupID), nudge = 0.02, outlier.color = NA) +
  geom_half_point(aes(x = newState, y = freq, color = groupID), size = 2, shape = 16)+ 
  scale_x_discrete(breaks = c("state1", "state2", "state3", "state4"),
                   labels = c("Rest", "Slow\nflight" ,"Move", "Commute"))+
  scale_fill_manual(values = mycols[2:3], 
                    breaks = c("brown", "yellow"),
                    labels = c("Group 2", "Group 3"))+
  scale_color_manual(values = mycols[2:3], 
                     breaks = c("brown", "yellow"),
                     labels = c("Group 2", "Group 3"))+
  labs(x = "", 
       y = "frequency", 
       subtitle = "Behavioral sychrony")+
  theme(legend.position = c(0.8, 0.8), 
        legend.title = element_blank())+
  ylim(-0.1, 1.1)

# How does proximity affect rest? ####
#This will start with a behavioral change point analysis to count up behavioral sequences.
# I think I can use the patchStay df from 9. Patch Residency Influence since it split up both the in/out of a patch and the behavior in the same loop by bat day. This only was done on the 'forage' behavior class (so no commuting included). This df has duplicates b/c all bat pairs are included for distance, patch effects, etc. I'll keep this analysis under 1000 m. 

load("./processed data/AllBatsSamePatchStateSeq 2020-07-06.Rdata") #patchStay

behavBoutTimeDist <- patchStay %>% filter(otherDist < 100) %>% 
  group_by(batID_day, stateBout, newState) %>% 
  dplyr::summarize(timeElapsed = difftime(max(timestamp), min(timestamp), units = "secs"),
                   meanDistNN = mean(otherDist, na.rm = TRUE),
                   minDistNN  = min(otherDist, na.rm = TRUE),
                   sdDistNN = min(otherDist))
stateNo <- c("state1", "state2", "state3", "state4")
behavName <- c("rest", "slow flight" ,"move", "commute")
names(behavName) <- stateNo
#At least with this figure, being in close proximity increases the duration of rest. Try to find an non-linear fit for rest
restProx <- behavBoutTimeDist %>% filter(newState == "state1")

library(aomisc)
model <- drm((as.numeric(timeElapsed)/60) ~ minDistNN, fct = DRC.powerCurve(),
             data = restProx)

# predictions and confidence intervals.
demo.fits <- expand.grid(minDistNN=seq(0.1, 100, length=100))
# new data with predictions
pm <- predict(model, newdata=demo.fits, interval="confidence") 
demo.fits$time <- pm
demo.fits$pmin <- pm[,2]
demo.fits$pmax <- pm[,3]

#Plot
distRestDur <- ggplot(restProx)+
  geom_point(aes(x = minDistNN, y = as.numeric(timeElapsed)/60), alpha = 0.6)+
  geom_line(aes(x = minDistNN, y = time), data = demo.fits, color = "red")+
  labs(x = "Distance (m)", 
       y = "Duration (min)", 
       subtitle = "Proximity influences on resting bout duration")

top_row <- plot_grid(actBudgOut, actBudgIn, labels = c('A', 'B'), label_size = 12, ncol = 2)
bottom_row <- plot_grid(behavSameFreq, distRestDur, labels = c('C', 'D'), label_size = 12, ncol = 2)

pdf("./output/Fig 2 - Distance Effects on Behavior.pdf", width = 8, height = 8)
plot_grid(top_row, bottom_row, label_size = 12, nrow = 2)
dev.off()



