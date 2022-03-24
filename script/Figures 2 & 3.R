#Figure 2. Effects of proximity on bat activity budgets. A) Bat activity budgets when animals occupied different patches, and B) when they occupied the same patch. C) the frequency of observations of synchronized behavior when bats occupied the same patch. Note that Group 1 members never occupied the same patch together. D) Duration of continuous resting bouts relative to the distance to a bat's nearest neighbor

library(tidyverse)
library(lme4)
library(cowplot)
theme_set(theme_cowplot())
library(gghalves)

load("./data/8_Hastatus_PatchPartnersIDd.Rdata") #partnersClose
load("./data/11_HastatusSegStateBiodatPwr.Rdata") #hastMorph

mycols <- viridisLite::plasma(9)[c(1, 5, 8)]

#What is the activity budget when individuals are in the same patch. SamePatch was defined in script 8, with 1 = same patch, 0 = other patch. Activity budgets should be calculated as _in the patch_ instead of total activity budget

inPatchActBudg <- partnersClose %>% filter(!is.na(samePatch)) %>% 
  filter(samePatch == 1) %>% 
  group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(stateLocs = n()) %>% 
  ungroup() %>% 
  dplyr::group_by(groupID, batID, batIDday) %>% 
  dplyr::mutate(Nobs = sum(stateLocs), 
                freq = stateLocs/Nobs)
inPatchActBudg$inOut <- "in"

outPatchActBudg <- partnersClose %>% filter(!is.na(samePatch)) %>% 
  filter(samePatch == 0) %>% 
  group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(stateLocs = n()) %>% 
  ungroup() %>% 
  dplyr::group_by(groupID, batID, batIDday) %>% 
  dplyr::mutate(Nobs = sum(stateLocs), 
                freq = stateLocs/Nobs)
outPatchActBudg$inOut <- "out"

comboActBudg <- inPatchActBudg %>% bind_rows(outPatchActBudg)
comboActBudg$stateInOut <- paste0(comboActBudg$newState, comboActBudg$inOut)
mAb <- glmer(freq~newState*inOut+(1|batID:groupID), family = "binomial", data = comboActBudg)
summary(mAb)
Anova(mAb, type = 2)
anova(mAb, type = 2, ddf = "Satterthwaite")
r.squaredGLMM(mAb)

greycols <- c("#808080", "#DCDCDC")

actBudgInOut <- comboActBudg %>% 
  ggplot()+
  geom_half_boxplot(aes(x = stateInOut, y = freq, fill = inOut), alpha = 0.6)+
  geom_half_point(aes(x = stateInOut, y = freq, fill = inOut, color = groupID), position=position_jitter(width = 0.1))+
  scale_x_discrete(breaks = c("state1in", "state1out", "state2in","state2out", "state3in","state3out","state4in","state4out"),
                   labels = c("         rest", "", "         slow\n          flight","","          move","","             commute",""))+
  scale_fill_manual(values = greycols,
                     name = "", 
                     breaks = c("in", "out"),
                     labels = c("Same\nPatch", "Different\nPatch"))+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  labs(x = "", 
       y = "proportion of time",
       subtitle = "Activity budgets of when nearest neighbors\nare in the same or a different patch")+
  theme(legend.position = c(0.75, 0.7), 
        legend.background = element_blank(),
        legend.title = element_blank(),
        axis.title.x = element_text(hjust=1),
        axis.ticks.x = element_blank())
actBudgInOut

m.restAb <- glmer(freq~inOut+(1|batID:groupID), family = "binomial", data = comboActBudg[comboActBudg$newState == "state1",])
summary(m.restAb)
Anova(m.restAb)


comboActBudg %>% dplyr::filter(newState == "state1") %>%
  dplyr::group_by(inOut) %>% 
  dplyr::summarize(median = median(freq),
                   mad = mad(freq))

m.slowAb <- lm(freq~samePatch, family = "binomial", data = comboActBudg[comboActBudg$newState == "state2",])
summary(m.slowAb)
anova(m.slowAb)

m.moveAb <- lm(freq~samePatch, family = "binomial", data = comboActBudg[comboActBudg$newState == "state3",])
summary(m.moveAb)
anova(m.moveAb)

m.commAb <- lm(freq~samePatch, family = "binomial", data = comboActBudg[comboActBudg$newState == "state4",])
summary(m.commAb)
anova(m.commAb)

# What about synchronization of behavior in a patch? ####
#Group slow flight & move into the same behavior
load("./data/12_AllBatsSamePatch.Rdata")
inPatch <- forage %>% filter(samePatch == 1)
inPatch$batStateRecode <- NA
inPatch$otherBatStateRecode <- NA
inPatch$sameBehav <- NA
for(i in 1:length(inPatch$newState)){
  inPatch$batStateRecode[i] <- ifelse(inPatch$newState[i] == "state2", "forage", 
                                      ifelse(inPatch$newState[i] == "state3", "forage", inPatch$newState[i]))
  inPatch$otherBatStateRecode[i] <- ifelse(inPatch$otherBatState[i] == "state2", "forage", 
                                          ifelse(inPatch$otherBatState[i] == "state3", "forage", inPatch$otherBatState[i]))
}

inPatch$sameBehav <- ifelse(inPatch$batStateRecode == inPatch$otherBatStateRecode,"sameBehav", "diffBehav")
inPatch$sameBehav_orig <- ifelse(inPatch$newState == inPatch$otherBatState,"sameBehav", "diffBehav")

#Calculate the frequency of synchronization for each behavior for each bat day.
sameSums <- inPatch %>% 
  dplyr::filter(sameBehav == "sameBehav") %>% 
  dplyr::group_by(groupID, batID, batIDday, newState) %>% 
  dplyr::summarize(synchStateObs = n()) %>% 
  ungroup() %>% 
  dplyr::group_by(groupID, batID, batIDday) %>% 
  dplyr::mutate(dayNobs = sum(synchStateObs), 
                freq = synchStateObs/dayNobs)

library(gghalves)
q1 <- glmer(freq~newState+(1|batID), family = "binomial", data = sameSums)
q1 <- lm(freq~newState, family = "binomial", data = sameSums)
Anova(q1)
summary(q1)


behavSameFreq <- sameSums %>% 
  ggplot()+
  geom_half_boxplot(aes(x = newState, y = freq), nudge = 0.02, outlier.color = NA) +
  geom_half_point(aes(x = newState, y = freq, color = groupID), size = 2, shape = 16, position=position_jitter(width = 0.1))+
  scale_x_discrete(breaks = c("state1", "state2", "state3"),
                   labels = c("Rest", "Slow\nflight" ,"Move"))+
  scale_fill_manual(values = mycols[2:3], 
                    breaks = c("brown", "yellow"),
                    labels = c("Group 2", "Group 3"))+
  scale_color_manual(values = mycols[2:3], 
                     breaks = c("brown", "yellow"),
                     labels = c("Group 2", "Group 3"))+
  labs(x = "", 
       y = "frequency", 
       subtitle = "Behavioral sychrony of nearest neighbors\n in a patch")+
  theme(legend.position = c(0.7, 0.9), 
        legend.title = element_blank())
behavSameFreq

# How does proximity affect rest? ####
#This will start with a behavioral change point analysis to count up behavioral sequences.
# I think I can use the patchStay df from 9. Patch Residency Influence since it split up both the in/out of a patch and the behavior in the same loop by bat day. This only was done on the 'forage' behavior class (so no commuting included). This df has duplicates b/c all bat pairs are included for distance, patch effects, etc. I'll keep this analysis under 1000 m. 

load("./data/12_AllBatsPatchStateSeq.Rdata") #patchStay_df

behavBoutTimeDist <- patchStay_df %>% filter(otherDist < 100) %>% 
  group_by(batIDday, stateBout, newState) %>% 
  dplyr::summarize(timeElapsed = difftime(max(timestamp), min(timestamp), units = "secs"),
                   meanDistNN = mean(otherDist, na.rm = TRUE),
                   minDistNN  = min(otherDist, na.rm = TRUE),
                   sdDistNN = min(otherDist),
                   medianDistNN = median(otherDist, na.rm = TRUE),
                   madDistNN = mad(otherDist, na.rm = TRUE))
stateNo <- c("state1", "state2", "state3", "state4")
behavName <- c("rest", "slow flight" ,"move", "commute")
names(behavName) <- stateNo
#At least with this figure, being in close proximity increases the duration of rest. Try to find an non-linear fit for rest
restProx <- behavBoutTimeDist %>% filter(newState == "state1")
batIDday_group <- patchStay_df %>% dplyr::select(batIDday, groupID) %>% 
  distinct()
restProx <- restProx %>% left_join(batIDday_group)

#library(remotes)
#install_github("OnofriAndreaPG/aomisc")
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
  geom_point(aes(x = minDistNN, y = as.numeric(timeElapsed)/60, color = groupID), alpha = 0.6)+
  scale_color_manual(values = mycols, 
                     name = "", 
                     breaks = c("blue", "brown", "yellow"),
                     labels = c("Group 1", "Group 2", "Group 3"))+
  geom_line(aes(x = minDistNN, y = time), data = demo.fits, color = "red")+
  labs(x = "distance (m)", 
       y = "duration of resting bout (min)", 
       subtitle = "Proximity influences on resting\nbout duration")+
  theme(legend.position = c(0.7, 0.8), 
        legend.background = element_blank(),
        legend.title = element_blank())
  

top_row <- plot_grid(actBudgOut, actBudgIn, labels = c('A', 'B'), label_size = 12, ncol = 2)
bottom_row <- plot_grid(behavSameFreq, distRestDur, labels = c('C', 'D'), label_size = 12, ncol = 2)

pdf("./output/Fig 2 - Patch & Distance Effects on Behavior.pdf", width = 8, height = 4)
plot_grid(actBudgInOut, distRestDur, labels = c('A', 'B'), label_size = 12, ncol = 2)
dev.off()



