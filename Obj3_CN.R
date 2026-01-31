#CN analyses

library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(viridis)
library(MuMIn)
library(vegan)

class(df.noUB$Horizon)

############## Data Loading and Tidying ##############################
#load tidied C/N data
df.noUB <- read.csv("C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/LAVO22_CN_metadata.csv")

df.noUB$Unit <- factor(df.noUB$Unit)
df.noUB$SampleSev <- factor(df.noUB$SampleSev, levels = c(1, 2, 3, 4))
df.noUB$Horizon <- factor(df.noUB$Horizon)

## Z... Do only once to tidy raw data for subsequent analyses
#load data
CN2 <- read.csv("C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/CN_051424.csv")
CN1 <- read.csv("C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/CN_032824.csv")

#merge two datasets from EA analysis runs 1 and 2
CN <- rbind(CN1, CN2)

#add columns removing data for samples with percent below or above detection limit
CN <- CN %>%
  mutate(cleanN = ifelse(Notes=="N below detection limit", NA, ifelse(Notes=="CN below detection limit",NA, percN))) %>%
  mutate(cleanC = ifelse(Notes=="C above detection limit", NA, ifelse(Notes=="CN below detection limit", NA, percC))) %>%
  mutate(CN.Notes = Notes) %>%
  mutate(CN.Notes = ifelse(is.na(CN.Notes), "", CN.Notes))

temp <- CN %>%
  group_by(CN.Notes) %>%
  summarise(count=n())
rm(temp)

#write.csv(CN, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/CN_merged.csv")

#add C/N concentrations to metadata
temp <- CN %>%
  select(1, 3:4, 6:8)
metadata <- read.csv("C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/03112024_LAVO22_Metadata.csv")

newmetadata <- left_join(metadata, temp, by=c("Full_Sample_Name" = "SampleName")) %>%
  mutate(CN.Notes=ifelse(is.na(CN.Notes), "", CN.Notes))

#save metadata
#write.csv(newmetadata, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/05312024_LAVO22_Metadata.csv")
rm(temp, CN1, CN2, CN, metadata)

#filter data so no duplicates and remove columns pertinent for sequencing
df <- newmetadata %>%
  filter(startsWith(Run_Submission_Name, "1")) %>% #removes 2nd sequencing run
  filter(!endsWith(Run_Submission_Name, "2")) %>% #removes duplicate samples from sequencing
  select(-(2:6), -(14:16), -(22:26)) %>% #removes extraneous columns
  mutate(Plot_Horizon=paste(Plot_Name, Horizon, sep="_")) %>% 
  mutate(unit=ifelse(Unit=="BL", "Butte Lake", ifelse(Unit=="H", "Hole", "Warner Valley"))) %>%
  mutate(FireHistory=ifelse(Trt=="LF", "less fire", "more fire"))

#filter data to remove unburned unit
df.noUB <- filter(df, Unit!="UB")

#set column types
df.noUB$Unit <- factor(df.noUB$Unit)
df.noUB$SampleSev <- factor(df.noUB$SampleSev, levels = c(1, 2, 3, 4))

#write.csv(df.noUB, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/LAVO22_CN_metadata.csv")

################ Data Summarization ####################
#sum samples for each grouping 
groupsamplecount <- df %>%
  group_by(unit, FireHistory, Horizon) %>%
  summarise(totalsamples=n())

unitsamplecount <- df %>%
  group_by(unit) %>%
  summarise(totalsamples=n())

fhsamplecount <- df %>%
  group_by(FireHistory) %>%
  summarise(totalsamples=n())

hsamplecount <- df %>%
  group_by(Horizon) %>%
  summarise(totalsamples=n())


# count samples with Notes from EA analyzer
# grouped by unit, fire history, and horizon
NotesSummary <- df %>%
  group_by(unit, FireHistory, Horizon, CN.Notes) %>%
  summarise(count=n())
NotesSummary <- left_join(NotesSummary, groupsamplecount, by=c("unit", "FireHistory", "Horizon")) 
NotesSummary <- NotesSummary %>%
  mutate(perc=(count/totalsamples)*100)

#grouped by unit
NotesSummary.unit <- df %>%
  group_by(unit, CN.Notes) %>%
  summarise(count=n())
NotesSummary.unit <- left_join(NotesSummary.unit, unitsamplecount, by="unit")
NotesSummary.unit <- mutate(NotesSummary.unit, perc=(count/totalsamples)*100)

#grouped by fire history
NotesSummary.fh <- df %>%
  group_by(FireHistory, CN.Notes) %>%
  summarise(count=n())
NotesSummary.fh <- left_join(NotesSummary.fh, fhsamplecount, by="FireHistory")
NotesSummary.fh <- mutate(NotesSummary.fh, perc=(count/totalsamples)*100)

#grouped by horizon
NotesSummary.h <- df %>%
  group_by(Horizon, CN.Notes) %>%
  summarise(count=n())
NotesSummary.h <- left_join(NotesSummary.h, hsamplecount, by="Horizon")
NotesSummary.h <- mutate(NotesSummary.h, perc=(count/totalsamples)*100)

#are there any systematic biases/patterns for N below detection, C above 
#detection, or No data based off of where the samples came from?
#visualize % of samples with EA analysis issues based on groupings
ggplot(NotesSummary.unit, aes(unit, perc)) + geom_bar(stat="identity") + facet_wrap(~CN.Notes)

ggplot(NotesSummary.fh, aes(FireHistory, perc)) + geom_bar(stat="identity") + facet_wrap(~CN.Notes)

ggplot(NotesSummary.h, aes(Horizon, perc)) + geom_bar(stat="identity") + facet_wrap(~CN.Notes)

# C above detection limit: O > M; WV > H > BL
# N below detection limit: M > O; LF > MF; BL > H & WV
# No data: O > M; H < BL & WV

rm(NotesSummary, NotesSummary.fh, NotesSummary.h, NotesSummary.unit, groupsamplecount, hsamplecount, fhsamplecount, unitsamplecount)


############## Analyses removing samples outside of detection limits ###############
################# Data Visualization ########################
#visualize C/N data by treatment faceted by unit and horizon
C.firehistory.M <- ggplot(df.noUB[df.noUB$Horizon=="M", ], aes(Trt, cleanC, fill=Trt)) 
C.firehistory.M <- C.firehistory.M + geom_boxplot() + theme_bw() + facet_wrap(~unit)
C.firehistory.M <- C.firehistory.M + scale_fill_manual(name="Fire History", values=c('#0097d3','#ffbb22'), labels = c("Less Fire", "More Fire"))
C.firehistory.M <- C.firehistory.M + scale_y_continuous(name="% Carbon")
C.firehistory.M <- C.firehistory.M + scale_x_discrete(name="Fire History", labels=c("Less Fire", "More Fire"))
C.firehistory.M

C.firehistory.O <- ggplot(df.noUB[df.noUB$Horizon=="O", ], aes(Trt, cleanC, fill=Trt)) 
C.firehistory.O <- C.firehistory.O + geom_boxplot() + theme_bw() + facet_wrap(~unit)
C.firehistory.O <- C.firehistory.O + scale_fill_manual(name="Fire History", values=c('#0097d3','#ffbb22'), labels = c("Less Fire", "More Fire"))
C.firehistory.O <- C.firehistory.O + scale_y_continuous(name="% Carbon")
C.firehistory.O <- C.firehistory.O + scale_x_discrete(name="Fire History", labels=c("Less Fire", "More Fire"))
C.firehistory.O

N.firehistory.M <- ggplot(df.noUB[df.noUB$Horizon=="M", ], aes(Trt, cleanN, fill=Trt)) 
N.firehistory.M <- N.firehistory.M + geom_boxplot() + theme_bw() + facet_wrap(~unit)
N.firehistory.M <- N.firehistory.M + scale_fill_manual(name="Fire History", values=c('#0097d3','#ffbb22'), labels = c("Less Fire", "More Fire"))
N.firehistory.M <- N.firehistory.M + scale_y_continuous(name="% Nitrogen")
N.firehistory.M <- N.firehistory.M + scale_x_discrete(name="Fire History", labels=c("Less Fire", "More Fire"))
N.firehistory.M

N.firehistory.O <- ggplot(df.noUB[df.noUB$Horizon=="O", ], aes(Trt, cleanN, fill=Trt)) 
N.firehistory.O <- N.firehistory.O + geom_boxplot() + theme_bw() + facet_wrap(~unit)
N.firehistory.O <- N.firehistory.O + scale_fill_manual(name="Fire History", values=c('#0097d3','#ffbb22'), labels = c("Less Fire", "More Fire"))
N.firehistory.O <- N.firehistory.O + scale_y_continuous(name="% Nitrogen")
N.firehistory.O <- N.firehistory.O + scale_x_discrete(name="Fire History", labels=c("Less Fire", "More Fire"))
N.firehistory.O

ggplot(df.noUB, aes(Trt, cleanC)) +geom_boxplot() + facet_grid(Horizon~unit)

ggplot(df.noUB, aes(Trt, cleanN)) +geom_boxplot() + facet_grid(Horizon~unit)



#C/N data by burn severity by unit and horizon
C.burnsev.M <- ggplot(df.noUB[df.noUB$Horizon=="M", ], aes(SampleSev, cleanC, fill=SampleSev)) 
C.burnsev.M <- C.burnsev.M + geom_boxplot() + theme_bw() + facet_wrap(~unit)
C.burnsev.M <- C.burnsev.M + scale_fill_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'),labels = c("Unburned", "Low", "Moderate", "Severe"))
C.burnsev.M <- C.burnsev.M + scale_y_continuous(name="% Carbon")
C.burnsev.M <- C.burnsev.M + scale_x_discrete(name="Burn Severity", labels=c("Unburned", "Low", "Moderate", "Severe"))
C.burnsev.M

C.burnsev.O <- ggplot(df.noUB[df.noUB$Horizon=="O", ], aes(SampleSev, cleanC, fill=SampleSev)) 
C.burnsev.O <- C.burnsev.O + geom_boxplot() + theme_bw() + facet_wrap(~unit)
C.burnsev.O <- C.burnsev.O + scale_fill_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'),labels = c("Unburned", "Low", "Moderate", "Severe"))
C.burnsev.O <- C.burnsev.O + scale_y_continuous(name="% Carbon")
C.burnsev.O <- C.burnsev.O + scale_x_discrete(name="Burn Severity", labels=c("Unburned", "Low", "Moderate", "Severe"))
C.burnsev.O

N.burnsev.M <- ggplot(df.noUB[df.noUB$Horizon=="M", ], aes(SampleSev, cleanN, fill=SampleSev)) 
N.burnsev.M <- N.burnsev.M + geom_boxplot() + theme_bw() + facet_wrap(~unit)
N.burnsev.M <- N.burnsev.M + scale_fill_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'),labels = c("Unburned", "Low", "Moderate", "Severe"))
N.burnsev.M <- N.burnsev.M + scale_y_continuous(name="% Nitrogen")
N.burnsev.M <- N.burnsev.M + scale_x_discrete(name="Burn Severity", labels=c("Unburned", "Low", "Moderate", "Severe"))
N.burnsev.M

N.burnsev.O <- ggplot(df.noUB[df.noUB$Horizon=="O", ], aes(SampleSev, cleanN, fill=SampleSev)) 
N.burnsev.O <- N.burnsev.O + geom_boxplot() + theme_bw() + facet_wrap(~unit)
N.burnsev.O <- N.burnsev.O + scale_fill_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'),labels = c("Unburned", "Low", "Moderate", "Severe"))
N.burnsev.O <- N.burnsev.O + scale_y_continuous(name="% Nitrogen")
N.burnsev.O <- N.burnsev.O + scale_x_discrete(name="Burn Severity", labels=c("Unburned", "Low", "Moderate", "Severe"))
N.burnsev.O


ggplot(df.noUB, aes(SampleSev, percN)) + geom_point() + theme_bw() + facet_grid(Horizon~unit)


################## Statistics ################################
#linear mixed model with random effects of plot and core 

########## %C Model Selection ##################
#Forwards model selection 
#Horizon
mod.hor <- lmer(percC ~ Horizon + (1|Plot_Name/Rep), data=df.noUB)
summary(mod.hor)
anova(mod.hor) #horizon is significant; p=2.2e-16

#Unit
mod.unit <- lmer(percC~Unit + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.unit) # *** p=1.6e-8

#Soil subgroup 
mod.soil <- lmer(percC~Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.soil) # *** p=0.0002097

#Sample severity 
mod.sev <-lmer(percC~SampleSev + (1|Plot_Name/Rep), data=df.noUB)
summary(mod.sev)
anova(mod.sev) #* p=0.030

#VegType
mod.veg <- lmer(percC~Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.veg) #*** p=8.1e-6

#Trt 
mod.trt <- lmer(percC~Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.trt) # Trt not sig p=0.5

#unit and soil subgroup
mod.u.sg <- lmer(percC~Unit+Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.sg) #when accounting for unit, soil type is no longer a sig predictor

#unit and veg type
mod.u.v <- lmer(percC~Unit+Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.v) #when accounting for unit, veg type is no longer a sig predictor

#unit and horizon
mod.u.h <-lmer(percC~Unit+Horizon+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h) #both unit and horizon sig predictors

#unit, horizon, trt
mod.u.h.t <- lmer(percC~Unit+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t) #trt not sig predictor

#unit, horizon, samplesev
mod.u.h.ss <- lmer(percC~Unit+Horizon+SampleSev+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss) #burn severity is sig predictor

#unit, horizon, samplesev, trt
mod.u.h.ss.t <- lmer(percC~Unit+Horizon+SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t) #trt still not sig predictor

##interactions 2-way
mod.u.h.2.ss.t <- lmer(percC~Unit*Horizon+SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.t) #2-way interaction between unit and horizon (p=0.00048)
mod.u.ss.2.h.t <- lmer(percC~Unit*SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.ss.2.h.t) #no 2-way interaction of unit and samplesev
mod.u.t.2.ss.h <- lmer(percC~Unit*Trt+SampleSev+Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.t.2.ss.h) #no Unit*Trt
mod.u.h.ss.2.t <- lmer(percC~Unit+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.2.t) #Horizon*SampleSev
mod.u.h.t.2.ss <- lmer(percC~Unit+Horizon*Trt+SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t.2.ss) #no Horizon*Trt
mod.u.h.ss.t.2 <- lmer(percC~Unit+Horizon+SampleSev*Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t.2) #no SampleSev*Trt

mod.u.h.2.ss.2.t <- lmer(percC~Unit*Horizon+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.2.t) # unit*horizon and horizon*SampleSev sig

##interactions 3-way
mod.u.h.ss.3.t <-lmer(percC~Unit*Horizon*SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.3.t) #no unit*horizon*samplesev

#from forward model selection appears that Unit*Horizon+Horizon*SampleSev+Trt is the best fit model
summary(mod.u.h.2.ss.2.t)


anova(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.2.ss.2.t, mod.u.h.ss.3.t)
#model with lowest AIC: Unit*Horizon+Horizon*SampleSev+Trt (1831.3)
#model with 2nd lowest AIC: Unit*Horizon*SampleSev+Trt (1834.3)
#3rd lowest AIC: Unit*Horizon+SampleSev+Trt (1834.7)

#%C model selected is Unit*Horizon+Horizon*SampleSev+Trt
finalmod.C <- mod.u.h.2.ss.2.t
rm(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.ss.3.t)


########## %N Model Selection #####################
#Forwards model selection 
#Horizon
mod.hor <- lmer(percN ~ Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.hor) #horizon is significant; p=2.2e-16

#Unit
mod.unit <- lmer(percN~Unit + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.unit) # *** p=2e-13

#Soil subgroup 
mod.soil <- lmer(percN~Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.soil) # *** p=0.1.38e-6

#Sample severity 
mod.sev <-lmer(percN~SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.sev) #* p=0.04812

#VegType
mod.veg <- lmer(percN~Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.veg) #*** p=2.0e-8

#Trt 
mod.trt <- lmer(percN~Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.trt) # Trt not sig p=0.608

#unit and soil subgroup
mod.u.sg <- lmer(percN~Unit+Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.sg) #when accounting for unit, soil type is no longer a sig predictor

#unit and veg type
mod.u.v <- lmer(percN~Unit+Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.v) #when accounting for unit, veg type is no longer a sig predictor

#unit and horizon
mod.u.h <-lmer(percN~Unit+Horizon+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h) #both unit and horizon sig predictors

#unit, horizon, trt
mod.u.h.t <- lmer(percN~Unit+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t) #trt not sig predictor

#unit, horizon, samplesev
mod.u.h.ss <- lmer(percN~Unit+Horizon+SampleSev+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss) #burn severity is sig predictor

#unit, horizon, samplesev, trt
mod.u.h.ss.t <- lmer(percN~Unit+Horizon+SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t) #trt still not sig predictor

##interactions 2-way
mod.u.h.2.ss.t <- lmer(percN~Unit*Horizon+SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.t) #2-way interaction between unit and horizon (p=1.88e-6)
mod.u.ss.2.h.t <- lmer(percN~Unit*SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.ss.2.h.t) #no 2-way interaction of unit and samplesev
mod.u.t.2.ss.h <- lmer(percN~Unit*Trt+SampleSev+Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.t.2.ss.h) #no Unit*Trt
mod.u.h.ss.2.t <- lmer(percN~Unit+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.2.t) #Horizon*SampleSev
mod.u.h.t.2.ss <- lmer(percN~Unit+Horizon*Trt+SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t.2.ss) #no Horizon*Trt
mod.u.h.ss.t.2 <- lmer(percN~Unit+Horizon+SampleSev*Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t.2) #no SampleSev*Trt

mod.u.h.2.ss.2.t <- lmer(percN~Unit*Horizon+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.2.t) # unit*horizon and horizon*SampleSev sig


##interactions 3-way
mod.u.h.ss.3.t <-lmer(percN~Unit*Horizon*SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.3.t) #no unit*horizon*samplesev

#from forward model selection appears that Unit*Horizon+Horizon*SampleSev+Trt is the best fit model

anova(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.2.ss.2.t, mod.u.h.ss.3.t)
#model with lowest AIC: Unit*Horizon+Horizon*SampleSev+Trt (-136.2)
#model with 2nd lowest AIC: Unit*Horizon*SampleSev+Trt (130.0)
#3rd lowest AIC: Unit*Horizon+SampleSev+Trt (-127.0)

#lowest BIC: Unit*Horizon+SampleSev+Trt (-79.3)
#2nd lowest BIC: Unit*Horizon+Horizon*SampleSev+Trt (-77.4)
#3rd lowest BIC: Unit+Horizon (-74.7)

#%N model selected is Unit*Horizon+Horizon*SampleSev+Trt
finalmod.N <- mod.u.h.2.ss.2.t
rm(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.ss.3.t)



#############Final models ################
#%C final model
finalmod.C <- mod.u.h.2.ss.2.t
aov3C <- anova(finalmod.C) 
#unit, horizon, samplesev, unit*horizon & horizon*samplesev are sig predictors of %C

#write.csv(aov3C, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Cmodel_aov3table.csv")
r.squaredGLMM(finalmod.C) #marginal (variance explained by fixed): 0.513; conditional (variance explained by both fixed and random): 0.562

#%N final model
finalmod.N <- mod.u.h.2.ss.2.t
aov3N <- anova(finalmod.N)
r.squaredGLMM(finalmod.N) #marginal: 0.549; conditional: 0.615

#write.csv(aov3N, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Nmodel_aov3table.csv")

########## Contrasts on final models ###################
############ %C contrasts ###############
###post hoc test on final model 
#fire history contrasts
con.trt <- emmeans(finalmod.C, pairwise~Trt)
con.trt
t <- as.data.frame(con.trt$contrasts) %>%
  mutate(., Unit=NA, Horizon=NA, SampleSev=NA) #save as a dataframe to merge with other contrasts

#burn severity contrasts
con.samplesevbyhorizon <-emmeans(finalmod.C, pairwise~SampleSev|Horizon)
con.samplesevbyhorizon #Mineral horizons no sig difference by burn severity
#Organic horizons: 1&2 < 3&4
ss.by.h <- as.data.frame(con.samplesevbyhorizon$contrasts) %>%
  mutate(., SampleSev=NA, Unit=NA)

###what if allowed for 3-way interaction: 
#mod.u.h.ss.3.t <- lmer(percC~Unit*Horizon*SampleSev+Trt+(1|Plot_Name/Rep), data=df.noUB)
#test <- emmeans(mod.u.h.ss.3.t, pairwise~SampleSev|Unit+Horizon)
#test #mineral horizons no sig diff by burn severity; 
#O: BL no sig diff; H: 1 >2 & 1>4, but 1=3; WV: 2>3&4

con.horbyunitsamplesev <- emmeans(mod.u.h.2.ss.2.t, pairwise~Horizon|Unit+SampleSev)
con.horbyunitsamplesev # O %C > M %C except in burn sev class 4 in Butte Lake where O & M %C are equal
h.by.ss <- as.data.frame(con.horbyunitsamplesev$contrasts) %>%
  mutate(., Horizon=NA)

con.unitbyhorizon <-emmeans(mod.u.h.2.ss.2.t, pairwise~Unit|Horizon)
con.unitbyhorizon # Mineral: BL & H << WV; Organic: BL < H < WV
u.by.h <- as.data.frame(con.unitbyhorizon$contrasts) %>%
  mutate(., Unit=NA, SampleSev=NA)

#bind together relevant contrasts
con.df <- rbind(ss.by.h, h.by.ss)
con.df <- rbind(con.df, u.by.h) 
con.df <- rbind(con.df, t)

#write.csv(con.df, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Cmodel_contrasts.csv")
rm(con.horbyunitsamplesev, con.all.df, con.trt.df, con.unit.df, con.samplesevbyunithor, con.trt, con.unitbyhorizon, t, ss.by.u.h, h.by.ss, u.by.h)

C.con.df <- con.df
rm(con.df)

############ %N Contrasts ####################
###post hoc test on final model 
con.trt <- emmeans(finalmod.N, pairwise~Trt)
con.trt
t <- as.data.frame(con.trt$contrasts) %>%
  mutate(., Unit=NA, Horizon=NA, SampleSev=NA)

con.samplesevbyhor <-emmeans(finalmod.N, pairwise~SampleSev|Horizon)
con.samplesevbyhor #Mineral horizons no sig difference by burn severity
#Organic horizons: 2 < 3&4; 1 is not different than 2, 3 or 4. This is strange.
ss.by.h <- as.data.frame(con.samplesevbyhor$contrasts) %>%
  mutate(., SampleSev=NA, Unit=NA)

####what if allowed for 3-way interaction: 
#test <- emmeans(mod.u.h.ss.3.t, pairwise~SampleSev|Unit+Horizon)
#test #mineral horizons no sig diff by burn severity; 
#O: BL no sig diff; H: no sig diff; WV 2 < 3 & 4
# what if included Unit*SampleSev interaction: 
#test <- lmer(percN~Unit*Horizon+Horizon*SampleSev+Unit*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
#test2 <- emmeans(test, pairwise~SampleSev|Unit+Horizon)
#test2 #Error in crossprod(nbasis, x) : "crossprod" is not a BUILTIN function

con.horbyunitsamplesev <- emmeans(finalmod.N, pairwise~Horizon|Unit+SampleSev)
con.horbyunitsamplesev 
# O %N > M %N 
# no difference in O & M %N: BL severity 1, BL severity 3,
#H severity 3, BL severity 4, H severity 4
h.by.ss <- as.data.frame(con.horbyunitsamplesev$contrasts) %>%
  mutate(., Horizon=NA)

con.unitbyhorizon <-emmeans(finalmod.N, pairwise~Unit|Horizon)
con.unitbyhorizon # Mineral: BL & H << WV; Organic: BL < H < WV
u.by.h <- as.data.frame(con.unitbyhorizon$contrasts) %>%
  mutate(., Unit=NA, SampleSev=NA)

#bind together relevant contrasts
con.df <- rbind(ss.by.h, h.by.ss)
con.df <- rbind(con.df, u.by.h) 
con.df <- rbind(con.df, t)

#write.csv(con.df, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Nmodel_contrasts.csv")
rm(con.horbyunitsamplesev, con.all.df, con.trt.df, con.unit.df, con.samplesevbyunithor, con.trt, con.unitbyhorizon, t, ss.by.u.h, h.by.ss, u.by.h)

N.con.df <- con.df
rm(con.df)


################ Z Analyses including all Samples (C above and N below detection)#####
################# Data Visualization ########################
#visualize C/N data by treatment faceted by unit and horizon
ggplot(df.noUB, aes(Trt, percC)) +geom_boxplot() + facet_grid(Horizon~unit)

ggplot(df.noUB, aes(Trt, percN)) +geom_boxplot() + facet_grid(Horizon~unit)


#removing the samples with C and N concentrations outside of the detection limits
#of the EA analyzer
ggplot(df.noUB, aes(Trt, cleanC)) +geom_boxplot() + facet_grid(Horizon~unit)

ggplot(df.noUB, aes(Trt, cleanN)) +geom_boxplot() + facet_grid(Horizon~unit)

#C/N data by burn severity by unit and horizon
ggplot(df.noUB, aes(SampleSev, percC)) + geom_point() + facet_grid(Horizon~unit)

ggplot(df.noUB, aes(SampleSev, percN)) + geom_point() + facet_grid(Horizon~unit)


################## Statistics ################################
#linear mixed model with random effects of plot and core 

########## %C Model Selection ##################
#Forwards model selection 
#Horizon
mod.hor <- lmer(percC ~ Horizon + (1|Plot_Name/Rep), data=df.noUB)
summary(mod.hor)
anova(mod.hor) #horizon is significant; p=2.2e-16

#Unit
mod.unit <- lmer(percC~Unit + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.unit) # *** p=1.6e-8

#Soil subgroup 
mod.soil <- lmer(percC~Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.soil) # *** p=0.0002097

#Sample severity 
mod.sev <-lmer(percC~SampleSev + (1|Plot_Name/Rep), data=df.noUB)
summary(mod.sev)
anova(mod.sev) #* p=0.030

#VegType
mod.veg <- lmer(percC~Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.veg) #*** p=8.1e-6

#Trt 
mod.trt <- lmer(percC~Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.trt) # Trt not sig p=0.5

#unit and soil subgroup
mod.u.sg <- lmer(percC~Unit+Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.sg) #when accounting for unit, soil type is no longer a sig predictor

#unit and veg type
mod.u.v <- lmer(percC~Unit+Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.v) #when accounting for unit, veg type is no longer a sig predictor

#unit and horizon
mod.u.h <-lmer(percC~Unit+Horizon+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h) #both unit and horizon sig predictors

#unit, horizon, trt
mod.u.h.t <- lmer(percC~Unit+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t) #trt not sig predictor

#unit, horizon, samplesev
mod.u.h.ss <- lmer(percC~Unit+Horizon+SampleSev+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss) #burn severity is sig predictor

#unit, horizon, samplesev, trt
mod.u.h.ss.t <- lmer(percC~Unit+Horizon+SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t) #trt still not sig predictor

##interactions 2-way
mod.u.h.2.ss.t <- lmer(percC~Unit*Horizon+SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.t) #2-way interaction between unit and horizon (p=0.00048)
mod.u.ss.2.h.t <- lmer(percC~Unit*SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.ss.2.h.t) #no 2-way interaction of unit and samplesev
mod.u.t.2.ss.h <- lmer(percC~Unit*Trt+SampleSev+Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.t.2.ss.h) #no Unit*Trt
mod.u.h.ss.2.t <- lmer(percC~Unit+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.2.t) #Horizon*SampleSev
mod.u.h.t.2.ss <- lmer(percC~Unit+Horizon*Trt+SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t.2.ss) #no Horizon*Trt
mod.u.h.ss.t.2 <- lmer(percC~Unit+Horizon+SampleSev*Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t.2) #no SampleSev*Trt

mod.u.h.2.ss.2.t <- lmer(percC~Unit*Horizon+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.2.t) # unit*horizon and horizon*SampleSev sig

##interactions 3-way
mod.u.h.ss.3.t <-lmer(percC~Unit*Horizon*SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.3.t) #no unit*horizon*samplesev

#from forward model selection appears that Unit*Horizon+Horizon*SampleSev+Trt is the best fit model
summary(mod.u.h.2.ss.2.t)


anova(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.2.ss.2.t, mod.u.h.ss.3.t)
#model with lowest AIC: Unit*Horizon+Horizon*SampleSev+Trt (1831.3)
#model with 2nd lowest AIC: Unit*Horizon*SampleSev+Trt (1834.3)
#3rd lowest AIC: Unit*Horizon+SampleSev+Trt (1834.7)

#%C model selected is Unit*Horizon+Horizon*SampleSev+Trt
finalmod.C <- mod.u.h.2.ss.2.t
rm(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.ss.3.t)


########## %N Model Selection #####################
#Forwards model selection 
#Horizon
mod.hor <- lmer(percN ~ Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.hor) #horizon is significant; p=2.2e-16

#Unit
mod.unit <- lmer(percN~Unit + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.unit) # *** p=2e-13

#Soil subgroup 
mod.soil <- lmer(percN~Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.soil) # *** p=0.1.38e-6

#Sample severity 
mod.sev <-lmer(percN~SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.sev) #* p=0.04812

#VegType
mod.veg <- lmer(percN~Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.veg) #*** p=2.0e-8

#Trt 
mod.trt <- lmer(percN~Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.trt) # Trt not sig p=0.608

#unit and soil subgroup
mod.u.sg <- lmer(percN~Unit+Soil_Subgroup + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.sg) #when accounting for unit, soil type is no longer a sig predictor

#unit and veg type
mod.u.v <- lmer(percN~Unit+Veg_Type + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.v) #when accounting for unit, veg type is no longer a sig predictor

#unit and horizon
mod.u.h <-lmer(percN~Unit+Horizon+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h) #both unit and horizon sig predictors

#unit, horizon, trt
mod.u.h.t <- lmer(percN~Unit+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t) #trt not sig predictor

#unit, horizon, samplesev
mod.u.h.ss <- lmer(percN~Unit+Horizon+SampleSev+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss) #burn severity is sig predictor

#unit, horizon, samplesev, trt
mod.u.h.ss.t <- lmer(percN~Unit+Horizon+SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t) #trt still not sig predictor

##interactions 2-way
mod.u.h.2.ss.t <- lmer(percN~Unit*Horizon+SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.t) #2-way interaction between unit and horizon (p=1.88e-6)
mod.u.ss.2.h.t <- lmer(percN~Unit*SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.ss.2.h.t) #no 2-way interaction of unit and samplesev
mod.u.t.2.ss.h <- lmer(percN~Unit*Trt+SampleSev+Horizon + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.t.2.ss.h) #no Unit*Trt
mod.u.h.ss.2.t <- lmer(percN~Unit+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.2.t) #Horizon*SampleSev
mod.u.h.t.2.ss <- lmer(percN~Unit+Horizon*Trt+SampleSev + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.t.2.ss) #no Horizon*Trt
mod.u.h.ss.t.2 <- lmer(percN~Unit+Horizon+SampleSev*Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.t.2) #no SampleSev*Trt

mod.u.h.2.ss.2.t <- lmer(percN~Unit*Horizon+Horizon*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.2.ss.2.t) # unit*horizon and horizon*SampleSev sig


##interactions 3-way
mod.u.h.ss.3.t <-lmer(percN~Unit*Horizon*SampleSev+Trt+ (1|Plot_Name/Rep), data=df.noUB)
anova(mod.u.h.ss.3.t) #no unit*horizon*samplesev

#from forward model selection appears that Unit*Horizon+Horizon*SampleSev+Trt is the best fit model

anova(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.2.ss.2.t, mod.u.h.ss.3.t)
#model with lowest AIC: Unit*Horizon+Horizon*SampleSev+Trt (-136.2)
#model with 2nd lowest AIC: Unit*Horizon*SampleSev+Trt (130.0)
#3rd lowest AIC: Unit*Horizon+SampleSev+Trt (-127.0)

#lowest BIC: Unit*Horizon+SampleSev+Trt (-79.3)
#2nd lowest BIC: Unit*Horizon+Horizon*SampleSev+Trt (-77.4)
#3rd lowest BIC: Unit+Horizon (-74.7)

#%N model selected is Unit*Horizon+Horizon*SampleSev+Trt
finalmod.N <- mod.u.h.2.ss.2.t
rm(mod.unit, mod.hor, mod.soil, mod.sev, mod.veg, mod.trt, mod.u.sg, mod.u.v, mod.u.h, mod.u.h.t, mod.u.h.ss, mod.u.h.ss.t, mod.u.h.2.ss.t, mod.u.ss.2.h.t, mod.u.t.2.ss.h, mod.u.h.ss.2.t, mod.u.h.t.2.ss, mod.u.h.ss.t.2, mod.u.h.ss.3.t)



#############Final models ################
#%C final model
finalmod.C <- mod.u.h.2.ss.2.t
aov3C <- anova(finalmod.C) 
#unit, horizon, samplesev, unit*horizon & horizon*samplesev are sig predictors of %C

#write.csv(aov3C, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Cmodel_aov3table.csv")
r.squaredGLMM(finalmod.C) #marginal (variance explained by fixed): 0.513; conditional (variance explained by both fixed and random): 0.562

#%N final model
finalmod.N <- mod.u.h.2.ss.2.t
aov3N <- anova(finalmod.N)
r.squaredGLMM(finalmod.N) #marginal: 0.549; conditional: 0.615

#write.csv(aov3N, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Nmodel_aov3table.csv")

########## Contrasts on final models ###################
############ %C contrasts ###############
###post hoc test on final model 
#fire history contrasts
con.trt <- emmeans(finalmod.C, pairwise~Trt)
con.trt
t <- as.data.frame(con.trt$contrasts) %>%
  mutate(., Unit=NA, Horizon=NA, SampleSev=NA) #save as a dataframe to merge with other contrasts

#burn severity contrasts
con.samplesevbyhorizon <-emmeans(finalmod.C, pairwise~SampleSev|Horizon)
con.samplesevbyhorizon #Mineral horizons no sig difference by burn severity
#Organic horizons: 1&2 < 3&4
ss.by.h <- as.data.frame(con.samplesevbyhorizon$contrasts) %>%
  mutate(., SampleSev=NA, Unit=NA)

###what if allowed for 3-way interaction: 
#mod.u.h.ss.3.t <- lmer(percC~Unit*Horizon*SampleSev+Trt+(1|Plot_Name/Rep), data=df.noUB)
#test <- emmeans(mod.u.h.ss.3.t, pairwise~SampleSev|Unit+Horizon)
#test #mineral horizons no sig diff by burn severity; 
#O: BL no sig diff; H: 1 >2 & 1>4, but 1=3; WV: 2>3&4

con.horbyunitsamplesev <- emmeans(mod.u.h.2.ss.2.t, pairwise~Horizon|Unit+SampleSev)
con.horbyunitsamplesev # O %C > M %C except in burn sev class 4 in Butte Lake where O & M %C are equal
h.by.ss <- as.data.frame(con.horbyunitsamplesev$contrasts) %>%
  mutate(., Horizon=NA)

con.unitbyhorizon <-emmeans(mod.u.h.2.ss.2.t, pairwise~Unit|Horizon)
con.unitbyhorizon # Mineral: BL & H << WV; Organic: BL < H < WV
u.by.h <- as.data.frame(con.unitbyhorizon$contrasts) %>%
  mutate(., Unit=NA, SampleSev=NA)

#bind together relevant contrasts
con.df <- rbind(ss.by.h, h.by.ss)
con.df <- rbind(con.df, u.by.h) 
con.df <- rbind(con.df, t)

#write.csv(con.df, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Cmodel_contrasts.csv")
rm(con.horbyunitsamplesev, con.all.df, con.trt.df, con.unit.df, con.samplesevbyunithor, con.trt, con.unitbyhorizon, t, ss.by.u.h, h.by.ss, u.by.h)

C.con.df <- con.df
rm(con.df)

############ %N Contrasts ####################
###post hoc test on final model 
con.trt <- emmeans(finalmod.N, pairwise~Trt)
con.trt
t <- as.data.frame(con.trt$contrasts) %>%
  mutate(., Unit=NA, Horizon=NA, SampleSev=NA)

con.samplesevbyhor <-emmeans(finalmod.N, pairwise~SampleSev|Horizon)
con.samplesevbyhor #Mineral horizons no sig difference by burn severity
#Organic horizons: 2 < 3&4; 1 is not different than 2, 3 or 4. This is strange.
ss.by.h <- as.data.frame(con.samplesevbyhor$contrasts) %>%
  mutate(., SampleSev=NA, Unit=NA)

####what if allowed for 3-way interaction: 
#test <- emmeans(mod.u.h.ss.3.t, pairwise~SampleSev|Unit+Horizon)
#test #mineral horizons no sig diff by burn severity; 
#O: BL no sig diff; H: no sig diff; WV 2 < 3 & 4
# what if included Unit*SampleSev interaction: 
#test <- lmer(percN~Unit*Horizon+Horizon*SampleSev+Unit*SampleSev+Trt + (1|Plot_Name/Rep), data=df.noUB)
#test2 <- emmeans(test, pairwise~SampleSev|Unit+Horizon)
#test2 #Error in crossprod(nbasis, x) : "crossprod" is not a BUILTIN function

con.horbyunitsamplesev <- emmeans(finalmod.N, pairwise~Horizon|Unit+SampleSev)
con.horbyunitsamplesev 
# O %N > M %N 
# no difference in O & M %N: BL severity 1, BL severity 3,
#H severity 3, BL severity 4, H severity 4
h.by.ss <- as.data.frame(con.horbyunitsamplesev$contrasts) %>%
  mutate(., Horizon=NA)

con.unitbyhorizon <-emmeans(finalmod.N, pairwise~Unit|Horizon)
con.unitbyhorizon # Mineral: BL & H << WV; Organic: BL < H < WV
u.by.h <- as.data.frame(con.unitbyhorizon$contrasts) %>%
  mutate(., Unit=NA, SampleSev=NA)

#bind together relevant contrasts
con.df <- rbind(ss.by.h, h.by.ss)
con.df <- rbind(con.df, u.by.h) 
con.df <- rbind(con.df, t)

#write.csv(con.df, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/CarbonNitrogen/Nmodel_contrasts.csv")
rm(con.horbyunitsamplesev, con.all.df, con.trt.df, con.unit.df, con.samplesevbyunithor, con.trt, con.unitbyhorizon, t, ss.by.u.h, h.by.ss, u.by.h)

N.con.df <- con.df
rm(con.df)


####### Additional graphs from model building ############
#visualize 2-way interactions
emmip(mod.u.h.2.ss.2.t, Unit ~ Horizon) #estimated mean %C x axis Horizon colored by Unit
emmip(mod.u.h.2.ss.2.t, Horizon ~ SampleSev) #biggest change in predicted %C occurred from low to moderate; particularly in organic horizon



