#Objective: model (ANOVA) and visualize how burn severity and fire history are 
#correlated with pH

install.packages("phyloseq", type = "source")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
library(emmeans)
library(viridis)
library(MuMIn)
library(vegan)
citation("lmerTest")
citation("vegan")

#load metadata
df <- read.csv("~/Desktop/03112024_LAVO22_Metadata.csv")
df <- read.csv("C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/03112024_LAVO22_Metadata.csv")
pH.df <- df %>%
  filter(startsWith(Run_Submission_Name, "1")) %>%
  filter(!endsWith(Run_Submission_Name, "2")) %>%
  select(-(2:6), -(14:16), -(22:26)) %>%
  mutate(Plot_Horizon=paste(Plot_Name, Horizon, sep="_")) %>%
  mutate(unit=ifelse(Unit=="BL", "Butte Lake", ifelse(Unit=="H", "Hole", "Warner Valley"))) %>%
  mutate(FireHistory=ifelse(Trt=="LF", "Less Fire", "More Fire"))

pH.df.noUB <- filter(pH.df, Unit!="UB")

pH.df.noUB$Unit <- factor(pH.df.noUB$Unit)
pH.df.noUB$SampleSev <- factor(pH.df.noUB$SampleSev, levels = c(1, 2, 3, 4))

write.csv(pH.df.noUB, "~/Desktop/LAVO22_pH_metadata.csv")

#### paired wilcoxon test for horizon #####
# make dataframe containing only samples that have both O & M horizons
pH.pair.hor <- filter(pH.df.noUB, pH!="NA") %>%
  .[.$Sample_Core_Name %in% .$Sample_Core_Name[duplicated(.$Sample_Core_Name)], ] %>%
  mutate(., Horizon=ifelse(Horizon=="O", "Organic", "Mineral")) %>%
  mutate(., unit=ifelse(Unit=="BL", "Butte Lake", ifelse(Unit=="H", "Hole", "Warner Valley")))
  

#test for normality 
hor.diff <- with(pH.pair.hor, 
                 pH[Horizon=="M"] - pH[Horizon=="O"]) 
d <- as.data.frame(hor.diff)
shapiro.test(hor.diff) #p=0.02; distribution of pH differences between O & M horizons is not normally distributed

ggplot(d, aes(hor.diff)) + geom_histogram()


# paired wilcoxon test
hor.wilcox <- wilcox.test(pH ~ Horizon, data=pH.pair.hor, paired=TRUE)
hor.wilcox # p=0.91 -> median pH in O horizon is not significantly different from the median pH in M horizon

# plot paired data
ggplot(pH.pair.hor, aes(x = Horizon, y = pH)) + 
  geom_boxplot(aes()) +
  geom_line(aes(group = Sample_Core_Name)) + 
  geom_point(size = 0) + 
  facet_wrap(~ unit)


########## Visualizing pH ##############
#is pH normally distributed? 
shapiro.test(pH.df.noUB$pH) #p<0.0000001 -> pH is not normally distributed!
#stat consultant: you should be testing whether the residuals are normally distributed
#LMMs are robust - some deviation from normality may still be ok 


###Full dataset (no UB)
#matching model
pH.df.noUB$offset <- as.numeric(factor(pH.df.noUB$Horizon)) + ifelse()
pH.df.noUB$FireHistory
pH.model <- ggplot(data=pH.df.noUB, aes(x=Horizon, y=pH, fill=SampleSev)) + geom_boxplot(position="dodge") 
pH.model <- pH.model + theme_bw() + scale_fill_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'),labels = c("unburned", "low", "moderate", "severe"))
pH.model <- pH.model + facet_grid(FireHistory~unit) 
pH.model <- pH.model + labs() + scale_x_discrete(labels=c('Mineral', 'Organic'))
pH.model

#pH by soil subgroup
pH.df.noUB$Soil_Subgroup = factor(pH.df.noUB$Soil_Subgroup, levels=c("Humic Haploxerands", "Typic Haploxerands-Lithic Haploxerands", "Humic Vitrixerands", "Typic Vitrixerands", "Vitrixerands", "Typic Xerorthent", "Vitrandic Xerofluvents"))

pH.bysoil <- ggplot(data=pH.df.noUB, aes(x=Soil_Subgroup, y=pH, color=SampleSev)) + geom_point(position=position_dodge(width=0.7)) 
pH.bysoil <- pH.bysoil + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
pH.bysoil <- pH.bysoil  + labs(x="Soil Subgroup") + facet_grid(Horizon~unit)
pH.bysoil

#pH by soil greatgroup
pH.bysoil <- ggplot(data=pH.df.noUB, aes(x=Soil_GreatGroup, y=pH, fill=FireHistory)) + geom_boxplot() 
pH.bysoil <- pH.bysoil + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
pH.bysoil <- pH.bysoil  + labs(x="Soil Greatgroup") + facet_grid(unit~Horizon) + scale_fill_manual(values=c('#0097d3','#ffbb22')) 
pH.bysoil

#pH by vegetation type
pH.byveg <- ggplot(data=pH.df.noUB, aes(x=Veg_Type, y=pH)) + geom_boxplot() 
pH.byveg <- pH.byveg + theme_bw() 
pH.byveg <- pH.byveg  + labs(x="Vegetation Type") + facet_grid(unit~Horizon)
pH.byveg

#pH by unit
pH.byunit <- ggplot(data=pH.df.noUB, aes(x=Unit, y=pH)) + geom_boxplot() 
pH.byunit <- pH.byunit + theme_bw() + scale_x_discrete(labels=c("Butte Lake", "Hole", "Warner Valley"))
pH.byunit <- pH.byunit  + labs(x="Unit")
pH.byunit

#pH by Sample burn severity 
pH.bybs <- ggplot(data=pH.df.noUB, aes(x=SampleSev, y=pH)) + geom_boxplot() 
pH.bybs <- pH.bybs + theme_bw() + scale_x_discrete(labels=c("unburned", "low", "moderate", "severe"))
pH.bybs <- pH.bybs  + labs(x="Sample Burn Severity")
pH.bybs

#pH by Fire history 
pH.bytrt <- ggplot(data=pH.df.noUB, aes(x=Trt, y=pH)) + geom_boxplot() 
pH.bytrt <- pH.bytrt + theme_bw() + scale_x_discrete(labels=c("less fire", "more fire"))
pH.bytrt <- pH.bytrt  + labs(x="Fire History")
pH.bytrt

#pH by Fire History Facet by unit and horizon 
pH.full <- ggplot(data=pH.df.noUB, aes(x=Horizon, y=pH, fill=FireHistory)) + geom_boxplot(position="dodge")
pH.full <- pH.full + facet_wrap(~unit) + scale_x_discrete(labels=c("Mineral", "Organic"))
pH.full <- pH.full + theme_bw() + scale_fill_manual(values=c('#0097d3','#ffbb22')) 
pH.full 


#pH by Treatment
ggplot(data=pH.df[pH.df$Unit!="UB", ], aes(x=Horizon, y=pH, fill=Trt)) + geom_boxplot() 

#pH by Sample Severity
ggplot(data=pH.df.noUB, aes(x=SampleSev, y=pH)) + geom_point() + facet_wrap(~Unit)

###Subset by Unit
#pH by soil subgroup
a <- ggplot(data=pH.BL, aes(x=Soil_Subgroup, y=pH)) + geom_boxplot() + labs(title="Butte Lake")
b <- ggplot(data=pH.H, aes(x=Soil_Subgroup, y=pH)) + geom_boxplot() + labs(title="Hole")
c <- ggplot(data=pH.WV, aes(x=Soil_Subgroup, y=pH)) + geom_boxplot() + labs(title="Warner Valley")
a + b + c


######## Linear Mixed Model, Full Dataset, #########
#this model accounts for Nestedness of the experimental design by treating the 
#plot as a random effect (addresses subsampling error) and core as a random 
#effect (addresses the subsampling error associated with the horizons)

#Backwards model selection 

#All together, with Unit interactions#
modAllInt <- lmer(pH ~ Unit*Trt*Horizon*Soil_Subgroup + (1|Plot_Name/Rep), data = pH.df.noUB)
summary(modAllInt)
anova(modAllInt) #only Trt is significant with a marginal 3-way interaction of Unit:Trt:Horizon

#remove interactions 
modnoInt <- lmer(pH ~ Unit+Trt+Horizon+Soil_Subgroup + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(modnoInt)
anova(modnoInt) #still only treatment is significant

#remove horizon
mod.3 <- lmer(pH ~ Unit + Trt + Soil_Subgroup + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(mod.3) 
anova(mod.3) # still only treatment is significant

#remove soil 
mod.2 <- lmer(pH ~ Unit + Trt + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(mod.2)
anova(mod.2) # now unit and treatment are significant



#Forwards model selection 

#Horizon
mod.hor <- lmer(pH ~ Horizon + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(mod.hor)
anova(mod.hor) #not significant

#Unit
mod.unit <- lmer(pH~Unit + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.unit) # ** p=0.0021

#Soil subgroup 
mod.soil <- lmer(pH~Soil_Subgroup + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.soil) # ** p=0.004

#Sample severity 
mod.sev <-lmer(pH~SampleSev + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(mod.sev)
anova(mod.sev) #***

#VegType
mod.veg <- lmer(pH~Veg_Type + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.veg) #**

#Trt 
mod.trt <- lmer(pH~Trt + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.trt) # **

#Unit and Samplesev
mod.u.s <- lmer(pH~Unit+SampleSev + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s) #unit *, samplesev ***

#Unit SampleSev and trt
mod.u.s.t <- lmer(pH~Unit+SampleSev+Trt+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.t) #unit *, samplesev **, trt *

#Unit Samplesev trt soil subgroup
mod.u.s.t.s <-lmer(pH~Unit+SampleSev+Trt+Soil_Subgroup+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.t.s) #unit and soil subgroup are correlated, unit not sig, soil subgroup not sig

#unit samplesev trt veg
mod.u.s.t.v <- lmer(pH~Unit+SampleSev+Trt+Veg_Type+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.t.v) #unit and veg type are correlated, unit marginally sig, veg type not

##Unit, Treatment, Severity, Horizon 
mod.u.s.h.t <- lmer(pH~Unit+SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=pH.df.noUB)
summary(mod.u.s.h.t)
anova(mod.u.s.h.t) 

##interactions 2-way
mod.u.s.2.h.t <- lmer(pH~Unit*SampleSev+Horizon+Trt + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.2.h.t) #no two-way interaction of unit and sample severity
mod.u.h.2.s.t <-lmer(pH~Unit*Horizon+Trt+SampleSev + (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.2.h.t) #no 2-way int of unit and horizon
mod.u.t.2.s.h <-lmer(pH~Unit*Trt+SampleSev+Horizon+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.t.2.s.h) #no 2 way int of unit and treatment
mod.u.s.h.2.t <-lmer(pH~Unit+SampleSev*Horizon+Trt+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.h.2.t) #no 2 way int of severity and horizon
mod.u.s.t.2.h <-lmer(pH~Unit+SampleSev*Trt+Horizon+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.t.2.h) #no 2 way int of severity and treatment
mod.u.s.t.h.2<-lmer(pH~Unit+SampleSev+Trt*Horizon+ (1|Plot_Name/Rep), data=pH.df.noUB)
anova(mod.u.s.t.h.2) #no 2 way int of treatment and horizon

anova(mod.u.s.h.t, mod.hor, mod.unit, mod.soil, mod.sev, mod.veg, mod.u.s.t, mod.u.s.t.s, mod.u.s.t.v, mod.u.h.2.s.t, mod.u.s.2.h.t, mod.u.s.h.2.t, mod.u.s.t.2.h, mod.u.s.t.h.2, mod.u.t.2.s.h)
#mod.u.s.h.t is selected as the final model.


#conditional and marginal R2 values
r.squaredGLMM(mod.u.s.h.t)


###post hoc test on final model 
con.all <- emmeans(mod.u.s.h.t, pairwise ~ Unit+SampleSev+Horizon+Trt, adjust="tukey")
con.all.df <- as.data.frame(con.all$contrasts)

con.test <- emmeans(mod.u.t.2.s.h, pairwise ~Unit+Trt, adjust='tukey')
con.test.df <- as.data.frame(con.test$contrasts)

con.unit <- emmeans(mod.u.s.h.t, pairwise ~ Unit, adjust="tukey")
con.ind.df <- as.data.frame(con.unit$contrasts)

con.horizon <- emmeans(mod.u.s.h.t, pairwise~Horizon, adjust="tukey")
con.ind.df <- rbind(con.ind.df, as.data.frame(con.horizon$contrasts))

con.sev <-emmeans(mod.u.s.h.t, pairwise~SampleSev, adjust="tukey")
con.ind.df <- rbind(con.ind.df, as.data.frame(con.sev$contrasts))

con.trt <-emmeans(mod.u.s.h.t, pairwise~Trt, adjust="tukey")
con.ind.df <- rbind(con.ind.df, as.data.frame(con.trt$contrast))

con.trtbyunit <- emmeans(mod.u.s.h.t, pairwise~Trt|Unit)
con.trtbyunit

write.csv(con.ind.df, "~/Desktop/LAVO22_GG2/individualgroupcontrasts.csv")
write.csv(con.all.df, "C:/Users/iemus/Box/MSWhitman/LAVO_FireHistory/Data/R/pH/pH_allcontrasts.csv")

#visualize 
emmip(mod.u.s.h.t, Trt ~ Horizon) #estimated mean pH x axis Horizon colored by Trt
emmip(mod.u.s.h.t, Horizon ~ SampleSev) #shows no difference in horizon and predicted
#pH is most different in the step from unburned to low severity, increases from low to mod
# likely significantly and then from moderate to severe is small increase
emmip(mod.u.s.h.t, Unit ~ SampleSev) #Hole is lower pH overall than other 2 units
#

################ Z Standardized effect sizes is meaningless ##############
#(meaningless to statisticians)
#standardized effect sizes for each variable - Cohen's d
#https://rdrr.io/cran/emmeans/man/eff_size.html
### Mixed model example:
if (require(nlme)) withAutoprint({
  Oats.lme <- lme(yield ~ Variety + factor(nitro), 
                  random = ~ 1 | Block / Variety,
                  data = Oats)
  # Combine variance estimates
  VarCorr(Oats.lme)
  (totSD <- sqrt(214.4724 + 109.6931 + 162.5590))
  # I figure edf is somewhere between 5 (Blocks df) and 51 (Resid df)
  emmV <- emmeans(Oats.lme, ~ Variety)
  eff_size(emmV, sigma = totSD, edf = 5)
  eff_size(emmV, sigma = totSD, edf = 51)
}, spaced = TRUE)

df.residual(mod.u.s.h.t)
sigma(mod.u.s.h.t)
eff_size(mod.u.s.h.t, sigma=sigma(mod.u.s.h.t), edf=df.residual(mod.u.s.h.t))



######## Linear Mixed Model, each unit separately##########
pH.BL <- pH.df.noUB[pH.df.noUB$Unit=="BL", ]
pH.H <- pH.df.noUB[pH.df.noUB$Unit=="H", ]
pH.WV <- pH.df.noUB[pH.df.noUB$Unit=="WV", ]

########### Butte Lake ####
#stepwise selection BL
#soil subgroup 
BL.soil <- lmer(pH ~ Soil_Subgroup + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.soil) # *
#visualise
BL.bysoil <- ggplot(data=pH.BL, aes(x=Soil_Subgroup, y=pH)) + geom_boxplot() 
BL.bysoil <- BL.bysoil + theme_bw() #+ theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
BL.bysoil <- BL.bysoil  + labs(x="Soil Subgroup") +labs(title ="Butte Lake pH by Soil Subgroup")
BL.bysoil

#trt
BL.trt <-lmer(pH ~ Trt + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.trt) # **

#sampleseverity 
BL.trt <-lmer(pH ~ SampleSev + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.trt) # **

#veg type
BL.veg <-lmer(pH ~ Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.veg) # *
BL.byveg <- ggplot(data=pH.BL, aes(x=Veg_Type, y=pH)) + geom_boxplot() 
BL.byveg <- BL.byveg + theme_bw() #+ theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
BL.byveg <- BL.byveg  + labs(x="Vegetation Type") +labs(title ="Butte Lake pH by Vegetation Type")
BL.byveg

#horizon
BL.hor <-lmer(pH ~ Horizon + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.hor) # *

##trt+sev
Bl.t.s <- lmer(pH ~ Trt+SampleSev + (1|Plot_Name/Rep), data = pH.BL)
anova(Bl.t.s) #trt & severity are correlated so trt not sig, samplesev is marginally

##sev+soil
BL.s.s <-lmer(pH~SampleSev+Soil_Subgroup+(1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.s) #sev **, soil .

##sev+soil+veg
BL.s.s.v <-lmer(pH~SampleSev+Soil_Subgroup+Veg_Type+(1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.s.v) #sev **, no other sig 

##sev+veg
BL.s.v <- lmer(pH ~ SampleSev+Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.v) #sev **, veg*

##two way interactions
BL.s.v.2 <- lmer(pH~ SampleSev*Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.v.2) #only sev is sig no 2 way effect

anova(BL.hor, BL.trt, Bl.t.s, BL.veg, BL.soil, BL.s.v.2, BL.s.v, BL.s.s.v, BL.s.s)
#lowest AIC: BL.s.s, 2nd lowest BL.s.v, 3rd lowest BL.s.s.v, 4th lowest BL.trt


########### Hole ####
#stepwise selection H
#soil subgroup 
H.soil <- lmer(pH ~ Soil_Subgroup + (1|Plot_Name/Rep), data = pH.H)
anova(H.soil) # not sig
#visualise
H.bysoil <- ggplot(data=pH.H, aes(x=Soil_Subgroup, y=pH, fill=SampleSev)) + geom_boxplot() 
H.bysoil <- H.bysoil + theme_bw() #+ theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
H.bysoil <- H.bysoil  + labs(x="Soil Subgroup") +labs(title ="Butte Lake pH by Soil Subgroup")
H.bysoil

#trt
BL.trt <-lmer(pH ~ Trt + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.trt) # **

#sampleseverity 
BL.trt <-lmer(pH ~ SampleSev + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.trt) # **

#veg type
BL.veg <-lmer(pH ~ Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.veg) # *
BL.byveg <- ggplot(data=pH.BL, aes(x=Veg_Type, y=pH)) + geom_boxplot() 
BL.byveg <- BL.byveg + theme_bw() #+ theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
BL.byveg <- BL.byveg  + labs(x="Vegetation Type") +labs(title ="Butte Lake pH by Vegetation Type")
BL.byveg

#horizon
BL.hor <-lmer(pH ~ Horizon + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.hor) # *

##trt+sev
Bl.t.s <- lmer(pH ~ Trt+SampleSev + (1|Plot_Name/Rep), data = pH.BL)
anova(Bl.t.s) #trt & severity are correlated so trt not sig, samplesev is marginally

##sev+soil
BL.s.s <-lmer(pH~SampleSev+Soil_Subgroup+(1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.s) #sev **, soil .

##sev+soil+veg
BL.s.s.v <-lmer(pH~SampleSev+Soil_Subgroup+Veg_Type+(1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.s.v) #sev **, no other sig 

##sev+veg
BL.s.v <- lmer(pH ~ SampleSev+Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.v) #sev **, veg*

##two way interactions
BL.s.v.2 <- lmer(pH~ SampleSev*Veg_Type + (1|Plot_Name/Rep), data = pH.BL)
anova(BL.s.v.2) #only sev is sig no 2 way effect

anova(BL.hor, BL.trt, Bl.t.s, BL.veg, BL.soil, BL.s.v.2, BL.s.v, BL.s.s.v, BL.s.s)
#lowest AIC: BL.s.s, 2nd lowest BL.s.v, 3rd lowest BL.s.s.v, 4th lowest BL.trt


################# Z WRONG TREATMENT OF DATA STRUCTURE ##########################
########################TREATING SAMPLES AS INDEPENDENT#######################
#visualize pH by horizon 
ggplot(data=pH.df, aes(x=Horizon, y=pH, fill=Trt)) + geom_boxplot(position="dodge") + facet_wrap(~Unit) + scale_color_manual(values=c('#56B4E9', '#E69F00'))

#interestingly there does not appear to be a dramatic difference in pH by horizon
#we would expect that horizons high in organic matter would have a lower pH than 
#the mineral horizon.

### separate pH by horizon 
pH.O <- filter(pH.df, Horizon=="O")
pH.M <- filter(pH.df, Horizon=="M")


ggplot(data=pH.df, aes(x=Unit, y=pH)) + geom_boxplot() + facet_wrap(~Horizon)
summary(aov(pH~Unit*Horizon, data=pH.df)) 
TukeyHSD(aov(pH~Unit, data=pH.df)) #units are different from each other except 
#WV-BL & UB-Hole

pH.BL <- filter(pH.df, Unit=="BL")
pH.H <- filter(pH.df, Unit=="H")
pH.WV <- filter(pH.df, Unit=="WV")

###is fire history a predictor of pH? 
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup, data=pH.BL)) #horizon is marginally significant (p<0.1) 
#*** fire history & sample burn severity are significant p<0.0001, soilsubgroup 
#** is no interaction of treatment and horizon 
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup, data=pH.H)) #horizon is not sig (p>0.1), fire history
#is marginally significant (p<0.1), no interaction of treatment and horizon
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup, data=pH.WV)) #horizon is not sig (p>0.1), fire history 
#is marginally significant (p<0.1), no interaction of treatment and horizon


#since horizon is not significant (p>0.05) remove as a variable and test for 
#fire history effects
d <- aov(pH~Trt, data=pH.BL) # *** p<0.0001
summary(aov(pH~Trt, data=pH.H)) # . p<0.1
summary(aov(pH~Trt, data=pH.WV)) # . p<0.1

###is burn severity a predictor of pH? 
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup*Veg_Type*Unit, data=pH.df.noUB)) 
# Trt SampleSev, and Soil subgroup are sig predictors, Trt:Veg_Type interaction (p<0.05)

###final model 
levels(pH.df[pH.df$Unit!="UB", ]$Unit)
pH.df$Unit <- as.factor(pH.df$Unit)
summary(aov(pH~Trt+SampleSev+Soil_Subgroup, data=pH.df[pH.df$Unit!="UB", ]))
TukeyHSD(aov(pH~Trt+SampleSev+Soil_Subgroup, data=pH.df[pH.df$Unit!="UB", ]))

#Are samples within a plot more similar to each other than to other plots?
subsample.M <- ggplot(data=pH.df.noUB[pH.df.noUB$Horizon=="M", ], aes(x=Plot_Name, y=pH)) + geom_point() + labs(title="Subsample Variation in pH by Plot - Mineral")
subsample.O <- ggplot(data=pH.df.noUB[pH.df.noUB$Horizon=="O", ], aes(x=Plot_Name, y=pH)) + geom_point() + labs(title="Subsample Variation in pH by Plot - Organic")
subsample.M
subsample.O

subsample.range <- pH.df.noUB %>%
  group_by(Horizon, Plot_Horizon) %>%
  summarise(min=min(pH), max=max(pH), range=max-min) %>%
  filter(range!="NA")

subsample.range.sum <- subsample.range %>%
  group_by(Horizon) %>%
  summarise(mean.range = mean(range), sd.range=sd(range))


#################### Select one Rep (no for loop) ###############################
##A reps
A.M.df <- filter(pH.df.noUB, Horizon=="M" & Rep=="A")
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup*Veg_Type*Unit, data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A"), ])) 
#significant variables *** Trt, SampleSev * Soil_Subgrou:Veg_Type

#visualize A reps
#all
ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A"), ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="A reps")
#by horizons 
a.M <- ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A") & pH.df.noUB$Horizon=="M", ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="A reps Mineral")
a.O <- ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A") & pH.df.noUB$Horizon=="O", ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="A reps Organic")
a.M + a.O
  
  
##B reps
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup*Veg_Type*Unit, data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "B"), ])) 
#significant variables *** SampleSev, **Soil_Subgroup, * Trt, Trt:SampleSev:Soil_Subgroup
ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "B"), ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="B reps")

#by horizons 
a.M <- ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A") & pH.df.noUB$Horizon=="M", ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="A reps Mineral")
a.O <- ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "A") & pH.df.noUB$Horizon=="O", ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="A reps Organic")
a.M + a.O


#C reps
summary(aov(pH~Horizon*Trt*SampleSev*Soil_Subgroup*Veg_Type*Unit, data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "C"), ])) 
#significant variables *** SampleSev
ggplot(data=pH.df.noUB[endsWith(pH.df.noUB$Sample_Core_Name, "C"), ], aes(x=Trt, y=pH)) + geom_boxplot() + facet_wrap(~Unit) + labs(title="C reps")
