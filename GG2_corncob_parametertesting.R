#Identifying taxa that are differentially abundant using corncob

#load libraries
library(corncob)
library(phyloseq)
library(patchwork)
library(dplyr)
library(tidyr)

################ LOAD AND PRUNE PHYLOSEQ OBJECTS ###############################

#load phyloseq object
##all 3 units, no positive controls, blanks, or duplicated samples
#ASUS laptop
ps.3 = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")

#MiniMac
#ps.3 = readRDS("~/Desktop/LAVO22_GG2/ps.3units")


#clean up taxa names so OTU does not start with a number 
ps.3 <- clean_taxa_names(ps.3)

#check phyloseq object
##check metadata
sample_data(ps.3)[1:5, ]
##check taxa table
tax_table(ps.3)[1:5, ]

############ ENTIRE DATASET MFvLF ######################
ps.norm = transform_sample_counts(ps.3, function(x) x/sum(x)) 

###ID taxa present meanAbund>0.0001 in entire dataset
ps.prune = taxa_names(filter_taxa(ps.norm, function(x) mean(x)>0.0001, TRUE))
length(ps.prune) #1103 (compared to 16070 of full dataset)

###prune phyloseq object to include only meanAbund>0.0001
ps.prune = prune_taxa(ps.prune, ps.3)
rm(ps.norm)

#Differential Test
set.seed(1)
DA.MFvLF <- differentialTest(formula = ~ Trt,
                               phi.formula = ~ Trt,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Trt,
                               test = "Wald", boot = FALSE,
                               data = ps.prune,
                               fdr_cutoff = 0.05)
DA.MFvLF

length(DA.MFvLF$significant_taxa) #WOO! 364 taxa significant differential abundance

plot(DA.MFvLF, levels=c("Family", "Genus"))

#add to relevant taxa dataframe
# Create empty dataframes to fill
z = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- DA.MFvLF
for (i in 1:length(H.DA.lib$significant_taxa)) {
  z[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  z[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  z[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  z[i,5] = 'MF.vs.LF'
  z[i,6] = 'all'
  z[i,7] = 'meanAbund>0.0001'
  z[i,8] = 'none'
}

#add to relevant taxa Hole dataframe
z 
z.concat = z

####Controlling for Unit and soil horizon 
set.seed(1)
DA.MFvLF.UnitHorizon <- differentialTest(formula = ~ Trt + Unit + Horizon,
                             phi.formula = ~ Trt,
                             formula_null = ~ Unit + Horizon,
                             phi.formula_null = ~ Trt,
                             test = "Wald", boot = FALSE,
                             data = ps.prune,
                             fdr_cutoff = 0.05)
DA.MFvLF.UnitHorizon

length(DA.MFvLF.UnitHorizon$significant_taxa) # 397 taxa significant differential abundance

plot(DA.MFvLF.UnitHorizon, levels=c("Family", "Genus"))


#add to relevant taxa dataframe
# Create empty dataframes to fill; how to fill mu of unit and horizon
z = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- DA.MFvLF.UnitHorizon
for (i in 1:length(H.DA.lib$significant_taxa)) {
  z[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  z[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  z[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  z[i,5] = 'MF.vs.LF'
  z[i,6] = 'all'
  z[i,7] = 'meanAbund>0.0001'
  z[i,8] = 'unit and horizon'
}

#add to relevant taxa Hole dataframe
z 
z.concat=rbind(z.concat, z)

#bring back taxonomy
SigOTUs = levels(as.factor(z.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.3)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(z.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output as csv
#write.csv(df.joined, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2_corncob/fulldataset_DA.csv")
rm(DA.MFvLF, DA.MFvLF.UnitHorizon, z, z.concat, df.joined, SigOTUs, pruned, taxtab)

full_par_DA <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2_corncob/fulldataset_DA.csv")

##summarize output as count of phyla 
full.par.DA.summary <- full_par_DA %>%
  group_by(covariate, Phylum) %>%
  summarise(n=n())

#save output as csv
#write.csv(full.par.DA.summary, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2_corncob/summary_fulldataset_DA.csv")

##how many OTUs show up in both models? 
temp <- full_par_DA %>%
  group_by(OTU) %>%
  filter(n()>1) %>%
  distinct(OTU, .keep_all = TRUE) %>% 
  group_by(Phylum) %>%
  summarise(n=n()) #345 OTUs exist in both models

### summarise phyla by mu (+ or -)
fullDA <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2_corncob/fulldataset_DA.csv")

fullDA$Phylum = as.factor(fullDA$Phylum)

fullDA_cntrlunithorizon <- filter(fullDA, covariate=="unit and horizon")
fullDA_nocntrl <- filter(fullDA, covariate=="none")

sum_fullDA_cntrlunithorizon <- fullDA_cntrlunithorizon %>%
  mutate(mu.pos=ifelse(mu.MF.trtmt<0, 0, 1), mu.neg=ifelse(mu.MF.trtmt>0, 0, 1)) %>%
  group_by(Phylum) %>%
  summarise(mu.pos.n=sum(mu.pos), mu.neg.n=sum(mu.neg))


##############SEPARATE ANALYSES BY UNIT#######################################
#create phyloseq object subset by Unit 
ps.3.H <- subset_samples(ps.3, Unit =="H")

ps.3.BL <- subset_samples(ps.3, Unit =="BL")

ps.3.WV <- subset_samples(ps.3, Unit =="WV")


#subset Unit phyloseq objects so that meanabundance cutoff of taxa > 0.0001
ps.norm = transform_sample_counts(ps.3, function(x) x/sum(x)) 

rm(ps.3)



##Hole
###ID taxa present meanAbund>0.0001 in H
H.prune = taxa_names(filter_taxa(subset_samples(ps.norm, Unit=="H"), function(x) mean(x)>0.0001, TRUE))
length(H.prune) #1414 taxa (compared to 170 of full dataset)

###prune phyloseq object to include only meanAbund>0.0001
ps.H.prune = prune_taxa(H.prune, ps.3.H)
rm(ps.3.H, H.prune)


##Butte Lake
###ID taxa present meanAbund>0.0001 in BL
BL.prune = taxa_names(filter_taxa(subset_samples(ps.norm, Unit=="BL"), function(x) mean(x)>0.0001, TRUE))
length(BL.prune) #929 (compared to 16070 of full dataset)

###prune phyloseq object to include only meanAbund>0.0001
ps.BL.prune = prune_taxa(BL.prune, ps.3.BL)
rm(ps.3.BL, BL.prune)


##Warner Valley
###ID taxa present meanAbund>0.0001 in WV
WV.prune = taxa_names(filter_taxa(subset_samples(ps.norm, Unit=="WV"), function(x) mean(x)>0.0001, TRUE))
length(WV.prune) #913 (compared to 16070 of full dataset)

###prune phyloseq object to include only meanAbund>0.0001
ps.WV.prune = prune_taxa(WV.prune, ps.3.WV)
rm(ps.3.WV, WV.prune)





######################### HOLE ###############################################
#Differential Test
set.seed(1)
H.DA.MFvLF <- differentialTest(formula = ~ Trt,
                              phi.formula = ~ Trt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Trt,
                              test = "Wald", boot = FALSE,
                              data = ps.H.prune,
                              fdr_cutoff = 0.05)
H.DA.MFvLF

length(H.DA.MFvLF$significant_taxa) #WOO! 53 taxa significant differential abundance

plot(H.DA.MFvLF)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- H.DA.MFvLF
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'none'
}

#add to relevant taxa Hole dataframe
a
a.concat= a

####Controlling for horizon###
##What about adding the covariate horizon to control for, horizon
set.seed(1)
H.DA.MFvLF.Hor <- differentialTest(formula = ~ Trt + Horizon,
                               phi.formula = ~ Trt,
                               formula_null = ~ Horizon,
                               phi.formula_null = ~ Trt,
                               test = "Wald", boot = FALSE,
                               data = ps.H.prune,
                               fdr_cutoff = 0.05)
H.DA.MFvLF.Hor

length(H.DA.MFvLF.Hor$significant_taxa) # 53 taxa significant differential abundance

plot(H.DA.MFvLF.Hor)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- H.DA.MFvLF.Hor
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

sample_data(ps.H.prune)
####Controlling for vegetation type###
##What about adding the covariate veg_type to control for veg_type
set.seed(1)
H.DA.MFvLF.Veg <- differentialTest(formula = ~ Trt + Veg_Type,
                                   phi.formula = ~ Trt,
                                   formula_null = ~ Veg_Type,
                                   phi.formula_null = ~ Trt,
                                   test = "Wald", boot = FALSE,
                                   data = ps.H.prune,
                                   fdr_cutoff = 0.05)

length(H.DA.MFvLF.Veg$significant_taxa) # 10 taxa significant differential abundance

plot(H.DA.MFvLF.Veg)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- H.DA.MFvLF.Veg
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'vegtype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Controlling for soil type
set.seed(1)
H.DA.MFvLF.Soil <- differentialTest(formula = ~ Trt + Soil_mu,
                                   phi.formula = ~ Trt,
                                   formula_null = ~ Soil_mu,
                                   phi.formula_null = ~ Trt,
                                   test = "Wald", boot = FALSE,
                                   data = ps.H.prune,
                                   fdr_cutoff = 0.05)

length(H.DA.MFvLF.Soil$significant_taxa) #68 taxa significant differential abundance

plot(H.DA.MFvLF.Soil, level=c("Family", "Genus"))
?plot()

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- H.DA.MFvLF.Soil
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'soiltype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Control for both soil horizon and vegetation type
##vegetation type is not paired at the hole unit within the fire treatments, 
##while this may be a consequence of the history of fire in the unit, this is a 
##conservative approach to try to control for aboveground influences on taxa that 
##are differentially abundant with more or less fire. 
set.seed(1)
H.DA.MFvLF.HorizonVeg <- differentialTest(formula = ~ Trt + Horizon + Veg_Type,
                                    phi.formula = ~ Trt,
                                    formula_null = ~ Horizon + Veg_Type,
                                    phi.formula_null = ~ Trt,
                                    test = "Wald", boot = FALSE,
                                    data = ps.H.prune,
                                    fdr_cutoff = 0.05)

length(H.DA.MFvLF.HorizonVeg$significant_taxa) #13 taxa significant differential abundance

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

H.DA.lib <- H.DA.MFvLF.HorizonVeg
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon and veg'
}

a
a.concat = rbind(a.concat, a)

rm(H.DA.lib, i)

#bring back taxonomy
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.3)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output as csv
#write.csv(df.joined, "~/Desktop/LAVO22_GG2/Hole_parametertesting_DA.csv")
rm(a, a.concat, H.DA.MFvLF, H.DA.MFvLF.Hor, H.DA.MFvLF.HorizonVeg, H.DA.MFvLF.Soil, H.DA.MFvLF.Veg, SigOTUs, pruned, taxtab, df.joined)

######################### BUTTE LAKE ###########################################
#Differential Test
set.seed(1)
BL.DA.MFvLF <- differentialTest(formula = ~ Trt,
                               phi.formula = ~ Trt,
                               formula_null = ~ 1,
                               phi.formula_null = ~ Trt,
                               test = "Wald", boot = FALSE,
                               data = ps.BL.prune,
                               fdr_cutoff = 0.05)

length(BL.DA.MFvLF$significant_taxa) #443 taxa significant differential abundance

plot(BL.DA.MFvLF)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- BL.DA.MFvLF
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'none'
}

#add to relevant taxa Hole dataframe
a
a.concat= a

####Controlling for horizon###
##What about adding the covariate horizon to control for, horizon
set.seed(1)
BL.DA.MFvLF.Hor <- differentialTest(formula = ~ Trt + Horizon,
                                   phi.formula = ~ Trt,
                                   formula_null = ~ Horizon,
                                   phi.formula_null = ~ Trt,
                                   test = "Wald", boot = FALSE,
                                   data = ps.BL.prune,
                                   fdr_cutoff = 0.05)

length(BL.DA.MFvLF.Hor$significant_taxa) #462 taxa significant differential abundance

plot(BL.DA.MFvLF.Hor)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- BL.DA.MFvLF.Hor
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)


####Controlling for vegetation type###
##What about adding the covariate veg_type to control for veg_type
set.seed(1)
BL.DA.MFvLF.Veg <- differentialTest(formula = ~ Trt + Veg_Type,
                                   phi.formula = ~ Trt,
                                   formula_null = ~ Veg_Type,
                                   phi.formula_null = ~ Trt,
                                   test = "Wald", boot = FALSE,
                                   data = ps.BL.prune,
                                   fdr_cutoff = 0.05)

length(BL.DA.MFvLF.Veg$significant_taxa) #433  taxa significant differential abundance

plot(BL.DA.MFvLF.Veg)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- BL.DA.MFvLF.Veg
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'vegtype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Controlling for soil type
set.seed(1)
BL.DA.MFvLF.Soil <- differentialTest(formula = ~ Trt + Soil_mu,
                                    phi.formula = ~ Trt,
                                    formula_null = ~ Soil_mu,
                                    phi.formula_null = ~ Trt,
                                    test = "Wald", boot = FALSE,
                                    data = ps.BL.prune,
                                    fdr_cutoff = 0.05)

length(BL.DA.MFvLF.Soil$significant_taxa) #447 taxa significant differential abundance

plot(BL.DA.MFvLF.Soil, level=c("Family", "Genus"))

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- BL.DA.MFvLF.Soil
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'soiltype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Control for both soil horizon and vegetation type
##vegetation type is not paired at the hole unit within the fire treatments, 
##while this may be a consequence of the history of fire in the unit, this is a 
##conservative approach to try to control for aboveground influences on taxa that 
##are differentially abundant with more or less fire. 
set.seed(1)
BL.DA.MFvLF.HorizonVeg <- differentialTest(formula = ~ Trt + Horizon + Veg_Type,
                                          phi.formula = ~ Trt,
                                          formula_null = ~ Horizon + Veg_Type,
                                          phi.formula_null = ~ Trt,
                                          test = "Wald", boot = FALSE,
                                          data = ps.BL.prune,
                                          fdr_cutoff = 0.05)

length(BL.DA.MFvLF.HorizonVeg$significant_taxa) #447 taxa significant differential abundance

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- BL.DA.MFvLF.HorizonVeg
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon and veg'
}

a
a.concat = rbind(a.concat, a)

rm(DA.lib, i)

#bring back taxonomy
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.3)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output as csv
#write.csv(df.joined, "~/Desktop/LAVO22_GG2/ButteLake_parametertesting_DA.csv")
rm(a, a.concat, BL.DA.MFvLF, BL.DA.MFvLF.Hor, BL.DA.MFvLF.HorizonVeg, BL.DA.MFvLF.Soil, BL.DA.MFvLF.Veg, SigOTUs, pruned, taxtab, df.joined)



####################### WARNER VALLEY ##########################################

#Differential Test
set.seed(1)
WV.DA.MFvLF <- differentialTest(formula = ~ Trt,
                                phi.formula = ~ Trt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Trt,
                                test = "Wald", boot = FALSE,
                                data = ps.WV.prune,
                                fdr_cutoff = 0.05)

length(WV.DA.MFvLF$significant_taxa) #82 taxa significant differential abundance

plot(WV.DA.MFvLF, levels=c("Family", "Genus"))

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- WV.DA.MFvLF
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'none'
}

#add to relevant taxa dataframe
a
a.concat= a

####Controlling for horizon###
##What about adding the covariate horizon to control for, horizon
set.seed(1)
WV.DA.MFvLF.Hor <- differentialTest(formula = ~ Trt + Horizon,
                                    phi.formula = ~ Trt,
                                    formula_null = ~ Horizon,
                                    phi.formula_null = ~ Trt,
                                    test = "Wald", boot = FALSE,
                                    data = ps.WV.prune,
                                    fdr_cutoff = 0.05)

length(WV.DA.MFvLF.Hor$significant_taxa) #109 taxa significant differential abundance

plot(WV.DA.MFvLF.Hor)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- WV.DA.MFvLF.Hor
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)


####Controlling for vegetation type###
##What about adding the covariate veg_type to control for veg_type
set.seed(1)
WV.DA.MFvLF.Veg <- differentialTest(formula = ~ Trt + Veg_Type,
                                    phi.formula = ~ Trt,
                                    formula_null = ~ Veg_Type,
                                    phi.formula_null = ~ Trt,
                                    test = "Wald", boot = FALSE,
                                    data = ps.WV.prune,
                                    fdr_cutoff = 0.05)

length(WV.DA.MFvLF.Veg$significant_taxa) #81  taxa significant differential abundance

plot(WV.DA.MFvLF.Veg)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- WV.DA.MFvLF.Veg
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'vegtype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Controlling for soil type
set.seed(1)
WV.DA.MFvLF.Soil <- differentialTest(formula = ~ Trt + Soil_mu,
                                     phi.formula = ~ Trt,
                                     formula_null = ~ Soil_mu,
                                     phi.formula_null = ~ Trt,
                                     test = "Wald", boot = FALSE,
                                     data = ps.WV.prune,
                                     fdr_cutoff = 0.05)

length(WV.DA.MFvLF.Soil$significant_taxa) #81 taxa significant differential abundance

plot(WV.DA.MFvLF.Soil, level=c("Family", "Genus"))

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- WV.DA.MFvLF.Soil
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'soiltype'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

####Control for both soil horizon and vegetation type
##vegetation type is not paired at the hole unit within the fire treatments, 
##while this may be a consequence of the history of fire in the unit, this is a 
##conservative approach to try to control for aboveground influences on taxa that 
##are differentially abundant with more or less fire. 
set.seed(1)
WV.DA.MFvLF.HorizonVeg <- differentialTest(formula = ~ Trt + Horizon + Veg_Type,
                                           phi.formula = ~ Trt,
                                           formula_null = ~ Horizon + Veg_Type,
                                           phi.formula_null = ~ Trt,
                                           test = "Wald", boot = FALSE,
                                           data = ps.WV.prune,
                                           fdr_cutoff = 0.05)

length(WV.DA.MFvLF.HorizonVeg$significant_taxa) #113 taxa significant differential abundance

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covariate'='covar')

DA.lib <- WV.DA.MFvLF.HorizonVeg
for (i in 1:length(DA.lib$significant_taxa)) {
  a[i, 1:2] =DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA.lib$p_fdr[DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon and veg'
}

a
a.concat = rbind(a.concat, a)

rm(DA.lib, i)

#bring back taxonomy
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.3)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output as csv
#write.csv(df.joined, "~/Desktop/LAVO22_GG2/WarnerValley_parametertesting_DA.csv")
rm(a, a.concat, WV.DA.MFvLF, WV.DA.MFvLF.Hor, WV.DA.MFvLF.HorizonVeg, WV.DA.MFvLF.Soil, WV.DA.MFvLF.Veg, SigOTUs, pruned, taxtab, df.joined)




