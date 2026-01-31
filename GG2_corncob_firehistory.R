#Identifying taxa that are differentially abundant using corncob

#load libraries
library(corncob)
library(phyloseq)
library(patchwork)
library(dplyr)
library(tidyr)

#load phyloseq object
##all 3 units, no positive controls, blanks, or duplicated samples
#ps.3 = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")
ps.3 = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")


#clean up taxa names so OTU does not start with a number 
ps.3 <- clean_taxa_names(ps.3)

#check phyloseq object
##check metadata
sample_data(ps.3)[1:5, ]
##check taxa table
tax_table(ps.3)[1:5, ]

#create phyloseq object subset by Unit 
ps.3.H <- subset_samples(ps.3, Unit =="H")

ps.3.BL <- subset_samples(ps.3, Unit =="BL")

ps.3.WV <- subset_samples(ps.3, Unit =="WV")

rm(ps.3)


#subset Unit phyloseq objects so that meanabundance cutoff of taxa > 0.0001

##Hole
###ID taxa present meanAbund>0.0001 in H
ps.norm.H = transform_sample_counts(ps.3.H, function(x) x/sum(x)) 
H.prune = taxa_names(filter_taxa(ps.norm.H, function(x) mean(x)>0.0001, TRUE))
length(H.prune) #1414 taxa (compared to 16070 of full dataset)
rm(ps.norm.H)

###prune phyloseq object to include only meanAbund>0.0001
ps.H.prune = prune_taxa(H.prune, ps.3.H)
rm(ps.3.H)


##Butte Lake
###ID taxa present meanAbund>0.0001 in BL
ps.norm.BL = transform_sample_counts(ps.3.BL, function(x) x/sum(x))
BL.prune = taxa_names(filter_taxa(ps.norm.BL, function(x) mean(x)>0.0001, TRUE))
length(BL.prune) #929 (compared to 16070 of full dataset)
rm(ps.norm.BL)

###prune phyloseq object to include only meanAbund>0.0001
ps.BL.prune = prune_taxa(BL.prune, ps.3.BL)
rm(ps.3.BL)


##Warner Valley
###ID taxa present meanAbund>0.0001 in WV
ps.norm.WV = transform_sample_counts(ps.3.WV, function(x) x/sum(x))
WV.prune = taxa_names(filter_taxa(ps.norm.WV, function(x) mean(x)>0.0001, TRUE))
length(WV.prune) #913 (compared to 16070 of full dataset)
rm(ps.norm.WV)

###prune phyloseq object to include only meanAbund>0.0001
ps.WV.prune = prune_taxa(WV.prune, ps.3.WV)
rm(ps.3.WV)


######### Differential Tests (cutoff mean abund > 0.0001) ###########
###Hole Fire history ####
####control for horizon 
set.seed(1)
H.DA <- differentialTest(formula = ~ Trt + Horizon,
                              phi.formula = ~ Trt+Horizon,
                              formula_null = ~ Horizon,
                              phi.formula_null = ~ Trt+Horizon,
                              test = "Wald", boot = FALSE,
                              data = ps.H.prune,
                              fdr_cutoff = 0.05)
H.DA

length(H.DA$significant_taxa) # 47 OTUs differentially abundant with more fire
#(compare to 53 identified without controlling for Horizon) 

plot(H.DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'mu.Ohorizon'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:3] =H.DA.lib$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  a[i, 4] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,5] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,6] = 'MF.vs.LF'
  a[i,7] = 'Hole'
  a[i,8] = 'meanAbund>0.0001'
  a[i,9] = 'horizon'
}

H.DA$significant_models[[1]]$coefficients
a[1,]

# calculate log2 fold change from mu values
a$Rel.Abund.LF <- (invlogit(a$intercept) +
                     invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.MF <- (invlogit(a$intercept + a$mu.MF.trtmt) +
                      invlogit(a$intercept + a$mu.Ohorizon + a$mu.MF.trtmt))/2
a$FC <- a$Rel.Abund.MF/a$Rel.Abund.LF
a$L2FC <- log(a$FC, base=2)

a[1,]

#Hole relevant taxa
# Let's bring back in our taxonomy from the tax table
# and save the dataframe as a csv file for later reference

#factored results
SigOTUs = levels(as.factor(a$OTU))
pruned = prune_taxa(SigOTUs,ps.H.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob/corncob_Hole_DA.csv")

####fire history alone (no covars)
set.seed(1)
H.DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.H.prune,
                         fdr_cutoff = 0.05)
H.DA

length(H.DA$significant_taxa) # 53 OTUs differentially abundant with more fire
#(compare to 47 identified when controlling for Horizon) 

plot(H.DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA
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
a.concat= rbind(a.concat, a)

a.concat <- read.csv("~/Desktop/LAVO22_GG2/corncob_HoleDA.csv")

#should I subset by Horizon? 
ps.H.O <- subset_samples(ps.H.prune, Horizon=="O")
ps.H.M <- subset_samples(ps.H.prune, Horizon=="M")
#### O horizon only fire history alone (no covars)
set.seed(1)
H.DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.H.O,
                         fdr_cutoff = 0.05)
H.DA

length(H.DA$significant_taxa) # 28 (compare to 47 identified when controlling for Horizon) 

plot(H.DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'O horizon only'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

#### M horizon fire history alone (no covars)
set.seed(1)
H.DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.H.M,
                         fdr_cutoff = 0.05)
H.DA

length(H.DA$significant_taxa) # 23 (compare to 47 identified when controlling for Horizon, compare to 28 from O only) 

plot(H.DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'M horizon only'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)


#Hole relevant taxa
# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.H.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_HoleDA.csv")

### compare the results when including horizon as a covariate.
### not sure how to account for the nestedness here; i.e. haven't controlled for 
### replicates within a plot.

H.compare.covars <- df.joined %>%
  mutate(dup=duplicated(OTU))

O.M.only <- filter(df.joined, covars=='O horizon only' | covars=='M horizon only') %>%
  mutate(dup=duplicated(OTU))

hor <- filter(df.joined, covars!='none')

dup <- H.compare.covars[unlist(tapply(1:nrow(H.compare.covars), H.compare.covars$OTU, function(x) if(length(x)>1) x)), ]
#31 OTUs exist in both the differential tests (when controlling for horizon and not controlling for it....)

horizondup <- O.M.only[unlist(tapply(1:nrow(O.M.only), O.M.only$OTU, function(x) if(length(x)>1) x)), ]
#no duplicates; there are no OTUs that are differentially abundant in both the O & M horizon with treatment.

hordup <- hor[unlist(tapply(1:nrow(hor), hor$OTU, function(x) if(length(x)>1) x)), ]
#21 of the 51 OTUs identified as differentially abundant when testing for horizons separately were 
#also differentially abundant when including horizon as a covariate. (8 of O horizon and 13 of M horizon)





######## Butte Lake Fire History ########

####control for horizon 
set.seed(1)
DA <- differentialTest(formula = ~ Trt + Horizon,
                         phi.formula = ~ Trt+Horizon,
                         formula_null = ~ Horizon,
                         phi.formula_null = ~ Trt+Horizon,
                         test = "Wald", boot = FALSE,
                         data = ps.BL.prune,
                         fdr_cutoff = 0.05)

length(DA$significant_taxa) # 432 OTUs differentially abundant with more fire when controlling for horizon

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'mu.Ohorizon'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:3] =DA$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  a[i, 4] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,5] = DA$p_fdr[DA$significant_taxa][i]
  a[i,6] = 'MF.vs.LF'
  a[i,7] = 'Butte Lake'
  a[i,8] = 'meanAbund>0.0001'
  a[i,9] = 'horizon'
}

DA$significant_models[[1]]$coefficients
a[1,]

# calculate log2 fold change from mu values
a$Rel.Abund.LF <- (invlogit(a$intercept) +
                     invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.MF <- (invlogit(a$intercept + a$mu.MF.trtmt) +
                     invlogit(a$intercept + a$mu.Ohorizon + a$mu.MF.trtmt))/2
a$FC <- a$Rel.Abund.MF/a$Rel.Abund.LF
a$L2FC <- log(a$FC, base=2)

a[1,]

# Let's bring back in our taxonomy from the tax table
# and save the dataframe as a csv file for later reference

#factored results
SigOTUs = levels(as.factor(a$OTU))
pruned = prune_taxa(SigOTUs,ps.BL.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob/corncob_BLDA.csv")

#add to relevant taxa BL dataframe
a
a.concat=a
#a.concat= rbind(a.concat, a)

####fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.BL.prune,
                         fdr_cutoff = 0.05)

length(DA$significant_taxa) # 443 OTUs differentially abundant with more fire
#(compare to 432 identified when controlling for Horizon) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'none'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)


#should I subset by Horizon? 
ps.BL.O <- subset_samples(ps.BL.prune, Horizon=="O")
ps.BL.M <- subset_samples(ps.BL.prune, Horizon=="M")
#### O horizon only fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.BL.O,
                         fdr_cutoff = 0.05)

length(DA$significant_taxa) # 292 differentially abundant in O horizon (compare to 432 identified when controlling for Horizon) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')


for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'O horizon only'
}

#add to relevant taxa BL dataframe
a
a.concat= rbind(a.concat, a)

#### M horizon fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                         phi.formula = ~ Trt,
                         formula_null = ~ 1,
                         phi.formula_null = ~ Trt,
                         test = "Wald", boot = FALSE,
                         data = ps.BL.M,
                         fdr_cutoff = 0.05)

length(DA$significant_taxa) # 341 (compare to 432 identified when controlling for Horizon, compare to 292 from O only) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Butte Lake'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'M horizon only'
}

#add to relevant taxa BL dataframe
a
a.concat= rbind(a.concat, a)


#BL relevant taxa
# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.BL.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_BLDA.csv")


######## Warner Valley Fire History ########
####control for horizon 
set.seed(1)
DA <- differentialTest(formula = ~ Trt + Horizon,
                       phi.formula = ~ Trt+Horizon,
                       formula_null = ~ Horizon,
                       phi.formula_null = ~ Trt+Horizon,
                       test = "Wald", boot = FALSE,
                       data = ps.WV.prune,
                       fdr_cutoff = 0.05)

length(DA$significant_taxa) # 115 OTUs differentially abundant with more fire when controlling for horizon

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'mu.Ohorizon'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:3] =DA$significant_models[[i]]$coefficients[1:3] #Pull out the coefficients for each OTU
  a[i, 4] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,5] = DA$p_fdr[DA$significant_taxa][i]
  a[i,6] = 'MF.vs.LF'
  a[i,7] = 'Warner Valley'
  a[i,8] = 'meanAbund>0.0001'
  a[i,9] = 'horizon'
}

DA$significant_models[[1]]$coefficients
a[1,]

# calculate log2 fold change from mu values
a$Rel.Abund.LF <- (invlogit(a$intercept) +
                     invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.MF <- (invlogit(a$intercept + a$mu.MF.trtmt) +
                     invlogit(a$intercept + a$mu.Ohorizon + a$mu.MF.trtmt))/2
a$FC <- a$Rel.Abund.MF/a$Rel.Abund.LF
a$L2FC <- log(a$FC, base=2)

a[1,]

# Let's bring back in our taxonomy from the tax table
# and save the dataframe as a csv file for later reference

#factored results
SigOTUs = levels(as.factor(a$OTU))
pruned = prune_taxa(SigOTUs,ps.WV.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob/corncob_WVDA.csv")

#add to relevant taxa WV dataframe
a
a.concat=a
#a.concat= rbind(a.concat, a)

####fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                       phi.formula = ~ Trt,
                       formula_null = ~ 1,
                       phi.formula_null = ~ Trt,
                       test = "Wald", boot = FALSE,
                       data = ps.WV.prune,
                       fdr_cutoff = 0.05)

length(DA$significant_taxa) # 82 OTUs differentially abundant with more fire
#(compare to 115 identified when controlling for Horizon) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'none'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)


#should I subset by Horizon? 
ps.WV.O <- subset_samples(ps.WV.prune, Horizon=="O")
ps.WV.M <- subset_samples(ps.WV.prune, Horizon=="M")
#### O horizon only fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                       phi.formula = ~ Trt,
                       formula_null = ~ 1,
                       phi.formula_null = ~ Trt,
                       test = "Wald", boot = FALSE,
                       data = ps.WV.O,
                       fdr_cutoff = 0.05)

length(DA$significant_taxa) # 59 differentially abundant in O horizon (compare to 115 identified when controlling for Horizon) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')


for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'O horizon only'
}

#add to relevant taxa BL dataframe
a
a.concat= rbind(a.concat, a)

#### M horizon fire history alone (no covars)
set.seed(1)
DA <- differentialTest(formula = ~ Trt,
                       phi.formula = ~ Trt,
                       formula_null = ~ 1,
                       phi.formula_null = ~ Trt,
                       test = "Wald", boot = FALSE,
                       data = ps.WV.M,
                       fdr_cutoff = 0.05)

length(DA$significant_taxa) # 67 (compare to 115 identified when controlling for Horizon, compare to 59 from O only) 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:2] =DA$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = DA$p_fdr[DA$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Warner Valley'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'M horizon only'
}

#add to relevant taxa WV dataframe
a
a.concat= rbind(a.concat, a)


#WV relevant taxa
# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.WV.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_WVDA.csv")


############ Z Sensitivity testing to mean abundance cutoff #############
#absolute bare minimum test: test for differential abundance across Trt, without
#controlling for anything else

#looking just at Hole unit for speed
ps.3.H <- ps.3 %>%
  phyloseq::subset_samples(Unit =="H")



sample_data(ps.3.H)[1:5, ]
tax_table(ps.3.H)[1:5, ]

#let's just try on 1 OTU
#null hypothesis: OTU1 is not associated with fire history treatment
H.OTU1.null = bbdml(formula)

H.OTU1 = bbdml(formula = OTU1 ~ Trt, 
               phi.formula = ~1,
               data=ps.3.H)

summary(H.OTU1)
#order of factors is LF, MF

plot(H.OTU1, color="Trt")
#that works, not sure why when doing multiple OTUs at once it fails, see below

#set seed to make the random generator in differentialTest reproducible
set.seed(1)
da_analysis <- differentialTest(formula = ~ Trt,
                                phi.formula = ~ 1,
                                formula_null = ~ 1,
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps.3.H,
                                fdr_cutoff = 0.05)
da_analysis

da_analysis$significant_taxa

plot(da_analysis)

#received "Warning: Separation detected in abundance model! Likely one of your 
#covariates/experimental conditions is such that there are all zero counts within
#a group. The results of this model should be interpreted with care because there
#is insufficient data to distinguish between groups. 

#^could this be failing because I have not set a cutoff of relative abundance? 
#What is the mean relative abundance of each OTU within the Hole unit?


############## CUTOFF 1 sum(count) > 0 #########################################

otu_table(ps.norm.H)[1:3, ]

ps.H.notthere = filter_taxa(ps.3.H, function(x) sum(x)<1, TRUE)

#visualize otu table as dataframe to double check
ps.H.notthere.m = as(otu_table(ps.H.notthere), "matrix")
ps.H.notthere.df = as.data.frame(ps.H.notthere.m)

ps.H.notthere.df
#excellent! there are 7155 taxa that do not exist in the Hole unit. Let's remove
#them from the ps object being used to determine differential abundance. This 
#is the most liberal test of differential abundance. Often, people will set a 
#cutoff of mean abundance (Thea PyOM Neon cutoff > 0.0001, Dana MS cutoff > 0.000001, 
#Dana WoodBuffalo paper > 0.0001)
rm(ps.H.notthere.df, ps.H.notthere.m)

#remove OTUs that don't exist in the Hole unit
ps.H.DA = prune_taxa(taxa_sums(ps.3.H)>0, ps.3.H)

#dummy checks - how many taxa? how many samples? 
ntaxa(ps.H.DA) #8915 taxa
ntaxa(ps.3.H) #16070 taxa
ntaxa(ps.H.notthere) #7155 taxa
16070-7155 #good math! good dummy! it worked!
nsamples(ps.H.DA) #63
nsamples(ps.3.H) #63


#Now let's try a differential abundance test, simple only Trt as covariate

set.seed(1)
da_analysis <- differentialTest(formula = ~ Trt,
                                phi.formula = ~ Trt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Trt,
                                test = "Wald", boot = FALSE,
                                data = ps.H.DA,
                                fdr_cutoff = 0.05)
da_analysis 

da_analysis$significant_taxa #WOO! 31 taxa significant differential abundance
da_analysis$p_fdr

plot(da_analysis)

#what is the order of factors for Trt? See line 42 - order is LF, MF
#This means that for taxa with + values, they are sig more abundant in MF than LF
# - values -> taxa are sig more abundant in LF than MF

#From NEON_PyOM DA_Corncob_16S.iypnb

#save output of relevant taxa (those ID'ed to be significant in differential abundance)
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text')

H.DA.lib <- da_analysis
for (i in 1:length(da_analysis$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'abundance>0'
}

a.concat=a

################## CUTOFF 3 mean abundance > 0.0001 ##########################
H.cutoff.0001 = taxa_names(filter_taxa(ps.norm.H, function(x) mean(x)>0.0001, TRUE))
length(H.cutoff.0001) #1414 taxa (compared to 8915 taxa in most liberal cutoff)

ps.H.cutoff.0001 = prune_taxa(H.cutoff.0001, ps.3.H)

set.seed(1)
H.DA.0001 <- differentialTest(formula = ~ Trt,
                                phi.formula = ~ Trt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Trt,
                                test = "Wald", boot = FALSE,
                                data = ps.H.cutoff.0001,
                                fdr_cutoff = 0.05)
H.DA.0001 

length(H.DA.0001$significant_taxa) #WOO! 53 taxa significant differential abundance

plot(H.DA.0001)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text')

H.DA.lib <- H.DA.0001
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)

################## CUTOFF 2 mean abundance > 0.000001 ########################
H.cutoff.000001 = taxa_names(filter_taxa(ps.norm.H, function(x) mean(x)>0.000001, TRUE))
length(H.cutoff.000001) #8120 taxa (compared to 8915 taxa in most liberal cutoff)

ps.H.cutoff.000001 = prune_taxa(H.cutoff.000001, ps.3.H)
rm(H.cutoff.000001)

set.seed(1)
H.DA.000001 <- differentialTest(formula = ~ Trt,
                              phi.formula = ~ Trt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Trt,
                              test = "Wald", boot = FALSE,
                              data = ps.H.cutoff.000001,
                              fdr_cutoff = 0.05)
H.DA.000001 

length(H.DA.000001$significant_taxa) #WOO! 31  taxa significant differential abundance

plot(H.DA.000001)


#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text')

H.DA.lib <- H.DA.000001
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.000001'
}

#add to relevant taxa Hole dataframe
a.concat= rbind(a.concat, a)

#Hole relevant taxa
# Let's bring back in our taxonomy from the tax table
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.3.H)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)






### Step 7. Save results -----

write.csv(df.joined,"~/Box/MSWhitman/LAVO_FireHistory/Data/R/Hole_DA_taxa.csv", row.names = FALSE)


#let's look at which OTUs were differentially abundant in all of the cutoff tests
H.da.taxa = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Hole_DA_taxa.csv", header=TRUE)

H.da.taxa.allcutoff = H.da.taxa %>%
  mutate(instance = 1) %>%
  group_by(OTU) %>%
  summarise(sensitivity_test_occurrences=sum(instance))

H.da.taxa = left_join(H.da.taxa, H.da.taxa.allcutoff, by="OTU")

rm(H.da.taxa.allcutoff)

#save as new csv
write.csv(H.da.taxa,"~/Box/MSWhitman/LAVO_FireHistory/Data/R/Hole_DA_taxa.csv", row.names = FALSE)

  

############## CONTROL FOR HORIZON ############################################
H.cutoff.0001 = taxa_names(filter_taxa(ps.norm.H, function(x) mean(x)>0.0001, TRUE))
length(H.cutoff.0001) #1414 taxa (compared to 8915 taxa in most liberal cutoff)

ps.H.cutoff.0001 = prune_taxa(H.cutoff.0001, ps.3.H)

set.seed(1)
H.DA.0001 <- differentialTest(formula = ~ Trt,
                              phi.formula = ~ Trt,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Trt,
                              test = "Wald", boot = FALSE,
                              data = ps.H.cutoff.0001,
                              fdr_cutoff = 0.05)
H.DA.0001 

length(H.DA.0001$significant_taxa) #WOO! 53 taxa significant differential abundance

plot(H.DA.0001)

#add to relevant taxa dataframe
# Create empty dataframes to fill
a = data.frame('intercept'=0, 'mu.MF.trtmt'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text')

H.DA.lib <- H.DA.0001
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'MF.vs.LF'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
}

#add to relevant taxa Hole dataframe
a
a.concat= rbind(a.concat, a)










