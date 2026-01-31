#BURN SEVERITY: Identifying taxa that are differentially abundant using corncob

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
sample_data(ps.3)$SampleSev <- factor(sample_data(ps.3)$SampleSev, levels = c(1,2,3,4), ordered = FALSE)

#check phyloseq object
##check metadata
sample_data(ps.3)[1:5, ]
##check taxa table
tax_table(ps.3)[1:5, ]

#create phyloseq object subset by Unit 
temp <- subset_samples(ps.3, Unit =="H")
#remove the sample with severity = severe; there is only one sample in this class
#thus doing statistics on this unequal dataset will lead to many OTUs ID'ed as 
#differentially abundant in the severe burn severity class, thus inflating the 
#OTUs ID'ed as 'fire responders'
ps.3.H <- subset_samples(temp, SampleSev!='4')

ps.3.BL <- subset_samples(ps.3, Unit =="BL")

ps.3.WV <- subset_samples(ps.3, Unit =="WV")

rm(ps.3, temp)


#subset Unit phyloseq objects so that meanabundance cutoff of taxa > 0.0001

##Hole
###ID taxa present meanAbund>0.0001 in H
ps.norm.H = transform_sample_counts(ps.3.H, function(x) x/sum(x)) 
H.prune = taxa_names(filter_taxa(ps.norm.H, function(x) mean(x)>0.0001, TRUE))
length(H.prune) #1408 (compared to 1414 taxa when including sample with severe burn severity, compared to 16070 of full dataset)
H.prune 
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
###Hole Burn Severity####
###NOTE: comparisons to severe may be weird because only 1 sample at Hole was "severe"

####control for horizon 
set.seed(1)
H.DA.fac <- differentialTest(formula = ~ SampleSev + Horizon,
                         phi.formula = ~ SampleSev+Horizon,
                         formula_null = ~ Horizon,
                         phi.formula_null = ~ SampleSev+Horizon,
                         test = "Wald", boot = FALSE,
                         data = ps.H.prune,
                         fdr_cutoff = 0.05)

length(H.DA.fac$significant_taxa) # 225 OTUs ID'ed (compared to 319 OTUs ID'ed 
#when including the severe burn severity.
#when burn severity treated as ordinal variable. 319 OTUs also identified as diff. 
#abundant when treated as a factor.
#When treat burn severity as factor (not ordered) estimates are compared against 
#the first factor (in the case of the Hole unit, unburned samples) 

plot(H.DA.fac) #work on this so the label is tidied - Genus only? Family only? 

H.DA.fac$significant_models[[1]]$coefficients

#### z add relevant taxa to dataframe - when burn severity is ordinal ####
a = data.frame('intercept'=0, 'mu.burnsev'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA.fac
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:2] =H.DA.lib$significant_models[[i]]$coefficients[1:2] #Pull out the coefficients for each OTU
  a[i, 3] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,4] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,5] = 'BS ordinal'
  a[i,6] = 'Hole'
  a[i,7] = 'meanAbund>0.0001'
  a[i,8] = 'horizon'
}

a.concat=a

#ordinal results
SigOTUs = levels(as.factor(a.concat$OTU))
pruned = prune_taxa(SigOTUs,ps.H.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a.concat,taxtab,by=c("OTU"))
head(df.joined)

write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_H_BSord_DA.csv")

#### ####
#add relevant taxa dataframe - when burn severity is factored and compared against 
#lowest burn severity classification (for Hole unit, compared against unburned samples)
# Create empty dataframe to fill
a = data.frame('intercept'=0, 'mu.low'=0, 'mu.moderate'=0, 'mu.Ohorizon'=0, 't.test_p_low'=0, 't.test_p_moderate'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

H.DA.lib <- H.DA.fac
for (i in 1:length(H.DA.lib$significant_taxa)) {
  a[i, 1:4] =H.DA.lib$significant_models[[i]]$coefficients[1:4] #Pull out the coefficients (B1 ) for each OTU
  a[i, 5:6] = H.DA.lib$significant_models[[i]]$coefficients[26:27] #pull out the pvalues for each level comparison to unburned
  a[i, 7] = H.DA.lib$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,8] = H.DA.lib$p_fdr[H.DA.lib$significant_taxa][i]
  a[i,9] = 'BS factor'
  a[i,10] = 'Hole'
  a[i,11] = 'meanAbund>0.0001'
  a[i,12] = 'horizon'
}

# calculate log2 fold change from mu values
a$Rel.Abund.unburned <- (invlogit(a$intercept) +
                          invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.low <- (invlogit(a$intercept + a$mu.low) +
                     invlogit(a$intercept + a$mu.Ohorizon + a$mu.low))/2
a$Rel.Abund.moderate <- (invlogit(a$intercept + a$mu.moderate) +
                          invlogit(a$intercept + a$mu.Ohorizon + a$mu.moderate))/2
a$FC.low <- a$Rel.Abund.low/a$Rel.Abund.unburned
a$FC.moderate <- a$Rel.Abund.moderate/a$Rel.Abund.unburned
a$L2FC.low <- log(a$FC.low, base=2)
a$L2FC.moderate <- log(a$FC.moderate, base=2)

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

df.joined[2:4,]

#save output 
write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_Hole_BSfac_DA.csv")

### compare the results when burn severity is treated as ordinal vs. factored
### NOTE: these results have not addressed the subsampling - (i.e. no random effect
#of plot) 

comp <- merge(a, a.concat, by='OTU') # all but 2 OTUs were identified with both methods.
#indicates that either method is robust in identifying taxa BUT the information 
#gained from treating the levels as non ordinal factors is more precise - there are many 
#OTUs that are not equally diff. abundant with each level - indicating that diff. 
#abundance is not a linear relationship with burn severity.



######## Butte Lake Fire History ########

####control for horizon 
set.seed(1)
DA <- differentialTest(formula = ~ SampleSev + Horizon,
                             phi.formula = ~ SampleSev+Horizon,
                             formula_null = ~ Horizon,
                             phi.formula_null = ~ SampleSev+Horizon,
                             test = "Wald", boot = FALSE,
                             data = ps.BL.prune,
                             fdr_cutoff = 0.05)

length(DA$significant_taxa) # 351 OTUs ID'ed as diff. abundant when treated as a factor.
#When treat burn severity as factor (not ordered) estimates are compared against 
#the first factor (in the case of the Hole unit, unburned samples) at least one comparison 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

#add relevant taxa dataframe 
# Create empty dataframe to fill
a = data.frame('intercept'=0, 'mu.low'=0, 'mu.moderate'=0, 'mu.severe'=0, 'mu.Ohorizon'=0, 't.test_p_low'=0, 't.test_p_moderate'=0, 't.test_p_severe'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:5] =DA$significant_models[[i]]$coefficients[1:5] #Pull out the coefficients for each OTU
  a[i, 6:8] = DA$significant_models[[i]]$coefficients[32:34] #pull out the pvalues for each level comparison to unburned
  a[i, 9] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,10] = DA$p_fdr[DA$significant_taxa][i]
  a[i,11] = 'BS factor'
  a[i,12] = 'Butte Lake'
  a[i,13] = 'meanAbund>0.0001'
  a[i,14] = 'horizon'
}

# calculate log2 fold change from mu values
a$Rel.Abund.unburned <- (invlogit(a$intercept) +
                           invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.low <- (invlogit(a$intercept + a$mu.low) +
                      invlogit(a$intercept + a$mu.Ohorizon + a$mu.low))/2
a$Rel.Abund.moderate <- (invlogit(a$intercept + a$mu.moderate) +
                           invlogit(a$intercept + a$mu.Ohorizon + a$mu.moderate))/2
a$Rel.Abund.severe <- (invlogit(a$intercept + a$mu.severe) +
                         invlogit(a$intercept + a$mu.Ohorizon + a$mu.severe))/2
a$FC.low <- a$Rel.Abund.low/a$Rel.Abund.unburned
a$FC.moderate <- a$Rel.Abund.moderate/a$Rel.Abund.unburned
a$FC.severe <- a$Rel.Abund.severe/a$Rel.Abund.unburned
a$L2FC.low <- log(a$FC.low, base=2)
a$L2FC.moderate <- log(a$FC.moderate, base=2)
a$L2FC.severe <- log(a$FC.severe, base=2)

a[1,]

# Let's bring back in our taxonomy from the tax table
# and save the dataframe as a csv file for later reference

SigOTUs = levels(as.factor(a$OTU))
pruned = prune_taxa(SigOTUs,ps.BL.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a,taxtab,by=c("OTU"))
head(df.joined)

write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob/corncob_BL_BSfac_DA.csv")


######## Warner Valley Fire History ########
####control for horizon 
####NOTE: Here the comparisons are against low burn severity (bc no unburned samples
####existed at Warner Valley)
set.seed(1)
DA <- differentialTest(formula = ~ SampleSev + Horizon,
                       phi.formula = ~ SampleSev+Horizon,
                       formula_null = ~ Horizon,
                       phi.formula_null = ~ SampleSev+Horizon,
                       test = "Wald", boot = FALSE,
                       data = ps.WV.prune,
                       fdr_cutoff = 0.05)

length(DA$significant_taxa) #198  OTUs ID'ed as diff. abundant when treated as a factor.
#When treat burn severity as factor (not ordered) estimates are compared against 
#the first factor (in the case of the Hole unit, unburned samples) at least one comparison 

plot(DA) #work on this so the label is tidied - Genus only? Family only? 

DA$significant_models[[1]]$coefficients

#add relevant taxa dataframe 
# Create empty dataframe to fill
a = data.frame('intercept'=0, 'mu.lowtomod'=0, 'mu.lowtosevere'=0, 'mu.Ohorizon'=0,  't.test_p_moderate'=0, 't.test_p_severe'=0, 'OTU'='taxa','p_fdr'=0,
               'comparison'='compare','unit'='unit', 'cutoff'='text', 'covars'='text')

for (i in 1:length(DA$significant_taxa)) {
  a[i, 1:4] =DA$significant_models[[i]]$coefficients[1:4] #Pull out the coefficients for each OTU
  a[i, 5:6] = DA$significant_models[[i]]$coefficients[26:27] #pull out the pvalues for each level comparison to unburned
  a[i, 7] = DA$significant_taxa[[i]] # pull out OTU name pull out coefficients for each taxon
  # Also grab the p_fdr estimate for that taxon's model
  a[i,8] = DA$p_fdr[DA$significant_taxa][i]
  a[i,9] = 'BS factor'
  a[i,10] = 'Warner Valley'
  a[i,11] = 'meanAbund>0.0001'
  a[i,12] = 'horizon'
}

a[1,]
# calculate log2 fold change from mu values
a$Rel.Abund.low <- (invlogit(a$intercept) +
                      invlogit(a$intercept+a$mu.Ohorizon))/2
a$Rel.Abund.moderate <- (invlogit(a$intercept + a$mu.lowtomod) +
                           invlogit(a$intercept + a$mu.Ohorizon + a$mu.lowtomod))/2
a$Rel.Abund.severe <- (invlogit(a$intercept + a$mu.lowtosevere) +
                         invlogit(a$intercept + a$mu.Ohorizon + a$mu.lowtosevere))/2
a$FC.mod.v.low <- a$Rel.Abund.moderate/a$Rel.Abund.low
a$FC.sev.v.low <- a$Rel.Abund.severe/a$Rel.Abund.low
a$L2FC.mod.v.low <- log(a$FC.mod.v.low, base=2)
a$L2FC.sev.v.low <- log(a$FC.sev.v.low, base=2)

a[1,]

# Let's bring back in our taxonomy from the tax table
# and save the dataframe as a csv file for later reference

SigOTUs = levels(as.factor(a$OTU))
pruned = prune_taxa(SigOTUs,ps.WV.prune)
taxtab = data.frame(tax_table(pruned))
taxtab$OTU = c(taxa_names(pruned))
df.joined = merge(a,taxtab,by=c("OTU"))
head(df.joined)

write.csv(df.joined, "~/Desktop/LAVO22_GG2/corncob_WV_BSfac_DA.csv")

