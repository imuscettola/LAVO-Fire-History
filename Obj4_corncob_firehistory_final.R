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

#### Subset Phyloseq by Unit & select OTUs with mean abundance > 0.0001 ####
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
