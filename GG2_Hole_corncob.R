#Identifying taxa that are differentially abundant using corncob

#load libraries
library(corncob)
library(phyloseq)
library(patchwork)
library(dplyr)
library(tidyr)

####################### FIRST PASS #############################################

#load phyloseq object
#all 3 units, no positive controls, blanks, or duplicated samples
ps.3 = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")


#check metadata
sample_data(ps.3)[1:5, ]

#check taxa table
tax_table(ps.3)[1:5, ]
#some otus start with a number, will be an issue for using corncob -> clean_taxa_names(ps)
ps.3 <- clean_taxa_names(ps.3)
tax_table(ps.3)[1:5, ]

#absolute bare minimum test: test for differential abundance across Trt, without
#controlling for anything else

#looking just at Hole unit for speed
ps.3.H <- ps.3 %>%
  phyloseq::subset_samples(Unit =="H")

ps.norm.H = transform_sample_counts(ps.3.H, function(x) x/sum(x)) 

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










