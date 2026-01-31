# goal: dissimilartity PERMANOVA on Silva assigned taxonomy using fire history, 
# quadrat median burn severity, horizon, and unit as explanatory variables 

library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)

#import phyloseq object
ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")

# 1. remove resequenced samples
#this phyloseq object has duplicates of 6 samples. These duplicates
#(Submission_Sample_Name ending in 2)  were done with the correct BSA 
#concentration to test whether a more dilute BSA concentration from the first 3 
#PCR plates resulted in changes in sequencing depth (correct BSA concentration 
#here actually decreased sequencing depth of the samples). However, there did not
# appear to be a dramatic difference in community composition of the resequenced 
# samples (see . So, we will be removing the duplicates for subsequent analyses
#(LAVO22-BL-PIJE-09AO2, LAVO22-UB-01-AO-2, LAVO22-UB-02-BO-2, LAVO22-UB-02-CO-2,
#LAVO22-WV-ABCO-23CO2, LAVO22-WV-04-CM-2). We will accomplish this by subsetting 
#the samples with those that do not end in "2".

#?subset_samples
ps = subset_samples(ps, !endsWith(Submission_Sample_Name, "2"))
                           
# check that this actually got rid of the samples in question #should have only 
# 6 samples left in the Resequence column instead of 12
with(sample_data(ps), table(Resequence))


# 2. add median quadrat burn severity to metadata (sample_data dataframe)

# read in quadrat burn severity dataframe (created in "MedianBS_Quadrat.R")
bs = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/BurnSeverity/MedianBS_Quadrat.csv")

# create dataframe of Submission_Sample_Name and its corresponding median bs
df_merge = left_join(data.frame(sample_data(ps)), bs, by ="Sample_Core_Name")
#remove extreneous columns
df_merge = df_merge %>%
  select(-Quad_median_burnsev, -Run, -X)
row.names(df_merge) = df_merge$Submission_Sample_Name
sample_data(ps) = sample_data(df_merge)

#check that merge worked 
sample_variables(ps)
# it did! so let's remove the data frames we created
rm(bs)
rm(df_merge)

# 3. remove Unburned controls from object for subsequent modeling
ps = subset_samples(ps,Unit != "UB")

# now that we did all that work cleaning up the phyloseq object, let's save it
#saveRDS(ps,"~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3units")

############ PERMANOVA ############################

# Create ordination

# Get relative abundances
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5 )

# Get distance matrix
dist.hell = phyloseq::distance(ps.hell, method="bray")
#take a peak at the distance matrix
head(dist.hell)

# Get sample_data as dataframe
SamDat = data.frame(sample_data(ps))
# Question: how should I be treating burn severity, as numeric? character? what 
# does this mean for 0.5 values? For now we will leave it as a numeric value 

# test if Trt, horizon, unit, and median_bs are significant
hell.adonis = adonis2(dist.hell ~ Trt + Horizon + Unit + median_bs, SamDat)
# returns error bc missing burn severity for WV 01 C quadrat.

#Question: what should I do in the future for this plot? 

#For now, I will just drop this plot

ps.drop = subset_samples(ps, median_bs>0)
ps.hell2 = transform_sample_counts(ps.drop, function(x) (x / sum(x))^0.5 )
dist.hell2 = phyloseq::distance(ps.hell2, method="bray")
SamDat2 = data.frame(sample_data(ps.drop))

hell.adonis2 = adonis2(dist.hell2 ~ Trt + Horizon + Unit + median_bs, SamDat2)
hell.adonis2
#median_bs is a significant variable in explaining community dissimilarity; 
#however the R2 (or degree to which this variable explains variance) is lowest 
#of all the other explanatory variables (0.0394), even lower than the treatment
#(0.04678) all explanatory variables have a p-value <0.001

#Let's look at all of my relevant explanatory variables
hell.adonis.full = adonis2(dist.hell2 ~ Trt + Horizon + Unit + median_bs + Veg_Type + Overstory_mort, SamDat2)
hell.adonis.full

# 3. When running a PERMANOVA, GustaMe reports:
# "Anderson (2001) warns that groups of objects with different dispersions, 
# yet no significant differences in centres (centres are similar to means, 
# but may be non-Euclidean), may result in misleadingly low P-values. 
# It is thus recommended that the dispersion be evaluated and considered when 
# interpreting the results of NPMANOVA. See Anderson (2006) for a discussion 
# on tests of multivariate dispersion."

# Thus, we should also look at the dispersion of our samples.
b.Trt = betadisper(dist.hell2, SamDat2$Trt)
b.Trt
permutest(b.Trt) # here the results of this permutation test is significant; this means that we cannot reject the null hypothesis that our groups have the same dispersions; we need to be skeptical of the adonis....

b.Horizon = betadisper(dist.hell2, SamDat2$Horizon)
permutest(b.Horizon) #again this test shows significance

b.Unit = betadisper(dist.hell2, SamDat2$Unit)
permutest(b.Unit) #not significant (p = 0.052); we can reject the null that our Units have the same dispersions

b.median_bs = betadisper(dist.hell2, SamDat2$median_bs)
permutest(b.median_bs) #again this is significant

#Let's visualize the dispersion by boxplot of distance to centroid for each factor
par(mfrow = c(2, 2))
boxplot(b.Trt, xlab = "Treatment")
boxplot(b.Horizon, xlab = "Horizon")
boxplot(b.Unit, xlab = "Unit")
boxplot(b.median_bs, xlab = "Quadrat Median Burn Severity")


# Question: what does this test actually show? Read Anderson 2001 paper and look into PERMANOVA assumptions and ways to deal with this


#is there an interaction between burn severity and unit? 
hell.adonis3 = adonis2(dist.hell2 ~ Trt + Horizon + Unit + median_bs + Unit*median_bs, SamDat2)
hell.adonis3
# yes there is an interaction, so the strength of burn severity as an explanatory
# variable varies across the 3 units


################ Subsetting by Unit #######################################
#how does burn severity strength vary across the units? Let's look at the Units 
#separately
H.ps <- subset_samples(ps.drop,Unit == "H")
BL.ps <- subset_samples(ps.drop,Unit == "BL")
WV.ps <- subset_samples(ps.drop,Unit == "WV")

H.hell = transform_sample_counts(H.ps, function(x) (x / sum(x))^0.5 )
BL.hell = transform_sample_counts(BL.ps, function(x) (x / sum(x))^0.5 )
WV.hell = transform_sample_counts(WV.ps, function(x) (x / sum(x))^0.5 )

H.dist = phyloseq::distance(H.hell, method="bray")
BL.dist = phyloseq::distance(BL.hell, method="bray")
WV.dist = phyloseq::distance(WV.hell, method="bray")

H.SamDat = data.frame(sample_data(H.ps))
BL.SamDat = data.frame(sample_data(BL.ps))
WV.SamDat = data.frame(sample_data(WV.ps))

H.adonis = adonis2(H.hell ~ Trt + Horizon + median_bs, H.SamDat)
BL.adonis = adonis2(BL.hell ~ Trt + Horizon + median_bs, BL.SamDat)
WV.adonis = adonis2(WV.hell ~ Trt + Horizon + median_bs, WV.SamDat)
