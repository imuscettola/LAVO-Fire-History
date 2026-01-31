#Corncob Tutorial 
#https://statdivlab.github.io/corncob/articles/corncob-intro.html

#Corncob is a B-logistic regression model that uses abundance tables and sample
#data to determine taxa that are differentially abundant and/or differentially 
#variable across your samples. 

#Install corncob
#install.packages("remotes")
#remotes::install_github("statdivlab/corncob")


#load example data as three different data frames and then combine into a phyloseq
#object
library(corncob)
library(phyloseq)
library(magrittr)
library(patchwork)
data(soil_phylo_sample)
data(soil_phylo_otu)
data(soil_phylo_taxa)
soil_phylo <- phyloseq::phyloseq(phyloseq::sample_data(soil_phylo_sample),
                                 phyloseq::otu_table(soil_phylo_otu, taxa_are_rows = TRUE),
                                 phyloseq::tax_table(soil_phylo_taxa))

#check phyloseq object, should have otu_table, sample_data, and tax_table
soil_phylo

#what are the sample variables? - Plants, DayAmdmt, Amdmt, ID, and Day
sample_data(soil_phylo)[1:3, ]

#check taxa table, should have 7 taxonomic ranks Also make sure that the ranks
#don't start with a number, that could be problematic for corncob
tax_table(soil_phylo)[1:3, ]

################## FITTING A MODEL #########################################

#subset our samples to only those with DayAmdmt equal to 11 or 21 and then 
#collapse samples to the phylum level - this isn't necessary because you can 
#look for differential abundance at any taxonomic level but is useful to expedite
#the tutorial
soil <- soil_phylo %>%
  phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
  phyloseq::tax_glom("Phylum")

#check
soil
tax_table(soil)[1:5, ]


##### Without covariates
#fit the model demonstrating with Proteobacteria or OTU.1 and without covariates 
corncob <- bbdml(formula = OTU.1 ~1, 
                 phi.formula = ~1, 
                 data=soil)

#visualize the model
plot(corncob, B=50)

# x-axis = samples, y-axis = relative abundances, point = relative abundance of 
# the OTU for the sample
# bars = 95% prediction intervals for the observed relative abundance by sample
# B = # of bootstrap simulations used to approximate the prediction intervals
### recommend higher setting for more accurate intervals, default B=1000

#visualize model with y-axis being counts (quasi-absolute abundance) rather than
#relative abundance
plot(corncob, total=TRUE, B=50)

#visualize model coloring by DayAmdmt covariate
a <- plot(corncob,  color = "DayAmdmt", B=50)

###########Adding Covariates
corncob_da <- bbdml(formula = OTU.1 ~ DayAmdmt, 
                    phi.formula = ~ DayAmdmt,
                    data=soil)

b <- plot(corncob_da, color = "DayAmdmt", B=50)

a+b

#b appears to provide a better fit to the data

###################### MODEL SELECTION #####################################

#likelihood ratio test to select final model for the taxa
#null hypothesis that the likelihood of the model with covariates is equal to the 
#likelihood of the model without covariates

lrtest(mod_null = corncob, mod = corncob_da)

#p-value is 4.5e-05 -> conclude that there is a statistically significant 
#difference in the likelihood of the two models -> use the model with covariates
#for this taxon.

########################## PARAMETER INTERPRETATION ###########################
summary(corncob_da)

#dayAmdmt21 abundance coefficient is negative and statistically significant ->
#suggests that the taxon is differentially abundant across DayAmdmt, and that 
#samples with DayAmdmt = 21 are expected to have lower relative abundance.
#also has - and sig coefficient for dispersion -> taxon is differentially 
#variable across DayAmdmt and DayAmdmt =21 expected to have lower variability


##############################################################################
##############################################################################

################ MULTIPLE TAXA ANALYSES ######################################

# use differentialTest function
# perform the tests on all taxa and it will control the false discovery rate to 
# account for multiple comparisons

# specify the model using formula and phi.formula, but don't include the response
# term bc testing multiple taxa. Specify which covariates we want to test for 
# by removing them in the formual_null and phi.formula_null arguments


# test of differential abundance across the DayAmdmt coefficient
# fdr_cutoff is controlled false discovery rate

set.seed(1)
da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ DayAmdmt,
                                test = "Wald", boot = FALSE,
                                data = soil,
                                fdr_cutoff = 0.05)

da_analysis

#list of differentially-abundant taxa: 
da_analysis$significant_taxa









