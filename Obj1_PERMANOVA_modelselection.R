#Objective: To determine which covariates explain significant amounts of variance
#in bacterial community composition 

#################Loading libraries and phyloseq object########################
install.packages("dplyr")

library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(vegan)

#load phyloseq object
#ps.3 <- readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")
#minimac
ps.3 <- readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")

#cohort_df <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/cohort_df.csv")
#minimac
cohort_df <- read.csv("~/Desktop/LAVO22_GG2/cohort_df.csv")
cohort_df_long <- gather(cohort_df, type, Sample, O_Samples:M_Samples, factor_key=TRUE) %>%
  filter(., Sample!="NA")

#hellinger transformation for PERMANOVAs
hell.ps.3 = transform_sample_counts(ps.3, function(x) (x / sum(x))^0.5 )

#make Soil type a factored variable instead of integer
sample_data(hell.ps.3)$Soil_mu = as.factor(sample_data(hell.ps.3)$Soil_mu)

sample_data(hell.ps.3)

rm(ps.3, cohort_df)

#check sample data 
head(sample_data(hell.ps.3))

#removing the one sample with missing pH 
ps.pH <- subset_samples(hell.ps.3, pH!="NA")
cohort_df_long_pH <- cohort_df_long %>%
  filter(Sample!="BL_ABCO_50_BO")

#########Avg OTU table Separate Unit and Horizon: All Covariates alone #########
ps.BL.O <- subset_samples(hell.ps.3, Unit=="BL" & Horizon=="O")
ps.H.O <- subset_samples(hell.ps.3, Unit=="H" & Horizon=="O")
ps.WV.O <- subset_samples(hell.ps.3, Unit=="WV"& Horizon=="O")
ps.BL.M <- subset_samples(hell.ps.3, Unit=="BL" & Horizon=="M")
ps.H.M <- subset_samples(hell.ps.3, Unit=="H" & Horizon=="M")
ps.WV.M <- subset_samples(hell.ps.3, Unit=="WV"& Horizon=="M")

# Get sample_data as dataframe for each phyloseq object
SamDat.BL.O = data.frame(sample_data(ps.BL.O))
SamDat.H.O = data.frame(sample_data(ps.H.O))
SamDat.WV.O = data.frame(sample_data(ps.WV.O))
SamDat.BL.M = data.frame(sample_data(ps.BL.M))
SamDat.H.M = data.frame(sample_data(ps.H.M))
SamDat.WV.M = data.frame(sample_data(ps.WV.M))

#make dissimilarity matrix
dist.BL.O = phyloseq::distance(ps.BL.O, method="bray")
dist.H.O = phyloseq::distance(ps.H.O, method="bray")
dist.WV.O = phyloseq::distance(ps.WV.O, method="bray")
dist.BL.M = phyloseq::distance(ps.BL.M, method="bray")
dist.H.M = phyloseq::distance(ps.H.M, method="bray")
dist.WV.M = phyloseq::distance(ps.WV.M, method="bray")

###only Trt
#Butte Lake
BL.M = adonis2(dist.BL.M ~ Trt, SamDat.BL.M)
BL.M$term=row.names(BL.M)
BL.M = BL.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.O = adonis2(dist.BL.O ~ Trt, SamDat.BL.O)
BL.O$term=row.names(BL.O)
BL.O = BL.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.M = adonis2(dist.H.M ~ Trt, SamDat.H.M)
H.M$term=row.names(H.M)
H.M = H.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.O = adonis2(dist.H.O ~ Trt, SamDat.H.O)
H.O$term=row.names(H.O)
H.O = H.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.M = adonis2(dist.WV.M ~ Trt, SamDat.WV.M)
WV.M$term = row.names(WV.M)
WV.M = WV.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.O = adonis2(dist.WV.O ~ Trt, SamDat.WV.O)
WV.O$term = row.names(WV.O)
WV.O = WV.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###SampleSev
#Butte Lake
BL.bs.M = adonis2(dist.BL.M ~ SampleSev, SamDat.BL.M)
BL.bs.M$term=row.names(BL.bs.M)
BL.bs.M = BL.bs.M %>%
  mutate(cohort=i, covars="SampleSev", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.bs.O = adonis2(dist.BL.O ~ SampleSev, SamDat.BL.O)
BL.bs.O$term=row.names(BL.bs.O)
BL.bs.O = BL.bs.O %>%
  mutate(cohort=i, covars="SampleSev", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.bs.M = adonis2(dist.H.M ~ SampleSev, SamDat.H.M)
H.bs.M$term=row.names(H.bs.M)
H.bs.M = H.bs.M %>%
  mutate(cohort=i, covars="SampleSev", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.bs.O = adonis2(dist.H.O ~ SampleSev, SamDat.H.O)
H.bs.O$term=row.names(H.bs.O)
H.bs.O = H.bs.O %>%
  mutate(cohort=i, covars="SampleSev", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.bs.M = adonis2(dist.WV.M ~ SampleSev, SamDat.WV.M)
WV.bs.M$term = row.names(WV.bs.M)
WV.bs.M = WV.bs.M %>%
  mutate(cohort=i, covars="SampleSev", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.bs.O = adonis2(dist.WV.O ~ SampleSev, SamDat.WV.O)
WV.bs.O$term = row.names(WV.bs.O)
WV.bs.O = WV.bs.O %>%
  mutate(cohort=i, covars="SampleSev", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Avg Quad Sev
#Butte Lake
BL.quadbs.M = adonis2(dist.BL.M ~ Wt_Avg_QuadBS, SamDat.BL.M)
BL.quadbs.M$term=row.names(BL.quadbs.M)
BL.quadbs.M = BL.quadbs.M %>%
  mutate(cohort=i, covars="QuadBS", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.quadbs.O = adonis2(dist.BL.O ~ Wt_Avg_QuadBS, SamDat.BL.O)
BL.quadbs.O$term=row.names(BL.quadbs.O)
BL.quadbs.O = BL.quadbs.O %>%
  mutate(cohort=i, covars="QuadBS", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.quadbs.M = adonis2(dist.H.M ~ Wt_Avg_QuadBS, SamDat.H.M)
H.quadbs.M$term=row.names(H.quadbs.M)
H.quadbs.M = H.quadbs.M %>%
  mutate(cohort=i, covars="QuadBS", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.quadbs.O = adonis2(dist.H.O ~ Wt_Avg_QuadBS, SamDat.H.O)
H.quadbs.O$term=row.names(H.quadbs.O)
H.quadbs.O = H.quadbs.O %>%
  mutate(cohort=i, covars="QuadBS", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.quadbs.M = adonis2(dist.WV.M ~ Wt_Avg_QuadBS, SamDat.WV.M)
WV.quadbs.M$term = row.names(WV.quadbs.M)
WV.quadbs.M = WV.quadbs.M %>%
  mutate(cohort=i, covars="QuadBS", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.quadbs.O = adonis2(dist.WV.O ~ Wt_Avg_QuadBS, SamDat.WV.O)
WV.quadbs.O$term = row.names(WV.quadbs.O)
WV.quadbs.O = WV.quadbs.O %>%
  mutate(cohort=i, covars="QuadBS", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Soil_mu
#Butte Lake
BL.soil.mu.M = adonis2(dist.BL.M ~ Soil_mu, SamDat.BL.M)
BL.soil.mu.M$term=row.names(BL.soil.mu.M)
BL.soil.mu.M = BL.soil.mu.M %>%
  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.soil.mu.O = adonis2(dist.BL.O ~ Soil_mu, SamDat.BL.O)
BL.soil.mu.O$term=row.names(BL.soil.mu.O)
BL.soil.mu.O = BL.soil.mu.O %>%
  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.soil.mu.M = adonis2(dist.H.M ~ Soil_mu, SamDat.H.M)
H.soil.mu.M$term=row.names(H.soil.mu.M)
H.soil.mu.M = H.soil.mu.M %>%
  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.soil.mu.O = adonis2(dist.H.O ~ Soil_mu, SamDat.H.O)
H.soil.mu.O$term=row.names(H.soil.mu.O)
H.soil.mu.O = H.soil.mu.O %>%
  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.soil.mu.M = adonis2(dist.WV.M ~ Soil_mu, SamDat.WV.M)
WV.soil.mu.M$term = row.names(WV.soil.mu.M)
WV.soil.mu.M = WV.soil.mu.M %>%
  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.soil.mu.O = adonis2(dist.WV.O ~ Soil_mu, SamDat.WV.O)
WV.soil.mu.O$term = row.names(WV.soil.mu.O)
WV.soil.mu.O = WV.soil.mu.O %>%
  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Soil Order
#Butte Lake
BL.soil.order.M = adonis2(dist.BL.M ~ Soil_Order, SamDat.BL.M)
BL.soil.order.M$term=row.names(BL.soil.order.M)
BL.soil.order.M = BL.soil.order.M %>%
  mutate(cohort=i, covars="Soil_Order", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.soil.order.O = adonis2(dist.BL.O ~ Soil_Order, SamDat.BL.O)
BL.soil.order.O$term=row.names(BL.soil.order.O)
BL.soil.order.O = BL.soil.order.O %>%
  mutate(cohort=i, covars="Soil_Order", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.soil.order.M = adonis2(dist.H.M ~ Soil_Order, SamDat.H.M)
H.soil.order.M$term=row.names(H.soil.order.M)
H.soil.order.M = H.soil.order.M %>%
  mutate(cohort=i, covars="Soil_Order", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.soil.order.O = adonis2(dist.H.O ~ Soil_Order, SamDat.H.O)
H.soil.order.O$term=row.names(H.soil.order.O)
H.soil.order.O = H.soil.order.O %>%
  mutate(cohort=i, covars="Soil_Order", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.soil.order.M = adonis2(dist.WV.M ~ Soil_Order, SamDat.WV.M)
WV.soil.order.M$term = row.names(WV.soil.order.M)
WV.soil.order.M = WV.soil.order.M %>%
  mutate(cohort=i, covars="Soil_Order", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.soil.order.O = adonis2(dist.WV.O ~ Soil_Order, SamDat.WV.O)
WV.soil.order.O$term = row.names(WV.soil.order.O)
WV.soil.order.O = WV.soil.order.O %>%
  mutate(cohort=i, covars="Soil_Order", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Soil Great Group
#Butte Lake
BL.soil.greatgroup.M = adonis2(dist.BL.M ~ Soil_GreatGroup, SamDat.BL.M)
BL.soil.greatgroup.M$term=row.names(BL.soil.greatgroup.M)
BL.soil.greatgroup.M = BL.soil.greatgroup.M %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.soil.greatgroup.O = adonis2(dist.BL.O ~ Soil_GreatGroup, SamDat.BL.O)
BL.soil.greatgroup.O$term=row.names(BL.soil.greatgroup.O)
BL.soil.greatgroup.O = BL.soil.greatgroup.O %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.soil.greatgroup.M = adonis2(dist.H.M ~ Soil_GreatGroup, SamDat.H.M)
H.soil.greatgroup.M$term=row.names(H.soil.greatgroup.M)
H.soil.greatgroup.M = H.soil.greatgroup.M %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.soil.greatgroup.O = adonis2(dist.H.O ~ Soil_GreatGroup, SamDat.H.O)
H.soil.greatgroup.O$term=row.names(H.soil.greatgroup.O)
H.soil.greatgroup.O = H.soil.greatgroup.O %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.soil.greatgroup.M = adonis2(dist.WV.M ~ Soil_GreatGroup, SamDat.WV.M)
WV.soil.greatgroup.M$term = row.names(WV.soil.greatgroup.M)
WV.soil.greatgroup.M = WV.soil.greatgroup.M %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.soil.greatgroup.O = adonis2(dist.WV.O ~ Soil_GreatGroup, SamDat.WV.O)
WV.soil.greatgroup.O$term = row.names(WV.soil.greatgroup.O)
WV.soil.greatgroup.O = WV.soil.greatgroup.O %>%
  mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Veg Type
#Butte Lake
BL.veg.M = adonis2(dist.BL.M ~ Veg_Type, SamDat.BL.M)
BL.veg.M$term=row.names(BL.veg.M)
BL.veg.M = BL.veg.M %>%
  mutate(cohort=i, covars="Veg_Type", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.veg.O = adonis2(dist.BL.O ~ Veg_Type, SamDat.BL.O)
BL.veg.O$term=row.names(BL.veg.O)
BL.veg.O = BL.veg.O %>%
  mutate(cohort=i, covars="Veg_Type", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.veg.M = adonis2(dist.H.M ~ Veg_Type, SamDat.H.M)
H.veg.M$term=row.names(H.veg.M)
H.veg.M = H.veg.M %>%
  mutate(cohort=i, covars="Veg_Type", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.veg.O = adonis2(dist.H.O ~ Veg_Type, SamDat.H.O)
H.veg.O$term=row.names(H.veg.O)
H.veg.O = H.veg.O %>%
  mutate(cohort=i, covars="Veg_Type", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.veg.M = adonis2(dist.WV.M ~ Veg_Type, SamDat.WV.M)
WV.veg.M$term = row.names(WV.veg.M)
WV.veg.M = WV.veg.M %>%
  mutate(cohort=i, covars="Veg_Type", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.veg.O = adonis2(dist.WV.O ~ Veg_Type, SamDat.WV.O)
WV.veg.O$term = row.names(WV.veg.O)
WV.veg.O = WV.veg.O %>%
  mutate(cohort=i, covars="Veg_Type", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

### % fine soil
#Butte Lake
BL.percfine.M = adonis2(dist.BL.M ~ Perc_Fine, SamDat.BL.M)
BL.percfine.M$term=row.names(BL.percfine.M)
BL.percfine.M = BL.percfine.M %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.percfine.O = adonis2(dist.BL.O ~ Perc_Fine, SamDat.BL.O)
BL.percfine.O$term=row.names(BL.percfine.O)
BL.percfine.O = BL.percfine.O %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Hole
H.percfine.M = adonis2(dist.H.M ~ Perc_Fine, SamDat.H.M)
H.percfine.M$term=row.names(H.percfine.M)
H.percfine.M = H.percfine.M %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.percfine.O = adonis2(dist.H.O ~ Perc_Fine, SamDat.H.O)
H.percfine.O$term=row.names(H.percfine.O)
H.percfine.O = H.percfine.O %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

#Warner Valley
WV.percfine.M = adonis2(dist.WV.M ~ Perc_Fine, SamDat.WV.M)
WV.percfine.M$term = row.names(WV.percfine.M)
WV.percfine.M = WV.percfine.M %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.percfine.O = adonis2(dist.WV.O ~ Perc_Fine, SamDat.WV.O)
WV.percfine.O$term = row.names(WV.percfine.O)
WV.percfine.O = WV.percfine.O %>%
  mutate(cohort=i, covars="Perc_Fine", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)


#bind together outputs of PERMANOVAs
b = rbind(BL.M, BL.O, H.M, H.O, WV.M, WV.O, BL.bs.M, BL.bs.O, H.bs.M, H.bs.O, WV.bs.M, WV.bs.O, BL.quadbs.M, BL.quadbs.O, H.quadbs.M, H.quadbs.O, WV.quadbs.M, WV.quadbs.O, BL.soil.mu.M, BL.soil.mu.O, H.soil.mu.M, H.soil.mu.O, WV.soil.mu.M, WV.soil.mu.O, BL.soil.order.M, BL.soil.order.O, H.soil.order.M, H.soil.order.O, WV.soil.order.M, WV.soil.order.O, BL.soil.greatgroup.M, BL.soil.greatgroup.O, H.soil.greatgroup.M, H.soil.greatgroup.O, WV.soil.greatgroup.M, WV.soil.greatgroup.O, BL.veg.M, BL.veg.O, H.veg.M, H.veg.O, WV.veg.M, WV.veg.O, BL.percfine.M, BL.percfine.O, H.percfine.M, H.percfine.O, WV.percfine.M, WV.percfine.O)

final_df <- rbind(final_df, b, stringsAsFactors=FALSE)

#########500 Cohorts Strata=Horizon: All Covariates alone ###########

final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create phyloseq objects for each cohort
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  ps.M <- subset_samples(ps.cohort, Horizon=="M")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps.cohort))
  SamDat.M <- data.frame(sample_data(ps.M))
  
  #make dissimilarity matrix
  dist = phyloseq::distance(ps.cohort, method="bray")
  dist.M <- phyloseq::distance(ps.M, method="bray")
  
  ###only Unit
  Unit = adonis2(dist ~ Unit, SamDat, strata=SamDat$Horizon)
  Unit$term=row.names(Unit)
  Unit = Unit %>%
    mutate(cohort=i, covars="unit", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###only Trt
  Trt = adonis2(dist ~ Trt, SamDat, strata=SamDat$Horizon)
  Trt$term=row.names(Trt)
  Trt = Trt %>%
    mutate(cohort=i, covars="trt", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)

  
  ###SampleSev
  SampleSev = adonis2(dist ~ SampleSev, SamDat, strata=SamDat$Horizon)
  SampleSev$term=row.names(SampleSev)
  SampleSev = SampleSev %>%
    mutate(cohort=i, covars="SampleSev", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Avg Quad Sev
  #if(any(is.na(SamDat$Wt_Avg_QuadBS))){
  #  QuadBS = adonis2(dist ~ Wt_Avg_QuadBS, SamDat, strata=SamDat$Horizon)
  #QuadBS$term=row.names(QuadBS)
  #QuadBS = QuadBS %>%
  #  mutate(cohort=i, covars="QuadBS", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
    
  #}
  
  
  ###Soil_mu
  soil.mu = adonis2(dist ~ Soil_mu, SamDat, strata=SamDat$Horizon)
  soil.mu$term=row.names(soil.mu)
  soil.mu = soil.mu %>%
    mutate(cohort=i, covars="Soil_mu", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Order
  soil.order = adonis2(dist ~ Soil_Order, SamDat, strata=SamDat$Horizon)
  soil.order$term=row.names(soil.order)
  soil.order = soil.order %>%
    mutate(cohort=i, covars="Soil_Order", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Great Group
  soil.greatgroup = adonis2(dist ~ Soil_GreatGroup, SamDat, strata=SamDat$Horizon)
  soil.greatgroup$term=row.names(soil.greatgroup)
  soil.greatgroup = soil.greatgroup %>%
    mutate(cohort=i, covars="Soil_GreatGroup", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Subgroup
  soil.subgroup = adonis2(dist ~ Soil_GreatGroup, SamDat, strata=SamDat$Horizon)
  soil.subgroup$term=row.names(soil.subgroup)
  soil.subgroup = soil.subgroup %>%
    mutate(cohort=i, covars="Soil_Subgroup", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Veg Type
  veg = adonis2(dist ~ Veg_Type, SamDat, strata=SamDat$Horizon)
  veg$term=row.names(veg)
  veg = veg %>%
    mutate(cohort=i, covars="Veg_Type", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ### % fine soil
  #Butte Lake
  percfine = adonis2(dist.M ~ Perc_Fine, SamDat.M)
  percfine$term=row.names(percfine)
  percfine = percfine %>%
    mutate(cohort=i, covars="Perc_Fine", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  
  #bind together outputs of PERMANOVAs
  b = rbind(Unit, Trt, SampleSev, soil.mu, soil.order, soil.greatgroup, soil.subgroup, veg, percfine)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_fulldatasethorstrata_explorecovarsalone_rawresults.csv")

#cohor500_df_summary <- final_df %>%
#  filter(term=="Total") %>%
#  group_by(unit, horizon, covars) %>%
#  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_explorecovarsalone_summary.csv")

#########500 Cohorts: All Covariates alone ###########

final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create phyloseq objects for each cohort
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  ps.M <- subset_samples(ps.cohort, Horizon=="M")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps.cohort))
  SamDat.M <- data.frame(sample_data(ps.M))
  
  #make dissimilarity matrix
  dist = phyloseq::distance(ps.cohort, method="bray")
  dist.M <- phyloseq::distance(ps.M, method="bray")
  
  ###Unit
  Unit = adonis2(dist ~ Unit, SamDat)
  Unit$term=row.names(Unit)
  Unit = Unit %>%
    mutate(cohort=i, covars="unit", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###only Trt
  Trt = adonis2(dist ~ Trt, SamDat)
  Trt$term=row.names(Trt)
  Trt = Trt %>%
    mutate(cohort=i, covars="trt", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  
  ###SampleSev
  SampleSev = adonis2(dist ~ SampleSev, SamDat)
  SampleSev$term=row.names(SampleSev)
  SampleSev = SampleSev %>%
    mutate(cohort=i, covars="SampleSev", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Avg Quad Sev
  #if(any(is.na(SamDat$Wt_Avg_QuadBS))){
  #  QuadBS = adonis2(dist ~ Wt_Avg_QuadBS, SamDat, strata=SamDat$Horizon)
  #QuadBS$term=row.names(QuadBS)
  #QuadBS = QuadBS %>%
  #  mutate(cohort=i, covars="QuadBS", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #}
  
  
  ###Soil_mu
  soil.mu = adonis2(dist ~ Soil_mu, SamDat)
  soil.mu$term=row.names(soil.mu)
  soil.mu = soil.mu %>%
    mutate(cohort=i, covars="Soil_mu", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Order
  soil.order = adonis2(dist ~ Soil_Order, SamDat)
  soil.order$term=row.names(soil.order)
  soil.order = soil.order %>%
    mutate(cohort=i, covars="Soil_Order", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Great Group
  soil.greatgroup = adonis2(dist ~ Soil_GreatGroup, SamDat)
  soil.greatgroup$term=row.names(soil.greatgroup)
  soil.greatgroup = soil.greatgroup %>%
    mutate(cohort=i, covars="Soil_GreatGroup", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Soil Subgroup
  soil.subgroup = adonis2(dist ~ Soil_Subgroup, SamDat)
  soil.subgroup$term=row.names(soil.subgroup)
  soil.subgroup = soil.subgroup %>%
    mutate(cohort=i, covars="Soil_Subgroup", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Veg Type
  veg = adonis2(dist ~ Veg_Type, SamDat)
  veg$term=row.names(veg)
  veg = veg %>%
    mutate(cohort=i, covars="Veg_Type", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ### % fine soil
  percfine = adonis2(dist.M ~ Perc_Fine, SamDat.M)
  percfine$term=row.names(percfine)
  percfine = percfine %>%
    mutate(cohort=i, covars="Perc_Fine", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ### horizon
  horizon = adonis2(dist ~ Horizon, SamDat)
  horizon$term=row.names(horizon)
  horizon = horizon %>%
    mutate(cohort=i, covars="Horizon", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #bind together outputs of PERMANOVAs
  b = rbind(Unit, Trt, SampleSev, soil.mu, soil.order, soil.subgroup, soil.greatgroup, veg, percfine, horizon)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}


###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_explorecovarsalone_rawresults.csv")

#cohor500_df_summary <- final_df %>%
#  filter(term=="Total") %>%
#  group_by(covars) %>%
#  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_explorecovarsalone_summary.csv")

########## 500 cohorts subset by unit and horizon: all covars alone ###################
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  ps.BL.O <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="O")
  ps.H.O <- subset_samples(ps.cohort, Unit=="H" & Horizon=="O")
  ps.WV.O <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="O")
  ps.BL.M <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="M")
  ps.H.M <- subset_samples(ps.cohort, Unit=="H" & Horizon=="M")
  ps.WV.M <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="M")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat.BL.O = data.frame(sample_data(ps.BL.O))
  SamDat.H.O = data.frame(sample_data(ps.H.O))
  SamDat.WV.O = data.frame(sample_data(ps.WV.O))
  SamDat.BL.M = data.frame(sample_data(ps.BL.M))
  SamDat.H.M = data.frame(sample_data(ps.H.M))
  SamDat.WV.M = data.frame(sample_data(ps.WV.M))
  
  #make dissimilarity matrix
  dist.BL.O = phyloseq::distance(ps.BL.O, method="bray")
  dist.H.O = phyloseq::distance(ps.H.O, method="bray")
  dist.WV.O = phyloseq::distance(ps.WV.O, method="bray")
  dist.BL.M = phyloseq::distance(ps.BL.M, method="bray")
  dist.H.M = phyloseq::distance(ps.H.M, method="bray")
  dist.WV.M = phyloseq::distance(ps.WV.M, method="bray")
  
  ###only Trt
  #Butte Lake
  BL.M = adonis2(dist.BL.M ~ Trt, SamDat.BL.M)
  BL.M$term=row.names(BL.M)
  BL.M = BL.M %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.O = adonis2(dist.BL.O ~ Trt, SamDat.BL.O)
  BL.O$term=row.names(BL.O)
  BL.O = BL.O %>%
    mutate(cohort=i, covars="trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.M = adonis2(dist.H.M ~ Trt, SamDat.H.M)
  H.M$term=row.names(H.M)
  H.M = H.M %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.O = adonis2(dist.H.O ~ Trt, SamDat.H.O)
  H.O$term=row.names(H.O)
  H.O = H.O %>%
    mutate(cohort=i, covars="trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.M = adonis2(dist.WV.M ~ Trt, SamDat.WV.M)
  WV.M$term = row.names(WV.M)
  WV.M = WV.M %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  WV.O = adonis2(dist.WV.O ~ Trt, SamDat.WV.O)
  WV.O$term = row.names(WV.O)
  WV.O = WV.O %>%
    mutate(cohort=i, covars="trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###SampleSev
  #Butte Lake
  BL.bs.M = adonis2(dist.BL.M ~ SampleSev, SamDat.BL.M)
  BL.bs.M$term=row.names(BL.bs.M)
  BL.bs.M = BL.bs.M %>%
    mutate(cohort=i, covars="SampleSev", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.bs.O = adonis2(dist.BL.O ~ SampleSev, SamDat.BL.O)
  BL.bs.O$term=row.names(BL.bs.O)
  BL.bs.O = BL.bs.O %>%
    mutate(cohort=i, covars="SampleSev", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.bs.M = adonis2(dist.H.M ~ SampleSev, SamDat.H.M)
  H.bs.M$term=row.names(H.bs.M)
  H.bs.M = H.bs.M %>%
    mutate(cohort=i, covars="SampleSev", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.bs.O = adonis2(dist.H.O ~ SampleSev, SamDat.H.O)
  H.bs.O$term=row.names(H.bs.O)
  H.bs.O = H.bs.O %>%
    mutate(cohort=i, covars="SampleSev", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.bs.M = adonis2(dist.WV.M ~ SampleSev, SamDat.WV.M)
  WV.bs.M$term = row.names(WV.bs.M)
  WV.bs.M = WV.bs.M %>%
    mutate(cohort=i, covars="SampleSev", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  WV.bs.O = adonis2(dist.WV.O ~ SampleSev, SamDat.WV.O)
  WV.bs.O$term = row.names(WV.bs.O)
  WV.bs.O = WV.bs.O %>%
    mutate(cohort=i, covars="SampleSev", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###Avg Quad Sev
  #Butte Lake
  BL.quadbs.M = adonis2(dist.BL.M ~ Wt_Avg_QuadBS, SamDat.BL.M)
  BL.quadbs.M$term=row.names(BL.quadbs.M)
  BL.quadbs.M = BL.quadbs.M %>%
    mutate(cohort=i, covars="QuadBS", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.quadbs.O = adonis2(dist.BL.O ~ Wt_Avg_QuadBS, SamDat.BL.O)
  BL.quadbs.O$term=row.names(BL.quadbs.O)
  BL.quadbs.O = BL.quadbs.O %>%
    mutate(cohort=i, covars="QuadBS", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.quadbs.M = adonis2(dist.H.M ~ Wt_Avg_QuadBS, SamDat.H.M)
  H.quadbs.M$term=row.names(H.quadbs.M)
  H.quadbs.M = H.quadbs.M %>%
    mutate(cohort=i, covars="QuadBS", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.quadbs.O = adonis2(dist.H.O ~ Wt_Avg_QuadBS, SamDat.H.O)
  H.quadbs.O$term=row.names(H.quadbs.O)
  H.quadbs.O = H.quadbs.O %>%
    mutate(cohort=i, covars="QuadBS", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  #WV.quadbs.M = adonis2(dist.WV.M ~ Wt_Avg_QuadBS, SamDat.WV.M)
  #WV.quadbs.M$term = row.names(WV.quadbs.M)
  #WV.quadbs.M = WV.quadbs.M %>%
  #  mutate(cohort=i, covars="QuadBS", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #WV.quadbs.O = adonis2(dist.WV.O ~ Wt_Avg_QuadBS, SamDat.WV.O)
  #WV.quadbs.O$term = row.names(WV.quadbs.O)
  #WV.quadbs.O = WV.quadbs.O %>%
  #  mutate(cohort=i, covars="QuadBS", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  ###Soil_mu
  #Butte Lake
  #BL.soil.mu.M = adonis2(dist.BL.M ~ Soil_mu, SamDat.BL.M)
  #BL.soil.mu.M$term=row.names(BL.soil.mu.M)
  #BL.soil.mu.M = BL.soil.mu.M %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #BL.soil.mu.O = adonis2(dist.BL.O ~ Soil_mu, SamDat.BL.O)
  #BL.soil.mu.O$term=row.names(BL.soil.mu.O)
  #BL.soil.mu.O = BL.soil.mu.O %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #Hole
  #H.soil.mu.M = adonis2(dist.H.M ~ Soil_mu, SamDat.H.M)
  #H.soil.mu.M$term=row.names(H.soil.mu.M)
  #H.soil.mu.M = H.soil.mu.M %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="H", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #issue with number of contrasts <2 because of low sample size
  #H.soil.mu.O = adonis2(dist.H.O ~ Soil_mu, SamDat.H.O)
  #H.soil.mu.O$term=row.names(H.soil.mu.O)
  #H.soil.mu.O = H.soil.mu.O %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="H", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #Warner Valley
  #WV.soil.mu.M = adonis2(dist.WV.M ~ Soil_mu, SamDat.WV.M)
  #WV.soil.mu.M$term = row.names(WV.soil.mu.M)
  #WV.soil.mu.M = WV.soil.mu.M %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #WV.soil.mu.O = adonis2(dist.WV.O ~ Soil_mu, SamDat.WV.O)
  #WV.soil.mu.O$term = row.names(WV.soil.mu.O)
  #WV.soil.mu.O = WV.soil.mu.O %>%
  #  mutate(cohort=i, covars="Soil_mu", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  ###Soil Order
  #Butte Lake
  BL.soil.order.M = adonis2(dist.BL.M ~ Soil_Order, SamDat.BL.M)
  BL.soil.order.M$term=row.names(BL.soil.order.M)
  BL.soil.order.M = BL.soil.order.M %>%
    mutate(cohort=i, covars="Soil_Order", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.soil.order.O = adonis2(dist.BL.O ~ Soil_Order, SamDat.BL.O)
  BL.soil.order.O$term=row.names(BL.soil.order.O)
  BL.soil.order.O = BL.soil.order.O %>%
    mutate(cohort=i, covars="Soil_Order", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.soil.order.M = adonis2(dist.H.M ~ Soil_Order, SamDat.H.M)
  H.soil.order.M$term=row.names(H.soil.order.M)
  H.soil.order.M = H.soil.order.M %>%
    mutate(cohort=i, covars="Soil_Order", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  #issue with number of contrasts <2 because of low sample size
  #H.soil.order.O = adonis2(dist.H.O ~ Soil_Order, SamDat.H.O)
  #H.soil.order.O$term=row.names(H.soil.order.O)
  #H.soil.order.O = H.soil.order.O %>%
  #  mutate(cohort=i, covars="Soil_Order", horizon="O", unit="H", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #Warner Valley
  # WV.soil.order.M = adonis2(dist.WV.M ~ Soil_Order, SamDat.WV.M)
  #WV.soil.order.M$term = row.names(WV.soil.order.M)
  #WV.soil.order.M = WV.soil.order.M %>%
  #  mutate(cohort=i, covars="Soil_Order", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #WV.soil.order.O = adonis2(dist.WV.O ~ Soil_Order, SamDat.WV.O)
  #WV.soil.order.O$term = row.names(WV.soil.order.O)
  #WV.soil.order.O = WV.soil.order.O %>%
  #  mutate(cohort=i, covars="Soil_Order", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  ###Soil Great Group
  #Butte Lake
  BL.soil.greatgroup.M = adonis2(dist.BL.M ~ Soil_GreatGroup, SamDat.BL.M)
  BL.soil.greatgroup.M$term=row.names(BL.soil.greatgroup.M)
  BL.soil.greatgroup.M = BL.soil.greatgroup.M %>%
    mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.soil.greatgroup.O = adonis2(dist.BL.O ~ Soil_GreatGroup, SamDat.BL.O)
  BL.soil.greatgroup.O$term=row.names(BL.soil.greatgroup.O)
  BL.soil.greatgroup.O = BL.soil.greatgroup.O %>%
    mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.soil.greatgroup.M = adonis2(dist.H.M ~ Soil_GreatGroup, SamDat.H.M)
  H.soil.greatgroup.M$term=row.names(H.soil.greatgroup.M)
  H.soil.greatgroup.M = H.soil.greatgroup.M %>%
    mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  #issue with number of contrasts <2 because of low sample size
  #H.soil.greatgroup.O = adonis2(dist.H.O ~ Soil_GreatGroup, SamDat.H.O)
  #H.soil.greatgroup.O$term=row.names(H.soil.greatgroup.O)
  #H.soil.greatgroup.O = H.soil.greatgroup.O %>%
  #  mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="H", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #Warner Valley
  #WV.soil.greatgroup.M = adonis2(dist.WV.M ~ Soil_GreatGroup, SamDat.WV.M)
  #WV.soil.greatgroup.M$term = row.names(WV.soil.greatgroup.M)
  #WV.soil.greatgroup.M = WV.soil.greatgroup.M %>%
  #  mutate(cohort=i, covars="Soil_GreatGroup", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  #WV.soil.greatgroup.O = adonis2(dist.WV.O ~ Soil_GreatGroup, SamDat.WV.O)
  #WV.soil.greatgroup.O$term = row.names(WV.soil.greatgroup.O)
  #WV.soil.greatgroup.O = WV.soil.greatgroup.O %>%
  #  mutate(cohort=i, covars="Soil_GreatGroup", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  ###Veg Type
  #Butte Lake
  BL.veg.M = adonis2(dist.BL.M ~ Veg_Type, SamDat.BL.M)
  BL.veg.M$term=row.names(BL.veg.M)
  BL.veg.M = BL.veg.M %>%
    mutate(cohort=i, covars="Veg_Type", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.veg.O = adonis2(dist.BL.O ~ Veg_Type, SamDat.BL.O)
  BL.veg.O$term=row.names(BL.veg.O)
  BL.veg.O = BL.veg.O %>%
    mutate(cohort=i, covars="Veg_Type", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.veg.M = adonis2(dist.H.M ~ Veg_Type, SamDat.H.M)
  H.veg.M$term=row.names(H.veg.M)
  H.veg.M = H.veg.M %>%
    mutate(cohort=i, covars="Veg_Type", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  #issue with number of contrasts <2 because of low sample size
  #H.veg.O = adonis2(dist.H.O ~ Veg_Type, SamDat.H.O)
  #H.veg.O$term=row.names(H.veg.O)
  #H.veg.O = H.veg.O %>%
  #  mutate(cohort=i, covars="Veg_Type", horizon="O", unit="H", p=`Pr(>F)`) %>%
  #  select(-`Pr(>F)`)
  
  #Warner Valley
  WV.veg.M = adonis2(dist.WV.M ~ Veg_Type, SamDat.WV.M)
  WV.veg.M$term = row.names(WV.veg.M)
  WV.veg.M = WV.veg.M %>%
    mutate(cohort=i, covars="Veg_Type", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  WV.veg.O = adonis2(dist.WV.O ~ Veg_Type, SamDat.WV.O)
  WV.veg.O$term = row.names(WV.veg.O)
  WV.veg.O = WV.veg.O %>%
    mutate(cohort=i, covars="Veg_Type", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ### % fine soil
  #Butte Lake
  BL.percfine.M = adonis2(dist.BL.M ~ Perc_Fine, SamDat.BL.M)
  BL.percfine.M$term=row.names(BL.percfine.M)
  BL.percfine.M = BL.percfine.M %>%
    mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.percfine.M = adonis2(dist.H.M ~ Perc_Fine, SamDat.H.M)
  H.percfine.M$term=row.names(H.percfine.M)
  H.percfine.M = H.percfine.M %>%
    mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.percfine.M = adonis2(dist.WV.M ~ Perc_Fine, SamDat.WV.M)
  WV.percfine.M$term = row.names(WV.percfine.M)
  WV.percfine.M = WV.percfine.M %>%
    mutate(cohort=i, covars="Perc_Fine", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  
  #bind together outputs of PERMANOVAs
  b = rbind(BL.M, BL.O, H.M, H.O, WV.M, WV.O, BL.bs.M, BL.bs.O, H.bs.M, H.bs.O, WV.bs.M, WV.bs.O, BL.quadbs.M, BL.quadbs.O, H.quadbs.M, H.quadbs.O, BL.soil.order.M, BL.soil.order.O, H.soil.order.M, BL.soil.greatgroup.M, BL.soil.greatgroup.O, H.soil.greatgroup.M, BL.veg.M, BL.veg.O, H.veg.M, WV.veg.M, WV.veg.O, BL.percfine.M, H.percfine.M, WV.percfine.M)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_explorecovarsalone_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(unit, horizon, covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_explorecovarsalone_summary.csv")


########## 500 cohorts full dataset: final model ########################
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create phyloseq objects for each cohort
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps.cohort))
  
  #make dissimilarity matrix
  dist = phyloseq::distance(ps.cohort, method="bray")
  
  ###full no interaction
  full.noint = adonis2(dist ~ Horizon + Unit  + Soil_Subgroup + Veg_Type , SamDat)
  full.noint$term=row.names(full.noint)
  full.noint = full.noint %>%
    mutate(cohort=i, covars="full.noint", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###full plus interaction
  full.int = adonis2(dist ~ Horizon*Unit  + Soil_Subgroup*Unit + Veg_Type*Unit , SamDat)
  full.int$term=row.names(full.int)
  full.int = full.int %>%
    mutate(cohort=i, covars="full.int", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #bind together outputs of PERMANOVAs
  b = rbind(full.noint, full.int)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}


###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_finalmodels_rawresults.csv")

#cohor500_df_summary <- final_df %>%
#  filter(term=="Total") %>%
#  group_by(covars) %>%
#  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_finalmodels_summary.csv")

########## 500 cohorts full dataset: pH########################

final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long_pH$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long_pH[cohort_df_long_pH$Cohort==i, ]
  
  
  #create phyloseq objects for each cohort
  ps.cohort <- subset_samples(ps.pH, Full_Sample_Name %in% cohort_samples$Sample)
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps.cohort))
  
  #make dissimilarity matrix
  dist = phyloseq::distance(ps.cohort, method="bray")
  
  ###full no interaction
  adonis.ph = adonis2(dist ~ pH, SamDat)
  adonis.ph$term=row.names(adonis.ph)
  adonis.ph = adonis.ph %>%
    mutate(cohort=i, covars="pH", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)

  
  #bind together outputs of PERMANOVAs
  b = rbind(adonis.ph)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}


###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_pH_rawresults.csv")

#cohor500_df_summary <- final_df %>%
#  filter(term=="Total") %>%
#  group_by(covars) %>%
#  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_fulldataset_pH_summary.csv")


########## 500 cohorts subset by unit and horizon: pH ###################
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long_pH$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long_pH[cohort_df_long_pH$Cohort==i, ]
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps.cohort <- subset_samples(ps.pH, Full_Sample_Name %in% cohort_samples$Sample)
  ps.BL.O <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="O")
  ps.H.O <- subset_samples(ps.cohort, Unit=="H" & Horizon=="O")
  ps.WV.O <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="O")
  ps.BL.M <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="M")
  ps.H.M <- subset_samples(ps.cohort, Unit=="H" & Horizon=="M")
  ps.WV.M <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="M")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat.BL.O = data.frame(sample_data(ps.BL.O))
  SamDat.H.O = data.frame(sample_data(ps.H.O))
  SamDat.WV.O = data.frame(sample_data(ps.WV.O))
  SamDat.BL.M = data.frame(sample_data(ps.BL.M))
  SamDat.H.M = data.frame(sample_data(ps.H.M))
  SamDat.WV.M = data.frame(sample_data(ps.WV.M))
  
  #make dissimilarity matrix
  dist.BL.O = phyloseq::distance(ps.BL.O, method="bray")
  dist.H.O = phyloseq::distance(ps.H.O, method="bray")
  dist.WV.O = phyloseq::distance(ps.WV.O, method="bray")
  dist.BL.M = phyloseq::distance(ps.BL.M, method="bray")
  dist.H.M = phyloseq::distance(ps.H.M, method="bray")
  dist.WV.M = phyloseq::distance(ps.WV.M, method="bray")
  
  ###only Trt
  #Butte Lake
  BL.M = adonis2(dist.BL.M ~ pH, SamDat.BL.M)
  BL.M$term=row.names(BL.M)
  BL.M = BL.M %>%
    mutate(cohort=i, covars="pH", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.O = adonis2(dist.BL.O ~ Trt, SamDat.BL.O)
  BL.O$term=row.names(BL.O)
  BL.O = BL.O %>%
    mutate(cohort=i, covars="pH", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.M = adonis2(dist.H.M ~ pH, SamDat.H.M)
  H.M$term=row.names(H.M)
  H.M = H.M %>%
    mutate(cohort=i, covars="pH", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.O = adonis2(dist.H.O ~ Trt, SamDat.H.O)
  H.O$term=row.names(H.O)
  H.O = H.O %>%
    mutate(cohort=i, covars="pH", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.M = adonis2(dist.WV.M ~ pH, SamDat.WV.M)
  WV.M$term = row.names(WV.M)
  WV.M = WV.M %>%
    mutate(cohort=i, covars="pH", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  WV.O = adonis2(dist.WV.O ~ Trt, SamDat.WV.O)
  WV.O$term = row.names(WV.O)
  WV.O = WV.O %>%
    mutate(cohort=i, covars="pH", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  
  
  #bind together outputs of PERMANOVAs
  b = rbind(BL.M, BL.O, H.M, H.O, WV.M, WV.O)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_sepunithor_pH_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(unit, horizon, covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_sepunithor_pH_summary.csv")

################ 500 cohorts separated by unit and horizon: soil subgroup#######
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long_pH$Cohort)){
  #list of Full Name of Mineral samples for cohort 1
  cohort_samples <- cohort_df_long_pH[cohort_df_long_pH$Cohort==i, ]
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps.cohort <- subset_samples(ps.pH, Full_Sample_Name %in% cohort_samples$Sample)
  ps.BL.O <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="O")
  ps.H.O <- subset_samples(ps.cohort, Unit=="H" & Horizon=="O")
  ps.WV.O <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="O")
  ps.BL.M <- subset_samples(ps.cohort, Unit=="BL" & Horizon=="M")
  ps.H.M <- subset_samples(ps.cohort, Unit=="H" & Horizon=="M")
  ps.WV.M <- subset_samples(ps.cohort, Unit=="WV"& Horizon=="M")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat.BL.O = data.frame(sample_data(ps.BL.O))
  SamDat.H.O = data.frame(sample_data(ps.H.O))
  SamDat.WV.O = data.frame(sample_data(ps.WV.O))
  SamDat.BL.M = data.frame(sample_data(ps.BL.M))
  SamDat.H.M = data.frame(sample_data(ps.H.M))
  SamDat.WV.M = data.frame(sample_data(ps.WV.M))
  
  #make dissimilarity matrix
  dist.BL.O = phyloseq::distance(ps.BL.O, method="bray")
  dist.H.O = phyloseq::distance(ps.H.O, method="bray")
  dist.WV.O = phyloseq::distance(ps.WV.O, method="bray")
  dist.BL.M = phyloseq::distance(ps.BL.M, method="bray")
  dist.H.M = phyloseq::distance(ps.H.M, method="bray")
  dist.WV.M = phyloseq::distance(ps.WV.M, method="bray")
  
  ###only Subgroup
  #Butte Lake
  BL.M = adonis2(dist.BL.M ~ Soil_Subgroup, SamDat.BL.M)
  BL.M$term=row.names(BL.M)
  BL.M = BL.M %>%
    mutate(cohort=i, covars="soil.subgroup", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.O = adonis2(dist.BL.O ~ Soil_Subgroup, SamDat.BL.O)
  BL.O$term=row.names(BL.O)
  BL.O = BL.O %>%
    mutate(cohort=i, covars="soil.subgroup", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.M = adonis2(dist.H.M ~ Soil_Subgroup, SamDat.H.M)
  H.M$term=row.names(H.M)
  H.M = H.M %>%
    mutate(cohort=i, covars="soil.subgroup", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.O = adonis2(dist.H.O ~ Soil_subgroup, SamDat.H.O)
  H.O$term=row.names(H.O)
  H.O = H.O %>%
    mutate(cohort=i, covars="soil.subgroup", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  #can't run analysis because all samples came from same subgroup except for one
  
  
  #bind together outputs of PERMANOVAs
  b = rbind(BL.M, BL.O, H.M, H.O)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_explorecovarsalone_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(unit, horizon, covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1))

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_explorecovarsalone_summary.csv")

