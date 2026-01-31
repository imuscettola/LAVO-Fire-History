# 11/14/2023 Isabella Muscettola
# goal: to add median burn severity for each quadrat 
# to the metadata of my phyloseq object 
# (ultimate goal:
# to incorporate burn severity into my PERMANOVA to test
# whether burn severity explains the dissimilarity 
# between samples in my dataset).

# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

#######Calculating median burn severity for each quadrat (each replicate)#########
# read csv of veg and burn from field data 
soildata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/LAVO_Soil_Data_vegburn.csv")

# create column Sample_Core_Name (ex. BL_01_A), select relevant columns for burn severity at the quadrat level, remove rows with NA SevFrac values, replace severity with numeric value
soildata_tall <- soildata %>%
  mutate(Sample_Core_Name = paste(Plot, Quadrat, sep = "_")) %>%
  select(Sample_Core_Name, Plot, Quadrat, Severity, SevFrac) %>%
  filter(SevFrac>0) %>%
  mutate(Severity=ifelse(Severity=='unburn', 1, ifelse(Severity=='low', 2, ifelse(Severity=='moderate', 3, ifelse(Severity=='severe', 4, ''))))) 

#converting the Severity column from character to numeric seems to be necessary 
#for calculating median Severity
soildata_tall$Severity <- as.numeric(soildata_tall$Severity)

#calculate median severity for each quadrat.
soildata_bs <- soildata_tall %>%
  mutate(SevFrac=SevFrac*100) %>%
  group_by(Sample_Core_Name) %>%
  summarise(median_bs = median(rep(Severity, SevFrac)))
# 11/14/23 - this way of calculating the median results in 0.5 values for quads
# where it was 50/50 of two severities; is this appropriate? 

#save the dataframe as "MedianBS_Quadrat"
#write.csv(soildata_bs, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/BurnSeverity/MedianBS_Quadrat.csv")

#read dataframe created above
#bs <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/BurnSeverity/MedianBS_Quadrat.csv")

######add the median burn severity from this dataframe to metadata###########
metadata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/28042023_LAVO22_Metadata.csv")


