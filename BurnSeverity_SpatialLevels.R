# 11/25/2023 Isabella Muscettola
# goal: to create a dataframe with all the spatial levels of measured burn severity
# (% substrate, sampled soil level, weighted avg of quad, canopy mortality, & 
#dNBR from MTBS. This will then be added to the metadata of my phyloseq object 
# ultimate goal: to incorporate burn severity into my PERMANOVA to test
# whether burn severity explains the dissimilarity between samples in my dataset).

# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

#######Burn Severity from Raw Data
soildata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/LAVO_Soil_Data_vegburn.csv")

# create column Sample_Core_Name (ex. BL_01_A), select relevant columns for burn 
#severity, replace severity with numeric value
soildata_tall <- soildata %>%
  mutate(Sample_Core_Name = paste(Plot, Quadrat, sep = "_")) %>%
  select(Sample_Core_Name, Plot, Quadrat, SampleSev, Substrate, SubFrac, Severity, SevFrac) %>%
  mutate(Severity=ifelse(Severity=='unburn', 1, ifelse(Severity=='low', 2, ifelse(Severity=='moderate', 3, ifelse(Severity=='severe', 4, ''))))) %>%
  mutate(SampleSev=ifelse(SampleSev=='unburn', 1, ifelse(SampleSev=='low', 2, ifelse(SampleSev=='moderate', 3, ifelse(SampleSev=='severe', 4, ''))))) 

#pivot substrate data so that each category is a column with % quad as value
soildata_Substrate = soildata_tall %>%
  select(Sample_Core_Name, Substrate, SubFrac) %>%
  drop_na(SubFrac) %>%
  pivot_wider(names_from = "Substrate", values_from="SubFrac")
soildata_Substrate = replace(soildata_Substrate, is.na(soildata_Substrate), 0)

#replace Substrate & SubFrac w/ LITT, BARE, ROCK, BOLE in soildata_tall df.  
soildata_tall = left_join(soildata_tall, soildata_Substrate, by = "Sample_Core_Name")
soildata_tall = select(soildata_tall, 1:4, 7:13)
#cleanup
rm(soildata_Substrate)

#calculate weighted average of the quadrat burn severity for each quadrat

#change Severity & SevFrac from character & double, respectively to numeric 
#so we can calculate a weighted average
soildata_tall$Severity = as.numeric(soildata_tall$Severity)
soildata_tall$SevFrac = as.numeric(soildata_tall$SevFrac)

#summarise quad severity as weighted average
soildata_quadsev = soildata_tall %>%
  group_by(Sample_Core_Name) %>%
  drop_na(SevFrac) %>%
  summarise(Wt_Avg_QuadBS = weighted.mean(Severity, SevFrac))
# add to soildata_tall df and replace Severity & SevFrac columns
soildata_tall = left_join(soildata_tall, soildata_quadsev, by = "Sample_Core_Name")
soildata_tall = select(soildata_tall, 1:4, 7:12)

#unique rows
soildata_final = soildata_tall[!duplicated(soildata_tall$Sample_Core_Name), ]

###Only thing left to add is the dNBR metric of burn severity for each plot

#save dataframe as "burnseverity_metadata.csv"
#write.csv(soildata_final, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/burnseverity_metadata.csv")

bs_metadata = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/burnseverity_metadata.csv")
bs_metadata = bs_metadata[, -1]
bs_metadata = bs_metadata %>%
  mutate(BARE_ROCK = BARE + ROCK)

#update metadata for phyloseq objects with the new burn severity metrics
#load old metadata 
metadata = read.csv(file ="~/Box/MSWhitman/LAVO_FireHistory/Data/R/11222023_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)

#merge dataframes according to Sample_Core_Name
df_merge = left_join(metadata, bs_metadata, by ="Sample_Core_Name")
#remove extreneous columns
df_merge =  select(df_merge, -19, -29, -30)

#save edited metadata as csv
#write.csv(df_merge, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/11252023_LAVO22_Metadata.csv")





