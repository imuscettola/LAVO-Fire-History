#Objective: 1) identify samples missing burn severity data, 2) estimate bs from 
#photos and determine correlation of field vs. photo observations, and 3) 

library(phyloseq)
library(tidyverse)
library(dplyr)
library(ggplot2)

ps.3 <- readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")

############# extract plots missing data for photo correlation testing#########
#which samples are missing burn severity? 
missingburnsev <- filter(data.frame(sample_data(ps.3)), is.na(SampleSev)|is.na(Wt_Avg_QuadBS)|is.na(BARE_ROCK)|is.na(Overstory_mort))
missingburnsev <- missingburnsev %>%
  select(Full_Sample_Name, Plot_Name, SampleSev, Wt_Avg_QuadBS, BARE_ROCK, Overstory_mort)

missingsamplesev <- filter(missingburnsev, is.na(SampleSev))

#save results
#write.csv(missingburnsev, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/MissingBurnSev.csv")

#randomly select plots for field vs. photo SampleSev estimations
corr <- filter(data.frame(sample_data(ps.3)), Unit=="BL" & SampleSev!="NA") %>%
  sample_n(., 10) %>%
  select(., Full_Sample_Name, Plot_Name, SampleSev, Wt_Avg_QuadBS, BARE_ROCK, Overstory_mort)

#write.csv(corr, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/temp.csv")

corr <- rbind(missingsamplesev, corr)

#write.csv(corr, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/SampleSevCorrTest.csv")

##############Correlation of field vs. photo sample burn severity estimates#####
#How well was I able to estimate sample burn severity from photos? 
corr <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/SampleSevCorrTest.csv")
corr.temp <- corr %>%
  filter(SampleSev==SampleSev_photo) %>%
  summarise(n=n()) #14/20 had exact match of field v photo estimates
rm(corr.temp)

ggplot(corr[corr$SampleSev!="NA", ], aes(x=SampleSev, y=SampleSev_photo)) + geom_point() + geom_abline()
#linear regression
summary(lm(SampleSev_photo~SampleSev, data=corr[corr$SampleSev!="NA", ])) #decent correlation!
#y=0.259 + 0.99442x

############ Update the sample data for the samples missing SampleSev##########
missingsamplesev <- filter(data.frame(sample_data(ps.3)), is.na(SampleSev))

# load data and prepare for joining dataframes
corr <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/SampleSevCorrTest.csv")
corr <- subset(corr, Full_Sample_Name %in% missingsamplesev$Full_Sample_Name) %>%
  select(. , Full_Sample_Name, SampleSev_photo)
  
#join missingsamplesev and tidied corr dataframes
add <- left_join(missingsamplesev, corr, by="Full_Sample_Name")

#replace SampleSev with SampleSev_photo values
add <- add %>%
  select(!SampleSev)%>%
  mutate(SampleSev=SampleSev_photo)%>%
  select(!SampleSev_photo)

#add previously missing samples back into sample_data phyloseq object
notmissing <- filter(data.frame(sample_data(ps.3)), SampleSev!="NA")
SamDat <- rbind(notmissing, add)

#replace current sample data in ps.3 object with tidied data
otu <- otu_table(ps.3)
tax <- tax_table(ps.3)

ps.3.edit <- phyloseq(otu, tax, SamDat)

saveRDS(ps.3.edit, "~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3.02122024")
############################### DIDNT WORK ####################################
#this didn't work -> will re-make the phyloseq objects from the .biom files and 
#updated metadata csv.

#add the sample sev to the metadata csv; most recent and previously used to make phyloseq objects: 11252023_LAVO22_Metadata.csv
#load old metadata
oldmetadata <- read.csv(file ="~/Box/MSWhitman/LAVO_FireHistory/Data/R/11252023_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE) 
missingsamplesev <- filter(oldmetadata, is.na(SampleSev) & Unit=="BL")

# load data and prepare for joining dataframes
corr <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/SampleSevCorrTest.csv")
corr <- subset(corr, Full_Sample_Name %in% missingsamplesev$Full_Sample_Name) %>%
  select(. , Full_Sample_Name, SampleSev_photo)

#join missingsamplesev and tidied corr dataframes
add <- left_join(missingsamplesev, corr, by="Full_Sample_Name")

#replace SampleSev with SampleSev_photo values
add <- add %>%
  select(!SampleSev)%>%
  mutate(SampleSev=SampleSev_photo)%>%
  select(!SampleSev_photo)

#add previously missing samples back into sample_data phyloseq object
notmissing <- setdiff(oldmetadata, missingsamplesev)
SamDat <- rbind(notmissing, add)

#save output as new metadata 
write.csv(SamDat, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/02122024_LAVO22_Metadata.csv")

################## ADD MTBS data to metadata ##################################
#add mtbs sev to the metadata csv; most recent and previously used to make phyloseq objects: 

mtbs <- read.csv("~/Desktop/plots_22_mtbs.csv")
#mtbs data is "gridcode" column

oldmetadata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/02122024_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors = TRUE)

#join gridcode to oldmetadata
newmetadata <- merge(oldmetadata, mtbs[, c("Plot", 'gridcode')], by.x="Plot_Name", by.y='Plot', all.x=TRUE)
colnames(newmetadata)[colnames(newmetadata)=='gridcode'] <- 'mtbs_sev'
#save as new metadata
write.csv(newmetadata, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/02152024_LAVO22_Metadata.csv")










