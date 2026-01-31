# Objective: to summarize findings of differential abundance tests
#1) How many OTUs identified with different models in fire history analysis; 
#and which OTUs are identified in multiple analyses.
#2) ID OTUs that are present across units for each treatment (a) burn severity 
#and b) fire history)
#3) ID OTUs that are increasingly abundant with increasing burn severity
#4) ID OTUs that are differnetially abundant across both fire history & burn 
#severity treatments
#5) plot OTUs that are dif abundant with a) fire history and b) burn sev for each unit
#6) compare to OTUs ID as fire responsive from Dana's paper -> requires going 
#back to FASTA files and BLASTing against my fasta files. 
#

# load libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(corncob)

#load dataframes 
BL_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob_BL_BSfac_DA.csv")
H_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob_Hole_BSfac_DA.csv")
WV_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob_WV_BSfac_DA.csv")

BL_trt_all <- read.csv("~/Desktop/LAVO22_GG2/corncob_BLDA.csv")
H_trt_all <- read.csv("~/Desktop/LAVO22_GG2/corncob_HoleDA.csv")
WV_trt_all <- read.csv("~/Desktop/LAVO22_GG2/corncob_WVDA.csv")

#filter the unit dataframes for horizon as the covariate
BL_trt <- filter(BL_trt_all, covars=="horizon")
H_trt <- filter(H_trt_all, covars=="horizon")
WV_trt <- filter(WV_trt_all, covars=="horizon")

#extract OTU # and barcode from phyloseq object to reference against ID'ed fire 
#responders in Johnson et al., 2023
ps.3 = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")
ps.3 <- clean_taxa_names(ps.3)

tax_table(ps.3)[1:5, ]
tax.table[1:5, ]
tax.table <- as.data.frame(tax_table(ps.3))

#clean up taxa names so OTU does not start with a number 


###### burn severity dataframes from wide to long ####
df_bs_BL_H <- rbind(BL_bs, H_bs) %>%
  mutate(mu.unburned = 0) %>%
  select(!1)
df_bs_BL_H$X <- row.names(df_bs_BL_H)
#Butte Lake & Hole 
temp1 <- df_bs_BL_H %>%
  mutate(unburned=0) %>%
  select(X, unburned, 3:5) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="mu") %>%
  mutate(comp=ifelse(comparison=="unburned", "unburned", ifelse(comparison=="mu.low", "low", ifelse(comparison=="mu.moderate", "moderate", "severe")))) %>%
  select(1, 3:4) 

temp2 <- df_bs_BL_H %>%
  select(X, 6:8) %>%
  pivot_longer(cols=!1, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_low", "low", ifelse(comparison=="t.test_p_moderate", "moderate", "severe")))%>%
  select(1, 3:4) 


temp_merge <- left_join(temp1, temp2, by = c('X', 'comp'))

taxonomy <- df_bs_BL_H %>%
  select(X, OTU, unit, Domain:Species)

df_bs_BL_H_long <- left_join(temp_merge, taxonomy, by='X') %>%
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  mutate(., `p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE")) %>%
  ungroup(.) %>%
  arrange(Phylum, mu)
y_ord <- unique(df_bs_BL_H_long$Family_Genus_OTU)
df_bs_BL_H_long$Family_Genus_OTU <- factor(df_bs_BL_H_long$Family_Genus_OTU, levels=y_ord)

rm(temp1, temp2, temp_merge, taxonomy)


#Warner Valley
temp1 <- WV_bs %>%
  mutate(low=0) %>%
  select(X, 4:5, low) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="mu") %>%
  mutate(comp=ifelse(comparison=="low", "low", ifelse(comparison=="mu.lowtomod", "moderate", "severe"))) %>%
  select(1, 3:4) 

temp2 <- WV_bs %>%
  select(X, 6:7) %>%
  pivot_longer(cols=!1, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_moderate", "moderate", "severe")) %>%
  select(1, 3:4) 


temp_merge <- left_join(temp1, temp2, by = c('X', 'comp'))

taxonomy <- WV_bs %>%
  select(X, OTU, Domain:Species)

temp <- left_join(temp_merge, taxonomy, by='X') %>%
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) 
WV_bs_long <- temp %>%
  ungroup() %>%
  arrange(Phylum) %>%
  mutate(`p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE"))
y_ord <- unique(WV_bs_long$Family_Genus_OTU)
WV_bs_long$Family_Genus_OTU <- factor(WV_bs_long$Family_Genus_OTU, levels=y_ord)

rm(temp, temp1, temp2, temp_merge, taxonomy)


###### 1) which OTUs are ID'ed with different covariates?  ########

###Butte Lake
temp <- BL_trt_all[BL_trt_all$covars=="none" | BL_trt_all$covars=="horizon",]
BL_novhor <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
794/2
#397 OTUs exist in both the differential tests (when controlling for horizon and not controlling for it....)

temp <- BL_trt_all[BL_trt_all$covars=="none" | BL_trt_all$covars=="O horizon only",]
BL_novO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
490/2
#245 

temp <- BL_trt_all[BL_trt_all$covars=="none" | BL_trt_all$covars=="M horizon only",]
BL_novM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
626/2
#313

temp <- BL_trt_all[BL_trt_all$covars=="horizon" | BL_trt_all$covars=="O horizon only",]
BL_horvO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
478/2
#239

temp <- BL_trt_all[BL_trt_all$covars=="horizon" | BL_trt_all$covars=="M horizon only",]
BL_horvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
594/2
#297

temp <- BL_trt_all[BL_trt_all$covars=="M horizon only" | BL_trt_all$covars=="O horizon only",]
BL_OvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
298/2
#149

#which exist in all 4?
BL_all <- BL_trt_all[unlist(tapply(1:nrow(BL_trt_all), BL_trt_all$OTU, function(x) if(length(x)>3) x)), ]
588/4
#147


### Hole
temp <- H_trt_all[H_trt_all$covars=="none" | H_trt_all$covars=="horizon",]
H_novhor <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
62/2
#31 OTUs exist in both the differential tests (when controlling for horizon and not controlling for it....)

temp <- H_trt_all[H_trt_all$covars=="none" | H_trt_all$covars=="O horizon only",]
H_novO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
28/2
#14 

temp <- H_trt_all[H_trt_all$covars=="none" | H_trt_all$covars=="M horizon only",]
H_novM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
26/2
#13

temp <- H_trt_all[H_trt_all$covars=="horizon" | H_trt_all$covars=="O horizon only",]
H_horvO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
16/2
#8

temp <- H_trt_all[H_trt_all$covars=="horizon" | H_trt_all$covars=="M horizon only",]
H_horvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
26/2
#13

temp <- H_trt_all[H_trt_all$covars=="M horizon only" | H_trt_all$covars=="O horizon only",]
H_OvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
#0

#which exist in all 4 -> 0


### Warner Valley
temp <- WV_trt_all[WV_trt_all$covars=="none" | WV_trt_all$covars=="horizon",]
WV_novhor <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
144/2
#72 OTUs exist in both the differential tests (when controlling for horizon and not controlling for it....)

temp <- WV_trt_all[WV_trt_all$covars=="none" | WV_trt_all$covars=="O horizon only",]
WV_novO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
70/2
#35 

temp <- WV_trt_all[WV_trt_all$covars=="none" | WV_trt_all$covars=="M horizon only",]
WV_novM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
90/2
#45

temp <- WV_trt_all[WV_trt_all$covars=="horizon" | WV_trt_all$covars=="O horizon only",]
WV_horvO <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
78/2
#39

temp <- WV_trt_all[WV_trt_all$covars=="horizon" | WV_trt_all$covars=="M horizon only",]
WV_horvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
98/2
#49

temp <- WV_trt_all[WV_trt_all$covars=="M horizon only" | WV_trt_all$covars=="O horizon only",]
WV_OvM <- temp[unlist(tapply(1:nrow(temp), temp$OTU, function(x) if(length(x)>1) x)), ] 
24/2
#12

#which exist in all 4?
WV_all <- WV_trt_all[unlist(tapply(1:nrow(WV_trt_all), WV_trt_all$OTU, function(x) if(length(x)>3) x)), ]
48/4


########### 2) ID OTUs across units ##############
#moving forward we will be using models that include horizon as a covariate. 

### a) fire history 
#merge dataframes together
df_trt <- rbind(BL_trt, H_trt, WV_trt)
#how many total unique OTUs identified as differentially abundant? 
temp <- df_trt %>%
  distinct(OTU, .keep_all=TRUE)
length(temp$OTU) #543 total OTUs ID as differentially abundant

#found in all 3 units
all_units_trt <- df_trt[unlist(tapply(1:nrow(df_trt), df_trt$OTU, function(x) if(length(x)>2) x)), ] 
#only 1 OTU was identified as differentially abundant across all 3 units

#found in at least 2 units
two_units_trt <- df_trt[unlist(tapply(1:nrow(df_trt), df_trt$OTU, function(x) if(length(x)>1) x)), ] 
#101 -> 50 OTUs identified as differnetially abundant across at least 2 units


### b) burn severity 
#bind together all 3 units bs data
temp1 <- df_bs_BL_H %>%
  select(OTU, p_fdr, comparison, mu.low:t.test_p_severe, unit, Domain:Species)

temp2 <- WV_bs %>%
  mutate(mu.low=0, mu.moderate=mu.lowtomod, mu.severe=mu.lowtosevere, t.test_p_low='NA') %>%
  select(OTU, p_fdr, comparison, mu.low, mu.moderate, mu.severe, t.test_p_low, t.test_p_moderate, t.test_p_severe, unit, Domain:Species)

temp_merge <- rbind(temp1, temp2)

#how many unique OTUs exist with burn severity?
length(unique(temp_merge$OTU)) #673 total
length(BL_bs$OTU) #351 in Butte Lake
length(H_bs$OTU) #319 in Hole
length(WV_bs$OTU) #198 in Warner Valley

#which OTUs are diff. abundant across the units? 
# across all 3 units
all_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>2) x)), ]
length(all_units_bs$OTU) #23 OTUs 

#found in at least 2 units
two_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>1) x)), ] 
length(two_units_bs$OTU) #367-69=298/2=149 OTUs diff abundant across only 2 units


##### 3) OTUs linearly diff abundant with increasing burn sev (or decreasing) ####


##### 4) OTUS ID with both fire history and burn severity ######

##### 5) Plot OTUs and log2 fold change for a) fire history and b) burn sev for each unit ####
### a) fire history ####
##Butte Lake 
BL_trt <- BL_trt %>%
  mutate(`p<0.05`=ifelse(p_fdr>0.05, 'FALSE', 'TRUE')) %>%
  arrange(mu.MF.trtmt)

temp <- BL_trt %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, mu.MF.trtmt)
y_ord <- temp$Family_Genus_OTU
temp$Family_Genus_OTU <- factor(temp$Family_Genus_OTU, y_ord)

a <- ggplot(temp[1:108, ], aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")
a

b <- ggplot(temp[109:216, ], aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")

c <- ggplot(temp[217:324, ], aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")

d <- ggplot(temp[325:432, ], aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")

e <- ggplot(temp, aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")


a+b
c+d
e

##Hole
temp <- H_trt %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, mu.MF.trtmt)
y_ord <- temp$Family_Genus_OTU
temp$Family_Genus_OTU <- factor(temp$Family_Genus_OTU, y_ord)

f <- ggplot(temp, aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")

f

##Warner Valley 
temp <- WV_trt %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, mu.MF.trtmt)
y_ord <- temp$Family_Genus_OTU
temp$Family_Genus_OTU <- factor(temp$Family_Genus_OTU, y_ord)

g <- ggplot(temp, aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change")

g

###ID'ed across all 3 units with fire history
h <- ggplot(all_units_trt, aes(x=mu.MF.trtmt, y=OTU, color=Phylum)) + geom_point(aes(shape=unit)) + labs(y="OTU", x= "log2 fold change") 
h


###ID'ed at least 2 units with fire history 
temp <- two_units_trt %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", ifelse(is.na(Family), "unclassified", Family))) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum)
y_ord <- temp$Family_Genus_OTU
temp$Family_Genus_OTU <- factor(temp$Family_Genus_OTU, y_ord)

i <- ggplot(temp, aes(x=mu.MF.trtmt, y=Family_Genus_OTU, color=Phylum)) + geom_point(aes(shape=unit)) + labs(y="OTU", x= "log2 fold change") 
i





### b) burn severity ####
# need to transform dataframe from wide to long
###Butte Lake & Hole together
length(df_bs_BL_H_long$OTU)
2010/4
m <- ggplot(df_bs_BL_H_long, aes(x=mu, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + geom_point(aes(shape=comp)) +geom_vline(xintercept=0)+ labs(y="OTU", x="log2 fold change", title="Butte Lake & Hole OTUs differentially abundant with burn severity") + facet_wrap(~unit)
m
#[df_bs_BL_H_long$Phylum=="Acidobacteriota", ]

### Warner Valley alone
length(WV_bs_long$OTU)
594/2
k <- ggplot(WV_bs_long[1:297, ], aes(x=mu, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + geom_point(aes(shape=comp)) + labs(y="OTU", x= "log2 fold change", title="Warner Valley OTUs differentially abundant with burn severity (1)") 
k

l <- ggplot(WV_bs_long[298:594, ], aes(x=mu, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + geom_point(aes(shape=comp)) + labs(y="OTU", x= "log2 fold change", title="Warner Valley OTUs differentially abundant with burn severity (2)") 
l

k+l

### 





##### 6) 