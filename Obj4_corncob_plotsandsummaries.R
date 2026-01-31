# Objective: to summarize findings of differential abundance tests
#1) ID OTUs that are present across units for each treatment (a) burn severity 
#and b) fire history)
#2) plot OTUs that are dif abundant with a) fire history and b) burn sev for each unit


# load libraries
library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)
library(corncob)
library(viridis)

#load dataframes 
BL_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_BL_BSfac_DA.csv")
H_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_Hole_BSfac_DA.csv")
WV_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_WV_BSfac_DA.csv")

BL_trt <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_BLDA.csv") %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, L2FC)
y_ord <- BL_trt$Family
BL_trt$Family <- factor(BL_trt$Family, unique(y_ord))

H_trt <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_Hole_DA.csv") %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, L2FC)
y_ord <- H_trt$Family
H_trt$Family <- factor(H_trt$Family, unique(y_ord))

WV_trt <- read.csv("~/Desktop/LAVO22_GG2/corncob/corncob_WVDA.csv") %>%
  mutate(Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup() %>%
  arrange(Phylum, L2FC)
y_ord <- WV_trt$Family
WV_trt$Family <- factor(WV_trt$Family, unique(y_ord))


#### Count OTUs in each unit by Phylum separating by sign of log2fold change
df_trt_count <- df_trt %>%
  mutate(sign=ifelse(L2FC>0, 'pos', 'neg')) %>%
  group_by(unit, Phylum, sign) %>%
  summarise(n=n())

#write.csv(df_trt_count, '~/Desktop/LAVO22_GG2/corncob/firehistory_count_summary.csv')
  

#### ID Fire responders as OTUs with significant + l2fc across bs types ###
BL_firerespond <- BL_bs 
  
  


###### burn severity dataframes from wide to long ####

#Butte Lake  
temp1 <- BL_bs_filter %>%
  mutate(unburned=0) %>%
  select(X, unburned, 23:25) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="L2FC") %>%
  mutate(comp=ifelse(comparison=="unburned", "unburned", ifelse(comparison=="L2FC.low", "low", ifelse(comparison=="L2FC.moderate", "moderate", "severe")))) %>%
  select(1, 3:4) 

temp2 <- BL_bs_filter %>%
  select(X, 8:10) %>%
  pivot_longer(cols=!1, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_low", "low", ifelse(comparison=="t.test_p_moderate", "moderate", "severe")))%>%
  select(1, 3:4) 

temp3 <- BL_bs_filter %>%
  select(X, Rel.Abund.unburned:Rel.Abund.severe) %>%
  pivot_longer(cols=!1, names_to='comparison', values_to='Rel.Abund') %>%
  mutate(comp=ifelse(comparison=="Rel.Abund.unburned", 'unburned', ifelse(comparison=='Rel.Abund.low', 'low', ifelse(comparison=='Rel.Abund.moderate', 'moderate', 'severe')))) %>%
  select(1, 3:4)

temp_merge <- temp1 %>%
  left_join(temp2, by = c('X', 'comp')) %>%
  left_join(temp3, by = c('X', 'comp'))

taxonomy <- BL_bs_filter %>%
  select(X, OTU, unit, Domain:Species)

BL_bs_long <- left_join(temp_merge, taxonomy, by='X') %>%
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  mutate(., `p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE")) %>%
  mutate(., a = ifelse(t.test_p<0.05, 0.8, 0.3)) %>%
  ungroup(.) %>%
  arrange(Phylum, L2FC)
y_ord <- unique(BL_bs_long$Family_Genus_OTU)
BL_bs_long$Family_Genus_OTU <- factor(BL_bs_long$Family_Genus_OTU, levels=y_ord)

rm(temp1, temp2, temp_merge, taxonomy)

#Hole  
temp1 <- H_bs_filter %>%
  mutate(unburned=0) %>%
  select(X, unburned, 19:20) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="L2FC") %>%
  mutate(comp=ifelse(comparison=="unburned", "unburned", ifelse(comparison=="L2FC.low", "low", ifelse(comparison=="L2FC.moderate", "moderate", "severe")))) %>%
  select(1, 3:4) 

temp2 <- H_bs_filter %>%
  select(X, 7:8) %>%
  pivot_longer(cols=!1, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_low", "low", ifelse(comparison=="t.test_p_moderate", "moderate", "severe")))%>%
  select(1, 3:4) 

temp3 <- H_bs_filter %>%
  select(X, Rel.Abund.unburned:Rel.Abund.moderate) %>%
  pivot_longer(cols=!1, names_to='comparison', values_to='Rel.Abund') %>%
  mutate(comp=ifelse(comparison=="Rel.Abund.unburned", 'unburned', ifelse(comparison=='Rel.Abund.low', 'low', 'moderate'))) %>%
  select(1, 3:4)

temp_merge <- temp1 %>%
  left_join(temp2, by = c('X', 'comp')) %>%
  left_join(temp3, by = c('X', 'comp'))

taxonomy <- H_bs_filter %>%
  select(X, OTU, unit, Domain:Species)

H_bs_long <- left_join(temp_merge, taxonomy, by='X') %>%
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  mutate(., `p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE")) %>%
  mutate(., a = ifelse(t.test_p<0.05, 0.8, 0.3)) %>%
  ungroup(.) %>%
  arrange(Phylum, L2FC)
y_ord <- unique(H_bs_long$Family_Genus_OTU)
H_bs_long$Family_Genus_OTU <- factor(H_bs_long$Family_Genus_OTU, levels=y_ord)

rm(temp1, temp2, temp_merge, taxonomy)

#Warner Valley
temp1 <- WV_bs_filter %>%
  mutate(low=0) %>%
  select(X, 19:20, low) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="L2FC") %>%
  mutate(comp=ifelse(comparison=="low", "low", ifelse(comparison=="L2FC.mod.v.low", "moderate", "severe"))) %>%
  select(1, 3:4) 

temp2 <- WV_bs_filter %>%
  select(X, 7:8) %>%
  pivot_longer(cols=!1, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_moderate", "moderate", "severe")) %>%
  select(1, 3:4) 

temp3 <- WV_bs_filter %>%
  select(X, Rel.Abund.low:Rel.Abund.severe) %>%
  pivot_longer(cols=!1, names_to='comparison', values_to='Rel.Abund') %>%
  mutate(comp=ifelse(comparison=='Rel.Abund.low', 'low', ifelse(comparison=='Rel.Abund.moderate', 'moderate', 'severe'))) %>%
  select(1, 3:4)

temp_merge <- temp1 %>%
  left_join(temp2, by = c('X', 'comp')) %>%
  left_join(temp3, by = c('X', 'comp'))

taxonomy <- WV_bs_filter %>%
  select(X, OTU, unit, Domain:Species)

WV_bs_long <- left_join(temp_merge, taxonomy, by='X') %>% 
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  ungroup(.) %>%
  arrange(., Phylum, L2FC) %>%
  mutate(., `p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE")) %>%
  mutate(., a = ifelse(t.test_p<0.05, 0.8, 0.3))
y_ord <- unique(WV_bs_long$Family_Genus_OTU)
WV_bs_long$Family_Genus_OTU <- factor(WV_bs_long$Family_Genus_OTU, levels=y_ord)

rm(temp3, temp1, temp2, temp_merge, taxonomy)

#bind together long dataframes to 1
df_bs_long <- rbind(BL_bs_long, H_bs_long, WV_bs_long)

#summarise OTU table
df_bs_long_filter <- df_bs_long %>%
  filter(t.test_p < 0.05) %>%
  mutate(sign=ifelse(L2FC>0, 'pos', 'neg')) %>%
  group_by(unit, Phylum, sign) %>%
  summarise(n=n())

write.csv(df_bs_long_filter, '~/Desktop/LAVO22_GG2/corncob/bs_count_summary.csv')


########### 1) ID OTUs across units ##############
#moving forward we will be using models that include horizon as a covariate. 

### a) fire history 
#merge dataframes together
df_trt <- rbind(BL_trt, H_trt, WV_trt)

#write.csv(df_trt, "~/Desktop/LAVO22_GG2/corncob/all_firehistory.csv")

library(stringr)



fivenum(df_trt$Rel.Abund.MF)
fivenum(temp_merge$Rel.Abund.unburned)
fivenum(temp_merge$Rel.Abund.low)

Nova_trt <- df_trt %>%
  filter(str_detect(Genus, "^Noviherbaspirillum")) %>%
  filter(abs(L2FC)>2)
write.csv(Nova_trt, "~/Desktop/LAVO22_GG2/corncob/Noviherbaspirillum_trt.csv")

Nova_bs <- temp_merge %>%
  filter(str_detect(Genus, "^Noviherbaspirillum")) %>%
  filter(abs(L2FC)>2)
write.csv(Nova_bs, "~/Desktop/LAVO22_GG2/corncob/Noviherbaspirillum_bs.csv")

Massilia_trt <- filter(df_trt, Genus=='Massilia'& abs(L2FC)> 2)
Massilia_bs <- filter(temp_merge, Genus=='Massilia')
write.csv(Massilia_trt, "~/Desktop/LAVO22_GG2/corncob/Massilia_trt.csv")
write.csv(Massilia_bs, "~/Desktop/LAVO22_GG2/corncob/Massilia_bs.csv")

Mod_trt <- filter(df_trt, Genus=="Modestobacter" & abs(L2FC)> 2)
Mod_bs <-filter(temp_merge, Genus=="Modestobacter")
write.csv(Mod_trt, "~/Desktop/LAVO22_GG2/corncob/Modestobacter_trt.csv")
write.csv(Mod_bs, "~/Desktop/LAVO22_GG2/corncob/Modestobacter_bs.csv")

#how many total unique OTUs identified as differentially abundant? 
temp <- df_trt %>%
  distinct(OTU, .keep_all=TRUE)
length(temp$OTU) #543 total OTUs ID as differentially abundant

#found in all 3 units
all_units_trt <- df_trt[unlist(tapply(1:nrow(df_trt), df_trt$OTU, function(x) if(length(x)>2) x)), ] 
length(all_units_trt$OTU)
#only 1 OTU was identified as differentially abundant across all 3 units

ggplot(all_units_trt, aes(x=L2FC, y=Family_Genus_OTU)) + geom_point(aes(shape=unit, color=Phylum, size=Rel.Abund.MF)) +
  labs(title='OTU differentially abundant across all three units', y='Family, Genus, OTU', x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#800000')) +
  scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08))

#found in at least 2 units
two_units_trt <- df_trt[unlist(tapply(1:nrow(df_trt), df_trt$OTU, function(x) if(length(x)>1) x)), ] %>%
  arrange(., Phylum, L2FC)
length(two_units_trt$OTU) #101 -> 50 OTUs identified as differnetially abundant across at least 2 units

#save output
write.csv(two_units_trt, "~/Desktop/LAVO22_GG2/corncob/twounits_firehistory.csv")

#plot OTUs 
y_ord <- unique(two_units_trt$Family)
two_units_trt$Family <- factor(two_units_trt$Family, levels=y_ord)

ggplot(two_units_trt, aes(x=L2FC, y=reorder(Family, desc(Family)), color=Phylum)) + 
  geom_point(aes(shape=unit, size=Rel.Abund.MF)) + theme_bw() +
  labs(title='OTUs differentially abundant across at least two units', y='Family', size='More Fire Avg Relative Abundance', x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)', shape='Unit') + 
  scale_color_manual(values=c('#e6194B','#f58231','#bfef45', '#800000','#808000','#000075','#000000')) + 
  scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + geom_vline(xintercept=0, aes(alpha=0.7)) + xlim(c(-10,10))


### b) burn severity 
#bind together all 3 units bs data
temp1 <- BL_bs %>%
  select(OTU, p_fdr, comparison, L2FC.low:L2FC.severe, t.test_p_low:t.test_p_severe, unit, Domain:Species)

temp2 <- H_bs %>%
  mutate(L2FC.severe=0, t.test_p_severe='NA') %>%
  select(OTU, p_fdr, comparison, L2FC.low:L2FC.moderate, L2FC.severe, t.test_p_low:t.test_p_moderate, t.test_p_severe, unit, Domain:Species)

temp3 <- WV_bs %>%
  mutate(L2FC.low=0, L2FC.moderate=L2FC.mod.v.low, L2FC.severe=L2FC.sev.v.low, t.test_p_low='NA') %>%
  select(OTU, p_fdr, comparison, L2FC.low, L2FC.moderate:L2FC.severe, t.test_p_low, t.test_p_moderate:t.test_p_severe, unit, Domain:Species)

temp_merge <- rbind(temp1, temp2, temp3)

#how many unique OTUs exist with burn severity?
length(unique(temp_merge$OTU)) #619 total
length(BL_bs$OTU) #351 in Butte Lake
length(H_bs$OTU) #225 in Hole
length(WV_bs$OTU) #198 in Warner Valley

#which OTUs are diff. abundant across the units? 
# across all 3 units
all_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>2) x)), ]
length(all_units_bs$OTU) #54/3 = 18 OTUs 


#found in at least 2 units
two_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>1) x)), ] 
length(two_units_bs$OTU) #292-54=238/2=119 OTUs diff abundant across only 2 units


#### focus on OTUs with a significant contrast to the lowest burn severity classification
##filter out OTUs that do not have significant contrasts to lowest burn severity (unburn or low)
BL_bs_filter <- BL_bs %>%
  filter(t.test_p_low < 0.05 | t.test_p_moderate < 0.05 | t.test_p_severe < 0.05)

H_bs_filter <- H_bs %>%
  filter(t.test_p_low < 0.05 | t.test_p_moderate < 0.05)

WV_bs_filter <- WV_bs %>%
  filter(t.test_p_moderate < 0.05 | t.test_p_severe < 0.05)

length(BL_bs_filter$OTU) #305
length(H_bs_filter$OTU) #212
length(WV_bs_filter$OTU) #191

#bind together all 3 units bs data
temp1 <- BL_bs_filter %>%
select(OTU, p_fdr, comparison, L2FC.low:L2FC.severe, Rel.Abund.unburned:Rel.Abund.severe, t.test_p_low:t.test_p_severe, unit, Domain:Species)

temp2 <- H_bs_filter %>%
  mutate(L2FC.severe=0, Rel.Abund.severe=0, t.test_p_severe='NA') %>%
  select(OTU, p_fdr, comparison, L2FC.low:L2FC.moderate, L2FC.severe, Rel.Abund.unburned:Rel.Abund.moderate, Rel.Abund.severe, t.test_p_low:t.test_p_moderate, t.test_p_severe, unit, Domain:Species)

temp3 <- WV_bs_filter %>%
  mutate(L2FC.low=0, L2FC.moderate=L2FC.mod.v.low, L2FC.severe=L2FC.sev.v.low, t.test_p_low='NA', Rel.Abund.unburned=0) %>%
  select(OTU, p_fdr, comparison, L2FC.low, L2FC.moderate:L2FC.severe, Rel.Abund.unburned, Rel.Abund.low:Rel.Abund.severe, t.test_p_low, t.test_p_moderate:t.test_p_severe, unit, Domain:Species)

temp_merge <- rbind(temp1, temp2, temp3)

#total number of unique OTUs identified with sig contrast to lowest burn severity 
length(unique(temp_merge$OTU)) #565

write.csv(temp_merge, "~/Desktop/LAVO22_GG2/corncob/all_burnseverity.csv")

#which OTUs are diff. abundant across the units? 
# across all 3 units
all_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>2) x)), ]
length(all_units_bs$OTU) #54/3 = 18 OTUs 

#write.csv(all_units_bs, "~/Desktop/LAVO22_GG2/corncob/threeunits_burnseverity.csv")

#found in at least 2 units
two_units_bs <- temp_merge[unlist(tapply(1:nrow(temp_merge), temp_merge$OTU, function(x) if(length(x)>1) x)), ] 
length(two_units_bs$OTU) #268-54=214/2=107 OTUs diff abundant across only 2 units

write.csv(two_units_bs, "~/Desktop/LAVO22_GG2/corncob/twounits_burnseverity.csv")

two_units_trt$Family <- factor(two_units_trt$Family, levels=y_ord)
yord <- unique(two_units_bs$Genus)
two_units_bs$Genus <- factor(two_units_bs$Genus, levels=yord)

atleast2 <- ggplot(two_units_bs, aes(x=L2FC, y=Genus, color=Phylum, shape=comp)) + 
  theme_bw() +
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) +
  labs(x="Degree of Enrichment with Burn Severity (mean relative abundance log2-fold change)", y="Genus", size='Avg Relative Abundance', shape='Burn Severity') +
  facet_wrap(~unit) #+
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
   #+
theme(legend.position = 'none')

atleast2
#### 2) 





###### 2) plot OTUs that are enriched/depleted ######
#### a) fire history ####
color.ord <- unique(BL_trt$Phylum)
length(color.ord)
BL_trt$Phylum <- factor(BL_trt$Phylum, levels = color.ord)

### Butte Lake
a <- ggplot(BL_trt, aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Family, Genus, OTU", x="log2 fold change", size='More Fire Avg Relative Abundance') + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  xlim(c(-10,10))+ scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08))

a

z <- ggplot(BL_trt[1:216,], aes(x=L2FC, y=reorder(Genus, desc(Genus)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Genus", title='Butte Lake', size='More Fire Avg Relative Abundance', x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) +
  theme(legend.position = 'none') + geom_vline(xintercept = 0, aes(alpha=0.7))

z

z2 <- ggplot(BL_trt[217:432,], aes(x=L2FC, y=reorder(Genus, desc(Genus)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Genus", size='More Fire Avg Relative Abundance', x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) +
  theme(legend.position = "none") + geom_vline(xintercept = 0, aes(alpha=0.7))
z
z+z2

y <- ggplot(BL_trt, aes(x=L2FC, y=reorder(Family, desc(Family)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Family", size='More Fire Avg Relative Abundance', x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + 
  theme(legend.position = 'none') + labs(title='Butte Lake') + geom_vline(xintercept = 0, aes(alpha=0.7))

y

length(unique(BL_trt$Genus))

432/2

b <- ggplot(BL_trt[1:100,], aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point(stat='identity') + labs(y="Family, Genus, OTU", x="log2 fold change") + 
  scale_color_manual(values=c('#e6194B', '#f58231')) +xlim(c(-10,10)) +
  theme(legend.position = "none")
c <- ggplot(BL_trt[101:200,], aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point(stat='identity') + labs(y="Family, Genus, OTU", x="log2 fold change") +
  scale_color_manual(values=c('#f58231', '#ffe119', '#bfef45', '#3cb44b','#42d4f4')) + 
  xlim(c(-10,10)) + theme(legend.position = "none")
d <- ggplot(BL_trt[201:300,], aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change") +
  scale_colour_manual(values=c('#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000')) +
  xlim(c(-10,10)) + theme(legend.position = "none")
e <-ggplot(BL_trt[301:432,], aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change") +
  scale_colour_manual(values=c('#000000', '#fabed4','#ffd8b1')) + xlim(c(-10,10)) + theme(legend.position = "none")

b
c
d
e


### Hole
f <- ggplot(H_trt, aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point() + labs(y="Family, Genus, OTU", x="log2 fold change") + xlim(c(-10,10)) + 
  scale_color_manual(values=c('#e6194B', '#f58231', '#bfef45','#f032e6','#800000','#808000','#000075','#000000')) 

f

z.H <- ggplot(H_trt, aes(x=L2FC, y=reorder(Genus, desc(Genus)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Genus", x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)', title='Hole') + 
  scale_color_manual(values=c('#e6194B', '#f58231', '#bfef45','#f032e6','#800000','#808000','#000075','#000000')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + theme(legend.position = 'none') + 
  geom_vline(xintercept = 0, aes(alpha=0.7))
z.H

y.H <- ggplot(H_trt, aes(x=L2FC, y=reorder(Family, desc(Family)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Family", x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#e6194B', '#f58231', '#bfef45','#f032e6','#800000','#808000','#000075','#000000')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + 
  labs(size='More Fire Avg Relative Abundance', title='Hole') + theme(legend.position = 'none') + 
  geom_vline(xintercept = 0, aes(alpha=0.7))
y.H

### Warner Valley
g <- ggplot(WV_trt, aes(x=L2FC, y=reorder(Family_Genus_OTU, desc(Family_Genus_OTU)), color=Phylum)) + 
  geom_point() + labs(y="Family, Genus, OTU", x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + xlim(c(-10,10)) + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#42d4f4','#4363d8','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4'))
g

z.WV <- ggplot(WV_trt, aes(x=L2FC, y=reorder(Genus, desc(Genus)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Genus", x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)', title='Warner Valley') + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#42d4f4','#4363d8','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + theme(legend.position = 'none') + 
  geom_vline(xintercept = 0, aes(alpha=0.7))
z.WV

y.WV <- ggplot(WV_trt, aes(x=L2FC, y=reorder(Family, desc(Family)), color=Phylum)) + 
  geom_point(aes(size=Rel.Abund.MF)) + labs(y="Family", x='Enrichment in More Fire Plots (mean relative abundance log2-fold change)') + 
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#42d4f4','#4363d8','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4')) +
  xlim(c(-10,10)) + scale_size(breaks = c(0.0001, 0.001, 0.01, 0.05, 0.08)) + 
  labs(size='More Fire Avg Relative Abundance', title='Warner Valley') + theme(legend.position = 'none') + 
  geom_vline(xintercept = 0, aes(alpha=0.7))
  
y.WV

###units side by side
#Family
y + y.H + y.WV

#Genus
z +z2 
z.H+z.WV

#### b) burn severity ####

#all together but faceted by unit
df_bs_long$comp <- factor(df_bs_long$comp, levels=c('unburned', 'low', 'moderate', 'severe'))

x <- ggplot(df_bs_long, aes(x=comp, y=L2FC, color=Phylum, group=OTU)) + 
  theme_bw() +
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  geom_line() + scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) +
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  labs(y="Degree of Enrichment with Burn Severity (mean relative abundance log2-fold change)", x="Burn Severity", size='Avg Relative Abundance') +
  facet_wrap(~unit) #+
theme(legend.position = 'none')

x

#Butte Lake
BL_bs_long$OTU <- factor(BL_bs_long$OTU)
BL_bs_long$comp <- factor(BL_bs_long$comp, levels=c('unburned', 'low','moderate','severe'))

length(unique(BL_bs_long$Phylum)) #15 Phyla; missing Bdellovibrionata_E (dark green), Cyanobacteria (dark blue); Dormibacterota (purple)

#to determine the scale for relative abundance
fivenum(BL_bs_long$Rel.Abund) #min: 1.59e-6; Q1: 8.04e-05; median: 4.00e-04; Q3: 1.33e-03; max: 0.25

h <- ggplot(BL_bs_long, aes(x=comp, y=L2FC, color=`p<0.05`, group=OTU)) + geom_point() + geom_line() + labs(y="log2 fold change", x="burn severity", title="Butte Lake OTUs differentially abundant with burn severity") 
h

x.BL <- ggplot(BL_bs_long, aes(x=comp, y=L2FC, color=Phylum, group=OTU)) + 
  theme_bw() + 
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  geom_line() + scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) +
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#42d4f4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  labs(y="Log2 Fold Change", x="Burn Severity", title="Butte Lake", size='Avg Relative Abundance') +
  theme(legend.position = 'none')

x.BL

#Hole
H_bs_long$OTU <- factor(H_bs_long$OTU)
H_bs_long$comp <- factor(H_bs_long$comp, levels=c('unburned', 'low','moderate'))
i <- ggplot(H_bs_long, aes(x=comp, y=L2FC, color=`p<0.05`, group=OTU)) + geom_point() + geom_line() + labs(y="log2 fold change", x="burn severity", title="Hole OTUs differentially abundant with burn severity") 
i

x.H <- ggplot(H_bs_long, aes(x=comp, y=L2FC, color=Phylum, group=OTU)) + 
  theme_bw() +
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  geom_line() + scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) +
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#42d4f4','#911eb4','#f032e6','#800000','#808000','#469990','#000075','#000000', '#ffd8b1')) +
  labs(y="Log2 Fold Change", x="Burn Severity", title="Hole", size='Avg Relative Abundance') #+
  theme(legend.position = 'none')

x.H

#Warner Valley
WV_bs_long$OTU <- factor(WV_bs_long$OTU)
WV_bs_long$comp <- factor(WV_bs_long$comp, levels=c('low','moderate','severe'))
j <- ggplot(WV_bs_long, aes(x=comp, y=L2FC, color=`p<0.05`, group=OTU)) + geom_point() + geom_line() + labs(y="log2 fold change", x="burn severity", title="Warner Valley OTUs differentially abundant with burn severity") 
j

x.WV <- ggplot(WV_bs_long, aes(x=comp, y=L2FC, color=Phylum, group=OTU)) + 
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  geom_line() + scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) +
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#800000','#9A6324','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  labs(y="Log2 Fold Change", x="Burn Severity", title="Warner Valley", size='Avg Relative Abundance') +
  theme(legend.position = 'none')
x.WV

#### 3) plot OTUs that are enriched/depleted across units ####

#### 4) Massilia, Noviherbaspirilum, and Modestobacter 
MassModNov_bs <- read.csv("~/Desktop/LAVO22_GG2/corncob/MassModNov_bs.csv")
MassModNov_firehist <- read.csv("~/Desktop/LAVO22_GG2/corncob/MassModNov_firehist.csv")

#wide to long
#Butte Lake  
MassModNov_bs_clean <- select(MassModNov_bs, OTU:L2FC.low, L2FC.mod, L2FC.severe, Rel.Abund.unburned:Rel.Abund.low, Rel.Abund.moderate, Rel.Abund.severe, t.test_p_low:t.test_p_severe, Domain:Family, Genus, Species) %>%
  mutate(., X=row.names(MassModNov_bs_clean))

temp1 <- MassModNov_bs_clean %>%
  mutate(unburned=0) %>%
  select(X, unburned, L2FC.low:L2FC.severe) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="L2FC") %>%
  mutate(comp=ifelse(comparison=="unburned", "unburned", ifelse(comparison=="L2FC.low", "low", ifelse(comparison=="L2FC.mod", "moderate", "severe")))) %>%
  select(1, 3:4) 

temp2 <- MassModNov_bs_clean %>%
  select(X, t.test_p_low:t.test_p_severe) %>%
  pivot_longer(cols=!X, names_to="comparison", values_to="t.test_p") %>%
  mutate(comp=ifelse(comparison=="t.test_p_low", "low", ifelse(comparison=="t.test_p_moderate", "moderate", "severe")))%>%
  select(1, 3:4) 

temp3 <- MassModNov_bs_clean %>%
  select(X, Rel.Abund.unburned:Rel.Abund.severe) %>%
  pivot_longer(cols=!X, names_to='comparison', values_to='Rel.Abund') %>%
  mutate(comp=ifelse(comparison=="Rel.Abund.unburned", 'unburned', ifelse(comparison=='Rel.Abund.low', 'low', ifelse(comparison=='Rel.Abund.moderate', 'moderate', 'severe')))) %>%
  select(1, 3:4)

temp_merge <- temp1 %>%
  left_join(temp2, by = c('X', 'comp')) %>%
  left_join(temp3, by = c('X', 'comp'))

taxonomy <- MassModNov_bs_clean %>%
  select(X, OTU, unit, Domain:Species)

MassModNov_bs_long <- left_join(temp_merge, taxonomy, by='X') %>%
  mutate(., Genus = ifelse(is.na(Genus), "unclassified", ifelse(Genus=="", "unclassified", Genus))) %>%
  mutate(., Family=ifelse(Family=="", "unclassified", Family)) %>%
  mutate(., Family_Genus_OTU = paste(Family, Genus, OTU, sep=", ")) %>%
  mutate(., `p<0.05` = ifelse(t.test_p<0.05, "TRUE", "FALSE")) %>%
  mutate(., a = ifelse(t.test_p<0.05, 0.8, 0.3)) %>%
  ungroup(.) %>%
  arrange(Phylum, L2FC)
y_ord <- unique(MassModNov_bs_long$Family_Genus_OTU)
MassModNov_bs_long$Family_Genus_OTU <- factor(MassModNov_bs_long$Family_Genus_OTU, levels=y_ord)

rm(temp1, temp2, temp_merge, taxonomy)


ggplot(MassModNov_bs_long, aes(x=L2FC, y=OTU)) + geom_point()

x <- ggplot(MassModNov_bs_long, aes(x=comp, y=L2FC, color=Phylum, group=OTU)) + 
  theme_bw() +
  geom_point(aes(alpha=a, size=Rel.Abund)) + scale_alpha_continuous(guide=FALSE) +
  geom_line() + scale_size(breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.1)) #+
  scale_color_manual(values=c('#e6194B', '#f58231','#ffe119', '#bfef45', '#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#a9a9a9','#800000','#9A6324','#808000','#469990','#000075','#000000','#fabed4','#ffd8b1')) +
  labs(y="Degree of Enrichment with Burn Severity (mean relative abundance log2-fold change)", x="Burn Severity", size='Avg Relative Abundance') +
  facet_wrap(~unit) #+
theme(legend.position = 'none')

x

