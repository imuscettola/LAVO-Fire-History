############## Are microbial communities independent within a plot? ###########
#2/16/2024
#Objective: to build models describing community composition and compare the 
#results when 1) samples are treated as independent, 2) 1 sample per plot is 
#randomly selected, and 3) samples are aggregated within a plot by average a) 
#relative abundance and b) otu counts. Are the treatment (more v less fire) 
#effects still significant treating the samples this way? 
#note: because the models will use differing numbers of samples, direct comparisons
#of R2 values are not possible. Instead we will test whether the findings from 
#the community models are qualitatively the same by 
#2/21/2024 
#update: veg type was adjusted for 1 BL plot and 6 H plots to reflect dominant 
#vegetation immediately prior to Dixie fire (rather than the veg prior to Reading
#for Hole plots and prior to last Rx burn in BL unit)
#reclassified veg type changed adonis output: 
##independent p-value from 0.001 to 0.003 
##
install.packages("dplyr")
install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install("phyloseq")
install.packages("vegan")
install.packages("ggplot2")

#load libraries and data
library(dplyr)
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)

################### 1 Samples treated as independent #######################
###Is there a significant effect of fire history on microbial community similarity?
#load phyloseq object
ps.3 <- readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")
#hellinger transformation
hell.ps.3 <- transform_sample_counts(ps.3, function(x) (x / sum(x))^0.5 )
sample_data(hell.ps.3)$Soil_mu <- as.factor(sample_data(hell.ps.3)$Soil_mu)

#distance matrix
dist.hell <- distance(hell.ps.3, method="bray")
#sample data
SamDat <- data.frame(sample_data(hell.ps.3))

#test if Trt is a significant predictor of microbial community similarity when 
#controlling for horizon, soil type, and vegetation type
hell.adonis <- adonis2(dist.hell ~ Horizon + Soil_mu + Veg_Type + Trt, strata = SamDat$Unit, SamDat)
hell.adonis 

#save output
independent <- hell.adonis
independent$test <- "independent"
independent$row <- rownames(independent)

write.csv(independent, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/independentsamples_full.csv")

################## 2 One sample per plot ####################################
#####################Create dataframe of cohorts################################
#get number of plots 
temp <- data.frame(sample_data(hell.ps.3)) %>%
  distinct(Plot_Name, .keep_all=TRUE)
plots <- temp$Plot_Name
length(plots) #58 plots

rm(temp)

#how many cohorts (group of randomly selected plots) do you want to test
n <- 500

#create dataframe of the samples randomly selected for each cohort  
cohort_df <- data.frame(Plot_Name=character(), Rep=character(), Cohort=numeric())

for(i in 1:n) {
  set.seed(i)
  my_list = list(Plot_Name=plots, Rep = sample(c("A", "B", "C"), length(plots), replace = TRUE), Cohort=replicate(n=length(plots), i))
  cohort_df = rbind(cohort_df, my_list, stringsAsFactors=FALSE)
}
rm(i)

#create columns to subset out phyloseq object
cohort_df$Sample_Core_Name <- paste(cohort_df$Plot_Name, cohort_df$Rep, sep="_")
cohort_df$O_Samples <- paste(cohort_df$Sample_Core_Name, "O", sep="") 
cohort_df$O_Samples[!(cohort_df$O_Samples %in% as.data.frame(sample_data(hell.ps.3))$Full_Sample_Name)] = "NA" #replaces O samples that don't exist in the dataset with "NA"
cohort_df$M_Samples <- paste(cohort_df$Sample_Core_Name, "M", sep="")
cohort_df$M_Samples[!(cohort_df$M_Samples %in% as.data.frame(sample_data(hell.ps.3))$Full_Sample_Name)] = "NA" # replaces M samples that don't exist in the dataset with "NA"

#column with all samples (organic and mineral) for a cohort
cohort_df <- pivot_longer(cohort_df, cols="O_Samples":"M_Samples", names_to="Horizon", values_to="Full_Sample_Name") %>%
  filter(., Full_Sample_Name!="NA")

###############Subset phyloseq object and run permanovas######################
#empty dataframe to store coefficients
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), `Pr(>F)`=numeric(), test=character(), row=character(), cohort=numeric())


for(i in 1:max(cohort_df$Cohort)){
  
  ###Get ps object, sample data, and distance matrix for the cohort
  #list of Full Name of Mineral samples for cohort 1
  Samples <- cohort_df[cohort_df$Cohort==i, ]$Full_Sample_Name
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps <- subset_samples(hell.ps.3, Full_Sample_Name %in% Samples) 
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps))
  
  #make dissimilarity matrix
  dist = phyloseq::distance(ps, method="bray")
  
  ###PERMANOVAs
  #
  i.adonis = adonis2(dist ~ Horizon + Soil_mu + Veg_Type + Trt, strata = SamDat$Unit, SamDat)
  i.adonis$test='one.per.plot'
  i.adonis$row=rownames(i.adonis)
  i.adonis$cohort=i
  
  final_df <- rbind(final_df, i.adonis)
  
}

##################SAVE OUTPUT #############################################

random_reps_500 <- mutate(final_df, p.range=ifelse(`Pr(>F)`>0.05, 'p>0.05', ifelse(`Pr(>F)`==0.001, 'p<0.001', 'p<0.05')))

random_reps_500_summary <- random_reps_500 %>%
  group_by(row, p.range) %>%
  summarise(n=n(), mean.R2=mean(R2), sd.R2=sd(R2)) %>%
  mutate(n.perc = n/500) %>%
  filter(p.range!='NA')
#93.6% of cohorts had Trt p-value < 0.001; 6.4% of cohorts had Trt p-value < 0.05\

#save output
write.csv(random_reps_500_summary, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/independentsamples_cohorts.csv")



################## 3 Plot averages ###########################################
### a) averaging relative abundance within a plot

#ps.3.norm.avg <-readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3.norm.avg")
ps.3.norm.avg <- readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3.norm.avg")
sample_data(ps.3.norm.avg)
hell.ps.3.norm.avg <- transform_sample_counts(ps.3.norm.avg, function(x) x^0.5 )
sample_data(hell.ps.3.norm.avg)$Soil_mu <- as.factor(sample_data(hell.ps.3.norm.avg)$Soil_mu)

#distance matrix
dist.hell <- distance(hell.ps.3.norm.avg, method="bray")
#sample data
SamDat <- data.frame(sample_data(hell.ps.3.norm.avg))

#test if Trt is a significant predictor of microbial community similarity when 
#controlling for horizon, soil type, and vegetation type
hell.adonis.norm.avg <- adonis2(dist.hell ~ Horizon + Soil_mu + Veg_Type + Trt, strata = SamDat$Unit, SamDat)
hell.adonis.norm.avg 

#save output
avg.relabund <- hell.adonis.norm.avg
avg.relabund$test <- "avg relative abundance"
avg.relabund$row <- rownames(avg.relabund)

write.csv(avg.relabund, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/independentsamples_avgrelabund.csv")

### b) averaging OTU counts within a plot
ps.3.avg <- readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3.avg")
ps.3.avg <- readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3.avg")

sample_data(ps.3.avg)
hell.ps.3.avg <- transform_sample_counts(ps.3.avg, function(x) (x/sum(x))^0.5 )
sample_data(hell.ps.3.avg)$Soil_mu <- as.factor(sample_data(hell.ps.3.avg)$Soil_mu)

#distance matrix
dist.hell <- distance(hell.ps.3.avg, method="bray")
#sample data
SamDat <- data.frame(sample_data(hell.ps.3.avg))

#test if Trt is a significant predictor of microbial community similarity when 
#controlling for horizon, soil type, and vegetation type
hell.adonis.avg <- adonis2(dist.hell ~ Horizon + Soil_mu + Veg_Type + Trt, strata = SamDat$Unit, SamDat)
hell.adonis.avg 

#save output
avgOTU <- hell.adonis.avg
avgOTU$test <- "avg OTU"
avgOTU$row <- rownames(avgOTU)

write.csv(avgOTU, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/independentsamples_avgOTU.csv")




############# How are burn severity covariates spatially related? ###############
#2/8/2024
#Objective: To understand the the spatial autocorrelation of sample variables 
#including 1) microbial community dissimilarity, 2) sample burn severity, and 
#3) average quadrat burn severity. This will help to determine the independence
#of variables and whether replication within the plot should be treated as 
#true replicates or pseudo-replicates that need to be pooled or averaged for 
#subsequent analyses.

###load libraries and data
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#minimac
data.plot <- read.csv("~/Desktop/LAVO22_GG2/LAVO22_PlotMetadata.csv")
data.sample <- read.csv("~/Desktop/LAVO22_GG2/LAVO22_SampleMetadata.csv")
data.quad <- data.sample %>%
  distinct(Sample_Core_Name, .keep_all = TRUE)


###################### 1 Bacterial community dissimilarity ####################
### Spatial autocorrelation of microbial community dissimilarity: see 
#ordinations in Notion: 16S Data Analysis >> Ordinations "Are replicates within a 
#plot pseudo-replicates? 

#key finding: there does appear to be some clustering of microbial community by 
#plot although the clustering is not totally separate from other samples indicating
#that there is some spatial autocorrelation but generally microbial communities 
#are so diverse that there are large differences in commuity composition 

#caution: ordinations only show dissimilarity along a 2 dimensional space and are 
#not a statistical test. Also, the ordination takes the entire dataset together 
#rather than looking at the dataset subset by unit. 



##################### 2 Sample burn severity #################################
### Is the variation of measured sample burn severity within a plot equal to 
# the variation between plots? 
#Approach: ANOVA of burn severity by Plot Name subset by Unit and Treatment 
# & Boxplot vizualization 

### ANOVA
# null hypothesis: group means are all equal; Ha: at least one group mean is diff
# from the rest -> Tukey posthoc to determine which are different from each other

# BL MF 
aov.sample.bs.BL.MF <- aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="BL" & data.quad$Trt=="MF", ])
summary(aov.sample.bs.BL.MF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.BL.MF) #17 pairwise comparisons are diff out of 36

# BL LF
aov.sample.bs.BL.LF <-aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="BL" & data.quad$Trt=="LF", ])
summary(aov.sample.bs.BL.LF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.BL.LF) #29 pairwise comparisons are diff out of 36

# H MF
aov.sample.bs.H.MF <- aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="H" & data.quad$Trt=="MF", ])
summary(aov.sample.bs.H.MF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.H.MF) #10 pairwise comparisons are diff out of 21

# H LF
aov.sample.bs.H.LF <- aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="H" & data.quad$Trt=="LF", ])
summary(aov.sample.bs.H.LF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.H.LF) #7 pairwise comparisons are diff out of 21

#WV MF
aov.sample.bs.WV.MF <- aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="WV" & data.quad$Trt=="MF", ])
summary(aov.sample.bs.WV.MF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.WV.MF) #17 pairwise comparisons are diff out of 36

#WV LF
aov.sample.bs.WV.LF <- aov(SampleSev ~ Plot_Name, data=data.quad[data.quad$Unit=="WV" & data.quad$Trt=="LF", ])
summary(aov.sample.bs.WV.LF) #p<0.001 -> at least one group sample burn severity is diff
TukeyHSD(aov.sample.bs.WV.LF) #24 pairwise comparisons are diff out of 36

## Visualize 
#BL
a <- ggplot(data.quad[data.quad$Unit=="BL", ], aes(x=Plot_Name, y=SampleSev)) + geom_boxplot(aes(color=Trt)) + labs(title="Butte Lake Sample Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
a

#H
b <- ggplot(data.quad[data.quad$Unit=="H", ], aes(x=Plot_Name, y=SampleSev)) + geom_boxplot(aes(color=Trt)) + labs(title="Hole Sample Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
b

#WV
c <- ggplot(data.quad[data.quad$Unit=="WV", ], aes(x=Plot_Name, y=SampleSev)) + geom_boxplot(aes(color=Trt)) + labs(title="Warner Valley Sample Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
c

a+b+c

############# 3 Weighted Avg of Quadrat Burn Severity ##############
##same approach as above although skipping the anovas

##Visualize
#BL
a <- ggplot(data.quad[data.quad$Unit=="BL", ], aes(x=Plot_Name, y=Wt_Avg_QuadBS)) + geom_boxplot(aes(color=Trt)) + labs(title="Butte Lake Quadrat Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
a

#H
b <- ggplot(data.quad[data.quad$Unit=="H", ], aes(x=Plot_Name, y=Wt_Avg_QuadBS)) + geom_boxplot(aes(color=Trt)) + labs(title="Hole Quadrat Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
b

#WV
c <- ggplot(data.quad[data.quad$Unit=="WV", ], aes(x=Plot_Name, y=Wt_Avg_QuadBS)) + geom_boxplot(aes(color=Trt)) + labs(title="Warner Valley Quadrat Burn Severity by Plot") + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
c

a+b+c
