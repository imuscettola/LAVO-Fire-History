#Objective: To create an algorithm that randomly picks samples from a plot in 
#order to compare burn severity models vs. community dissimilarity. This algorithm
#will test objective 2.1: assessing which burn severity metric is more correlated
#with bacterial community dissimilarity, essentially asking the question of which
#spatial scale of burn severity is most informative of changes to bacterial communities
#induced by fire. 
####of course this assumes that the assumption here is that the communities started
#at roughly the same composition and only fire is explaining community composition
#unless I add covariates such as soil type, vegetation type, or soil pH. 

#################Loading libraries and phyloseq object########################

library(dplyr)
library(tidyr)
library(ggplot2)
library(phyloseq)
library(vegan)
library(viridisLite)

#load phyloseq object
#minimac
#ps.3 <- readRDS("~/Desktop/LAVO22_GG2/ps.3units")
ps.3 <- readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")

#remove samples (total of 4 samples (2 different quadrats)) that are missing 
#burn severity metrics; including these samples would result in an error in the 
#PERMANOVA for having "NA" values
ps.3 <- subset_samples(ps.3, BARE_ROCK!='NA' & Wt_Avg_QuadBS!='NA')

#hellinger transformation for PERMANOVAs
hell.ps.3 = transform_sample_counts(ps.3, function(x) (x / sum(x))^0.5 )

rm(ps.3)
sample_data(hell.ps.3)

sample_data(hell.ps.3)$mtbs_sev

#remove BL ABCO 04 bc mtbs value = 5
#2/23 not necessary because mtbs data downloaded directly from their site (not 
#from Google Earth) plots overlaid did not have values > 4.
#hell.ps.3 <- subset_samples(hell.ps.3, Plot_Name!='BL_ABCO_04')


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
cohort_df$O_Samples[!(cohort_df$O_Samples %in% as.data.frame(sample_data(hell.ps.3))$Full_Sample_Name)] = "NA" #replaces O samples that don't exist in the dataset with "NA"; works when n=3!
cohort_df$M_Samples <- paste(cohort_df$Sample_Core_Name, "M", sep="")


###############Subset phyloseq object and run permanovas######################
#empty dataframe to store coefficients
final_df = data.frame(cohort=numeric(), horizon=character(), bs_metric=character(), total_df=numeric(),R2=numeric(), F=numeric(), p_fdr=numeric())


for(i in 1:max(cohort_df$Cohort)){
  

#list of Full Name of Mineral samples for cohort 1
M_Samples <- cohort_df[cohort_df$Cohort==i, ]$M_Samples
O_Samples <- cohort_df[cohort_df$Cohort==i, ]$O_Samples


#create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
ps.M <- subset_samples(hell.ps.3, Full_Sample_Name %in% M_Samples) 
ps.O <- subset_samples(hell.ps.3, Full_Sample_Name %in% O_Samples) 

# Get sample_data as dataframe for each phyloseq object
SamDat.M = data.frame(sample_data(ps.M))
SamDat.O = data.frame(sample_data(ps.O))

#make dissimilarity matrix
dist.M = phyloseq::distance(ps.M, method="bray")
dist.O = phyloseq::distance(ps.O, method="bray")

###PERMANOVAs
#substrate
substr.adonis.M = adonis2(dist.M ~  BARE_ROCK, strata = SamDat.M$Unit, SamDat.M)
substr.adonis.M$row = row.names(substr.adonis.M)
substr.adonis.O = adonis2(dist.O ~  BARE_ROCK, strata = SamDat.O$Unit, SamDat.O)
substr.adonis.O$row = row.names(substr.adonis.O)

#samplesev
samplesev.adonis.M = adonis2(dist.M ~  SampleSev, strata = SamDat.M$Unit, SamDat.M)
samplesev.adonis.M$row = row.names(samplesev.adonis.M)
samplesev.adonis.O = adonis2(dist.O ~  SampleSev, strata = SamDat.O$Unit, SamDat.O)
samplesev.adonis.O$row = row.names(samplesev.adonis.O)

#wt avg of quadrat
avgquadsev.adonis.M = adonis2(dist.M ~  Wt_Avg_QuadBS, strata = SamDat.M$Unit, SamDat.M)
avgquadsev.adonis.M$row = row.names(avgquadsev.adonis.M)
avgquadsev.adonis.O = adonis2(dist.O ~  Wt_Avg_QuadBS, strata = SamDat.O$Unit, SamDat.O)
avgquadsev.adonis.O$row = row.names(avgquadsev.adonis.O)

#Canopy mortality
canopymort.adonis.M = adonis2(dist.M ~  Overstory_mort, strata = SamDat.M$Unit, SamDat.M)
canopymort.adonis.M$row = row.names(canopymort.adonis.M)
canopymort.adonis.O = adonis2(dist.O ~  Overstory_mort, strata = SamDat.O$Unit, SamDat.O)
canopymort.adonis.O$row = row.names(canopymort.adonis.O)

#MTBS burn severity
mtbs.adonis.M = adonis2(dist.M ~  mtbs_sev, strata = SamDat.M$Unit, SamDat.M)
mtbs.adonis.M$row = row.names(mtbs.adonis.M)
mtbs.adonis.O = adonis2(dist.O ~  mtbs_sev, strata = SamDat.O$Unit, SamDat.O)
mtbs.adonis.O$row = row.names(mtbs.adonis.O)

#########Now that we have the permanovas made, let's store the coefficients

#for loop to create permanovas and store coefficients in dataframe

#for loop to make PERMANOVAs for a single cohort 
#(1:4 because 4 burn severity metrics to compare to start)

#make dataframe containing all of the outputs of the PERMANOVAs
#can't start for loop from 0 so need to make an empty dataframe size 5x3 to 
#include at the beginning 
temp = data.frame(replicate(6, c(0,0,0)))
names(temp) <- names(substr.adonis.M)

b = rbind(temp, substr.adonis.M, samplesev.adonis.M, avgquadsev.adonis.M, canopymort.adonis.M, mtbs.adonis.M, substr.adonis.O, samplesev.adonis.O, avgquadsev.adonis.O, canopymort.adonis.O, mtbs.adonis.O)

a = data.frame(cohort=numeric(), horizon=character(), bs_metric=character(), total_df=numeric(),R2=numeric(), F=numeric(), p_fdr=numeric())

for(x in 1:10) {
 a[x, 1] = i #cohort number
 a[x, 2] = ifelse(x<6, 'M', 'O')
 a[x, 3] = b$row[x*3+1]
 a[x, 4] = b$Df[x*3+2]
 a[x, 5] = b$R2[x*3+1]
 a[x, 6] = b$F[x*3+1]
 a[x, 7] = b$`Pr(>F)`[x*3+1]
 
}
final_df <- rbind(final_df, a, stringsAsFactors=FALSE)

}

##################SAVE OUTPUT #############################################
#first save
df <- mutate(final_df, total=n)

#save as csv file
#write.csv(df, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_500raw_nomtbs5.csv")

###################SUMMARIZE OUTPUT########################################

summary <- final_df %>%
  .[order(.$cohort, .$horizon, .$R2, decreasing = TRUE), ] %>%
  group_by(cohort, horizon) %>%
  mutate(R2_order = paste0(bs_metric, collapse=">")) %>%
  group_by(horizon, R2_order) %>%
  distinct(cohort, .keep_all = TRUE) %>%
  summarise(x=n()) %>%
  mutate(total=n, perc=x/total)

write.csv(summary, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_results.csv")

#summary_all <- summary

summary_all <- rbind(summary_all, summary)
#save as csv
#write.csv(summary_all, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_results.csv")



################### VISUALIZE RESULTS ##########################################
###plot results of iterative tests
#summary_all <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_results.csv")
summary_all <- read.csv("~/Desktop/LAVO22_GG2/2.1_rankingseveritymetrics_results.csv")
summary_all$R2_order <-as.factor(summary_all$R2_order)
ggplot(summary_all, aes(x=total, y=perc)) + geom_point(size=3, aes(color=R2_order)) + facet_wrap(~horizon) +theme_bw()

##plot results of ranks for 500 no mtbs5 dataset
#ggplot(raw500_mtbscompare[raw500_mtbscompare$mtbs5=="N", ], aes(x=rank, y=rank.perc)) + geom_point(size=3, aes(color=bs_metric)) + facet_wrap(~horizon) +theme_bw() + labs(title="Ranks of R2 values comparing burn severity metrics (no mtbs=5)")


###what is the effect size (difference in R2 for top 2 burn severity metrics?)
final_df <- read.csv("~/Desktop/LAVO22_GG2/2.1_rankingseveritymetrics_500raw_nomtbs5.csv")
rank_500 <- final_df %>%
  arrange(cohort, horizon, R2) %>%
  group_by(cohort, horizon) %>%
  mutate(rank=rank(desc(R2))) %>%
  mutate(p.range=ifelse(p_fdr>0.05, 'p>0.05', ifelse(p_fdr==0.001, 'p<0.001', 'p<0.05')))
rank_500$rank = factor(rank_500$rank)
rank_500$bs_metric=factor(rank_500$bs_metric, levels = c("SampleSev", "BARE_ROCK", "Wt_Avg_QuadBS", "Overstory_mort", "mtbs_sev"))
rank_500$hor <- ifelse(rank_500$horizon=="M", "Mineral", "Organic")

####################################FIGURE FOR THESIS##########################
#plot ranks
ggplot(rank_500, aes(x=rank, y=R2, fill=bs_metric)) + 
  geom_boxplot() + facet_wrap(~hor) + theme_bw() + 
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12,)) +
  scale_fill_discrete(name="Burn Severity Metric", labels=c("Core Substrate Burn Severity", "% Exposed Mineral Soil", "Quadrat Average Substrate Burn Severity", "Canopy Mortality", "MTBS thresholded dNBR")) #+
  theme(legend.position = "none")
  

#summarise mean R2 for each rank and burn severity metric
R2_sum <- rank_500 %>%
  ungroup() %>%
  group_by(horizon, bs_metric, rank, p.range) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))
#save output
#write.csv(R2_sum, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_rankR2summary.csv")

#summarise p-values for each burn severity metric
p_sum <- rank_500 %>%
  ungroup() %>%
  group_by(horizon, bs_metric, p.range) %>%
  summarise(n=n(), perc=n/500)

#difference in R2 for the top 2 ranked burn severity metrics (rank1 R2 - rank2 R2)
#df_sum <- filter(mtbs5N_rank, rank=="1"|rank=="2") %>%
#  mutate(., diff=R2-lag(R2))

#group by rank order of first two variables
df_sum <- df_sum %>%
  .[order(.$cohort, .$horizon, .$rank, decreasing = FALSE), ] %>%
  group_by(cohort, horizon) %>%
  mutate(R2_order = paste0(bs_metric, collapse=">"))

#summarize by groups of rank 1 & 2 bs metrics
diff_sum <- df_sum %>%
  filter(rank==1) %>%
  ungroup() %>%
  group_by(horizon, R2_order) %>%
  summarise(n=n(), meandiff=mean(diff), sd=sd(diff))

temp <- diff_sum %>%
  filter(rank==1) %>%
  ungroup() %>%
  summarise(n=n(), meandiff=mean(diff), sd=sd(diff))

#save output 
write.csv(df_sum, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_R2diffsummary_nomtbs5.csv")

###what is avg R2 for rank 1? 
#subset by bs_metric and horizon
R2_sum <- df_sum %>%
  filter(rank==1) %>%
  ungroup() %>%
  group_by(horizon, bs_metric) %>%
  summarise(n=n(), meanR2=mean(R2), sdR2=sd(R2))

#all together subset by horizon
temp <- df_sum %>%
  filter(rank==1) %>%
  ungroup() %>%
  group_by(horizon) %>%
  summarise(n=n(), mean.R2=mean(R2), sd.R2=sd(R2))


###########are these trends different than the results from a full model?######


#create 2 subsetted phyloseq objects for full data; one for mineral one for organic
ps.M <- subset_samples(hell.ps.3, Horizon=="M") 
ps.O <- subset_samples(hell.ps.3, Horizon=="O") 

# Get sample_data as dataframe for each phyloseq object
SamDat.M = data.frame(sample_data(ps.M))
SamDat.O = data.frame(sample_data(ps.O))

#make dissimilarity matrix
dist.M = phyloseq::distance(ps.M, method="bray")
dist.O = phyloseq::distance(ps.O, method="bray")

###PERMANOVAs
#substrate
substr.adonis.M = adonis2(dist.M ~  BARE_ROCK, strata = SamDat.M$Unit, SamDat.M)
substr.adonis.M$row = row.names(substr.adonis.M)
substr.adonis.O = adonis2(dist.O ~  BARE_ROCK, strata = SamDat.O$Unit, SamDat.O)
substr.adonis.O$row = row.names(substr.adonis.O)

#samplesev
samplesev.adonis.M = adonis2(dist.M ~  SampleSev, strata = SamDat.M$Unit, SamDat.M)
samplesev.adonis.M$row = row.names(samplesev.adonis.M)
samplesev.adonis.O = adonis2(dist.O ~  SampleSev, strata = SamDat.O$Unit, SamDat.O)
samplesev.adonis.O$row = row.names(samplesev.adonis.O)

#wt avg of quadrat
avgquadsev.adonis.M = adonis2(dist.M ~  Wt_Avg_QuadBS, strata = SamDat.M$Unit, SamDat.M)
avgquadsev.adonis.M$row = row.names(avgquadsev.adonis.M)
avgquadsev.adonis.O = adonis2(dist.O ~  Wt_Avg_QuadBS, strata = SamDat.O$Unit, SamDat.O)
avgquadsev.adonis.O$row = row.names(avgquadsev.adonis.O)

#Canopy mortality
canopymort.adonis.M = adonis2(dist.M ~  Overstory_mort, strata = SamDat.M$Unit, SamDat.M)
canopymort.adonis.M$row = row.names(canopymort.adonis.M)
canopymort.adonis.O = adonis2(dist.O ~  Overstory_mort, strata = SamDat.O$Unit, SamDat.O)
canopymort.adonis.O$row = row.names(canopymort.adonis.O)

#MTBS burn severity
mtbs.adonis.M = adonis2(dist.M ~  mtbs_sev, strata = SamDat.M$Unit, SamDat.M)
mtbs.adonis.M$row = row.names(mtbs.adonis.M)
mtbs.adonis.O = adonis2(dist.O ~  mtbs_sev, strata = SamDat.O$Unit, SamDat.O)
mtbs.adonis.O$row = row.names(mtbs.adonis.O)

###add PERMANOVA results together and pull coefficients
temp = data.frame(replicate(6, c(0,0,0)))
names(temp) <- names(substr.adonis.M)

b = rbind(temp, substr.adonis.M, samplesev.adonis.M, avgquadsev.adonis.M, canopymort.adonis.M, mtbs.adonis.M, substr.adonis.O, samplesev.adonis.O, avgquadsev.adonis.O, canopymort.adonis.O, mtbs.adonis.O)

a.concat = data.frame(cohort=character(), horizon=character(), bs_metric=character(), total_df=numeric(),R2=numeric(), F=numeric(), p_fdr=numeric())

for(x in 1:10) {
  a = data.frame(cohort=character(), horizon=character(), bs_metric=character(), total_df=numeric(),R2=numeric(), F=numeric(), p_fdr=numeric())
  a[x, 1] = 'full' #cohort number
  a[x, 2] = ifelse(x<6, 'M', 'O')
  a[x, 3] = b$row[x*3+1]
  a[x, 4] = b$Df[x*3+2]
  a[x, 5] = b$R2[x*3+1]
  a[x, 6] = b$F[x*3+1]
  a[x, 7] = b$`Pr(>F)`[x*3+1]
  a.concat <- rbind(a.concat, a)
}

summary.full <- a.concat %>%
  .[order(.$cohort, .$horizon, .$R2, decreasing = TRUE), ] %>%
  group_by(cohort, horizon) %>%
  mutate(R2_order = paste0(bs_metric, collapse=">")) %>%
  group_by(horizon, R2_order) %>%
  distinct(cohort, .keep_all = TRUE) %>%
  summarise(x=n()) %>%
  filter(horizon!="NA")
#M horizon: Wt_Avg_QuadBS>SampleSev>Overstory_mort>mtbs_sev>BARE_ROCK
#O horizon: SampleSev>Wt_Avg_QuadBS>Overstory_mort>mtbs_sev>BARE_ROCK

#R2 differences of top 2 R2 values
#M
0.10173775-0.10069109 #0.00104666
0.17813317-0.1743642 #0.00376897


###################OLD ##################################################
#############COMPARING RESULTS W/ & W/O MTBS=5 ################################
#compare removing the plot with mtbs = 5 -> does this change the results? 
#mtbs5Y <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_500raw_mtbs5incl.csv")
#mtbs5N <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/2.1_rankingseveritymetrics_500raw_nomtbs5.csv")

summary_all <- mutate(summary_all, mtbs5="Y")
summary <- mutate(summary, mtbs5="N")
summary_5compare <- rbind(summary_all, summary) %>%
  filter(., total==500)

raw500_mtbscompare <- rbind(mutate(mtbs5N, mtbs5="N"), mutate(mtbs5Y, mtbs5="Y")) %>%
  arrange(mtbs5, cohort, horizon, R2) %>%
  group_by(mtbs5, cohort, horizon) %>%
  mutate(rank=rank(desc(R2))) %>%
  group_by(mtbs5, horizon, rank, bs_metric) %>%
  summarise(rank.n = n()) %>%
  mutate(rank.perc=rank.n/500)

ggplot(raw500_mtbscompare, aes(x=rank, y=rank.perc)) + geom_point(size=3, aes(color=bs_metric, shape=mtbs5)) + facet_wrap(~horizon) +theme_bw() + labs(title="Ranks of R2 values comparing burn severity metrics (no mtbs=5)")
