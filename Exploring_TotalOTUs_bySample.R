#Isabella Muscettola
#11/22/23 
#objective: to visualize the OTU relative abundances for each sample to understand
#if there were any systematic biases during the library prep process that can 
#help to explain trends in the data

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)


##############Creating dataframes; only do once################################
###########Create total OTUs dataframe
# How many sequences are in each sample. 
#import phyloseq object
ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
#transform to Hellinger relative abundances
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5)
#pull metadata from ps object
df = sample_data(ps)

# defining variable d as the column sums of our otu table.
d = sample_sums(ps)
# Creating a dataframe of our sample names
d = data.frame(names(d),d)
# Naming the columns
colnames(d)=c("Submission_Sample_Name","Total")
# Add other columns to allow to look at groups of data at the same time
d = left_join(d, df, by="Submission_Sample_Name")
d = d %>%
  select(Submission_Sample_Name:Total, Full_Sample_Name:Sequal_Rd)
#Save dataframe as csv
#write.csv(d, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Sample_Total_OTUs.csv")



###################Total OTUs per Sample #########################
#import dataframe
d = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Sample_Total_OTUs.csv")


#Graph Total OTUs per Sample Histogram
p = ggplot(d, aes(Total)) + theme_bw()
p = p + geom_histogram(binwidth=1000)
p = p + labs(x="Total OTUs", y="Sample Count", title = "All Samples Total OTUs Silva")
p

#Do the same for all the raw sequences (including blanks, positive controls, and
#duplicates)
ps.raw = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.raw")
d.raw = sample_sums(ps.raw)
df.raw = sample_data(ps.raw)
d.raw = data.frame(names(d.raw), d.raw)
colnames(d.raw) = c("Submission_Sample_Name", "Total")
d.raw = left_join(d.raw, df.raw, by = "Submission_Sample_Name")
d.raw = d.raw %>%
  select(Submission_Sample_Name:Total, Full_Sample_Name:Sequal_Rd)
#ps.raw Graph Total OTUs per Sample Histogram
p = ggplot(d.raw, aes(Total)) + theme_bw()
p = p + geom_histogram(binwidth=1000)
p = p + labs(x="Total OTUs", y="Sample Count", title = "Raw All Samples Total OTUs Silva")
p

#Do the same for only the 3 units (no blanks, pos controls, duplicates, or UB)
ps.3 = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3units")
d.3 = sample_sums(ps.3)
df.3 = sample_data(ps.3)
d.3 = data.frame(names(d.3), d.3)
colnames(d.3) = c("Submission_Sample_Name", "Total")
d.3 = left_join(d.3, df.3, by = "Submission_Sample_Name")
d.3 = d.3 %>%
  select(Submission_Sample_Name:Total, Full_Sample_Name:Sequal_Rd)
#ps.3 Graph Total OTUs per Sample Histogram
p = ggplot(d.3, aes(Total)) + theme_bw()
p = p + geom_histogram(binwidth=1000)
p = p + labs(x="Total OTUs", y="Sample Count", title = "3 Units All Samples Total OTUs Silva")
p



###################Block  Effects of Total OTUs##########################
# now I'm interested to see if there is a significant effect of block at each 
#step of the library prep in high (or low) OTU totals. High OTUs for a block of 
#samples may indicate there was contamination. Low OTUs could suggest poor following 
#of protocol or some other issue.

#Change Rounds to Characters instead of numeric & reorder in ascending
d$Wt_Rd = as.character(d$Wt_Rd)
d$Extract_Rd = as.character(d$Extract_Rd)
d$PCR_Rd = as.character(d$PCR_Rd)
d$Sequal_Rd = as.character(d$Sequal_Rd)
d$Soil_mu = as.character(d$Soil_mu)
d = d %>%
  mutate(Wt_Rd = factor(Wt_Rd, levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17", "18"))) %>%
  mutate(Extract_Rd = factor(Extract_Rd, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23","24","25","26","27","28","29","30","31"))) %>%
  mutate(PCR_Rd = factor(PCR_Rd, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13"))) %>%
  mutate(Soil_mu = factor(Soil_mu, levels = c("100", "101","108","138","145","168","126","127"))) %>%
  mutate(Trt = factor(Trt, levels=c("LF","MF")))

###Distribution of Data (non-parametric)###
#From the above graphs, it is clear that the Total number of OTUs per sample is 
#NOT normally distributed, but let's statistically test it anyways
shapiro.test(d$Total) #p < 0.001; Total is not normally distributed
#what if the Total OTUs is log transformed? 
p = ggplot(d, aes(log(Total))) + theme_bw()
p = p + geom_histogram(binwidth=0.05)
p = p + labs(x="Log Total OTUs", y="Sample Count", title = "All Samples Log Total OTUs Silva")
p
#looks like a bimodal distribution so once again, not normally distributed
shapiro.test(log(d$Total)) #p<0.001; Log(Total) is not normally distributed

#We need to use non-parametric statistical tests moving forward! Instead of 
#one-way ANOVAs to test for block effects, we will use the Kruskal-Wallis test.
#For post-hoc tests we will use Wilcoxon test with a Bonferroni correction.


#######Is there a systematic bias of Unit that can explain variance of total OTUs? 
p = ggplot(d, aes(x=Unit, y=Total)) + theme_bw() 
p = p + geom_boxplot() + scale_x_discrete(labels = c("Butte Lake", "Hole", "Unburned", "Warner Valley"))
p = p + labs(x="Unit", y="Total OTUs per Sample", title = "Total OTUs by Unit")
p = p + annotate("text", x=4, y=250000, label = "p > 0.05")
p

#ask the question with statistical power
kruskal.test(Total~Unit, data=d) # p < 0.05
# To determine which units differ, use a Wilcoxon test w/ Bonferroni correction 
pairwise.wilcox.test(d$Total, d$Unit, p.adjust.method = "bonferroni")
# w/ Bonferroni correction, there are no pairwise comparisons that are statistically
# significant (p<0.05).


#########By Horizon
p = ggplot(d, aes(x=Horizon, y=Total)) + theme_bw() 
p = p + geom_boxplot()+ scale_x_discrete(labels = c("Mineral", "Organic"))
p = p + labs(x="Horizon", y="Total OTUs per Sample", title = "Total OTUs by Horizon")
p = p + annotate("text", x=0.8, y=250000, label = "* p < 0.05")
p

#ask the question with statistical power
kruskal.test(Total~Horizon, data=d) # p = 0.018

########By Vegetation Type
p = ggplot(d, aes(x=Veg_Type, y=Total)) + theme_bw() 
p = p + geom_boxplot()
p = p + labs(x="Vegetation Type", y="Total OTUs per Sample", title = "Total OTUs by Vegetation Type")
p = p + annotate("text", x=5, y=250000, label = " p > 0.05")
p

#ask the question with statistical power
kruskal.test(Total~Veg_Type, data=d) #p>0.05


#By Soil Type
p = ggplot(d, aes(x=Soil_mu, y=Total)) + theme_bw() 
p = p + geom_boxplot()
p = p + labs(x="Soil Type", y="Total OTUs per Sample", title = "Total OTUs by Soil Type")
p = p + annotate("text", x=7, y=250000, label = " p > 0.05")
p

#ask the question with statistical power
kruskal.test(Total~Soil_mu, data=d) # p>0.05

#######By Treatment (Less Fire v. More Fire)
p = ggplot(d, aes(x=Trt, y=Total)) + theme_bw() 
p = p + geom_boxplot() + scale_x_discrete(labels = c("Less Fire","More Fire"))
p = p + labs(x="Treatment", y="Total OTUs per Sample", title = "Total OTUs by Treatment")
p = p + annotate("text", x=1, y=250000, label = " p > 0.05")
p

#ask the question with statistical power
kruskal.test(Total~Trt, data=d) # p = 0.098


#########Weight Round
p = ggplot(d, aes(x=Wt_Rd, y=Total)) + theme_bw() + geom_boxplot()
p = p + labs(x="Weighing Round", y="Total OTUs per Sample", title = "Total OTUs by Weighing Round")
p = p + annotate("text", x=16, y=250000, label = "*** p < 0.001")
p

#statistical test
kruskal.test(Total~Wt_Rd, data=d) # p < 0.001
pairwise.wilcox.test(d$Total, d$Wt_Rd, p.adjust.method = "bonferroni") 
#seems to be driven by Rounds 15 & 16 being different than the rest


########DNA Extraction Round
#pay close attention to round 1 because negative control showed high total reads, 
#which suggests contamination
p = ggplot(d, aes(x=Extract_Rd, y=Total)) + theme_bw() + geom_boxplot()
p = p + labs(x="DNA Extraction Round", y="Total OTUs per Sample", title = "Total OTUs by DNA Extraction")
p = p + annotate("text", x=28, y=250000, label = "*** p < 0.001")
p

#kruskal test
kruskal.test(Total~Extract_Rd, data=d) # p < 0.001
#posthoc test
pairwise.wilcox.test(d$Total, d$Extract_Rd, p.adjust.method = "bonferroni")

##########PCR Round
p = ggplot(d, aes(x=PCR_Rd, y=Total)) + theme_bw() + geom_boxplot()
p = p + labs(x="PCR Round", y="Total OTUs per Sample", title = "Total OTUs by PCR Round")
p = p + annotate("text", x=12, y=250000, label = "*** p < 0.001")
p

#kruskal test
kruskal.test(Total~PCR_Rd, data=d) #p < 0.001
#post hoc pairwise comparison test
pairwise.wilcox.test(d$Total, d$PCR_Rd, p.adjust.method = "bonferroni")

#########Sequal Plate 
p = ggplot(d, aes(x=Sequal_Rd, y=Total)) + theme_bw() + geom_boxplot()
p = p + labs(x="Sequal Plate", y="Total OTUs per Sample", title = "Total OTUs by Sequal Plate")
p = p + annotate("text", x=4, y=250000, label = "*** p < 0.001")
p

#kruskal test
kruskal.test(Total~Sequal_Rd, data=d) # p<0.001
#post hoc test 
pairwise.wilcox.test(d$Total, d$Sequal_Rd, p.adjust.method = "bonferroni")