#Objective: To randomly select samples from the three units to include in 
#PERMANOVAs to test for significant effects of fire history while preventing
#pseudo-replication inflating significance values

#################Loading libraries and phyloseq object########################

library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(vegan)
library(patchwork)
library(viridis)

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

#make Soil type a factor and sample severity an ordinal variable
sample_data(hell.ps.3)$Soil_mu = as.factor(sample_data(hell.ps.3)$Soil_mu)
sample_data(hell.ps.3)$SampleSev <- factor(sample_data(hell.ps.3)$SampleSev, order=TRUE, levels=c(1, 2, 3, 4))
class(sample_data(hell.ps.3)$SampleSev)


rm(ps.3, cohort_df)

#check sample data 
head(sample_data(hell.ps.3))

# 3. When running a PERMANOVA, GustaMe reports:
# "Anderson (2001) warns that groups of objects with different dispersions, 
# yet no significant differences in centres (centres are similar to means, 
# but may be non-Euclidean), may result in misleadingly low P-values. 
# It is thus recommended that the dispersion be evaluated and considered when 
# interpreting the results of NPMANOVA. See Anderson (2006) for a discussion 
# on tests of multivariate dispersion."

# Thus, we should also look at the dispersion of our samples.
b.Trt = betadisper(dist.hell, SamDat$Trt)
b.Trt
permutest(b.Trt) 

################# PCoA ORDINATIONS #######################################
############Create separate ordinations per Unit
######BL
# Create ordination
ps.BL = subset_samples(hell.ps.3, Unit == "BL")
hell.pcoa.ord.BL = ordinate(ps.BL, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(hell.pcoa.ord.BL$vectors)$Axis.1
y = data.frame(hell.pcoa.ord.BL$vectors)$Axis.2
df.ord.hell.BL = data.frame(sample_data(ps.BL))
df.ord.hell.BL$PCoA1 = x
df.ord.hell.BL$PCoA2 = y
df.ord.hell.BL$Horizon = factor(df.ord.hell.BL$Horizon, levels = c("O", "M"))


#get values for axes
Var.Axis.1 = hell.pcoa.ord.BL$values$Relative_eig[1]
Var.Axis.2 = hell.pcoa.ord.BL$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

#rename factors for the key
levels(df.ord.hell.BL$Horizon) <- list(Organic = "O", Mineral="M")
levels(df.ord.hell.BL$SampleSev) <-list(`Unburned` = "1", `Low`="2", `Moderate`="3", `Severe`=4)

#Color by treatment
a = ggplot(df.ord.hell.BL)
a = a + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon), size=3)
a = a + theme_bw() + labs(title="Butte Lake", x=Axis1.label, y=Axis2.label)
a = a + scale_color_manual(values=c('#0097d3','#ffbb22')) 
a = a + scale_shape_manual(values = c(15,17)) + theme(legend.position = "none")
a

#Color by burn severity (sample severity)
b = ggplot(df.ord.hell.BL)
b = b + geom_point(aes(x=PCoA1,y=PCoA2,color=SampleSev,shape=Horizon), size=3)
b = b + theme_bw() + labs(title="Butte Lake", x=Axis1.label, y=Axis2.label)
b = b + scale_color_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'))
b = b + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none")
b

#Color by pH
h = ggplot(df.ord.hell.BL)
h = h + geom_point(aes(x=PCoA1,y=PCoA2,color=pH,shape=Horizon), size=3)
h = h + theme_bw() + labs(title="Butte Lake", x=Axis1.label, y=Axis2.label)
h = h + scale_color_viridis(limits=c(3, 9), option="D", direction=-1)
h = h + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none") 
h

a + b


########### H
# Create ordination
ps.H = subset_samples(hell.ps.3, Unit == "H")
hell.pcoa.ord.H = ordinate(ps.H, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(hell.pcoa.ord.H$vectors)$Axis.1
y = data.frame(hell.pcoa.ord.H$vectors)$Axis.2
df.ord.hell.H = data.frame(sample_data(ps.H))
df.ord.hell.H$PCoA1 = x
df.ord.hell.H$PCoA2 = y
df.ord.hell.H$Horizon = factor(df.ord.hell.H$Horizon, levels = c("O", "M"))


#get values for axes
Var.Axis.1 = hell.pcoa.ord.H$values$Relative_eig[1]
Var.Axis.2 = hell.pcoa.ord.H$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

#Color by treatment
c = ggplot(df.ord.hell.H)
c = c + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon), size=3)
c = c + theme_bw() + labs(title="Hole", x=Axis1.label, y=Axis2.label)
c = c + scale_color_manual(values=c('#0097d3','#ffbb22')) 
c = c + scale_shape_manual(values = c(15,17)) + theme(legend.position = "none")
c

d = ggplot(df.ord.hell.H)
d = d + geom_point(aes(x=PCoA1,y=PCoA2,color=SampleSev,shape=Horizon), size=3)
d = d + theme_bw() + labs(title="Hole", x=Axis1.label, y=Axis2.label)
d = d + scale_color_manual(name="Burn Severity", values=c('#f5db4c', '#f98e09', '#bc3754', '#57106e'))
d = d + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none") 
d

#color by pH
i = ggplot(df.ord.hell.H)
i = i + geom_point(aes(x=PCoA1,y=PCoA2,color=pH,shape=Horizon), size=3)
i = i + theme_bw() + labs(title="Hole", x=Axis1.label, y=Axis2.label)
i = i + scale_color_viridis(limits=c(3, 9), option="D", direction=-1)
i = i + scale_shape_manual(values=c(15,17)) + theme(legend.position="none")
i

c+d

########### WV
# Create ordination
ps.WV = subset_samples(hell.ps.3, Unit == "WV")
hell.pcoa.ord.WV = ordinate(ps.WV, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(hell.pcoa.ord.WV$vectors)$Axis.1
y = data.frame(hell.pcoa.ord.WV$vectors)$Axis.2
df.ord.hell.WV = data.frame(sample_data(ps.WV))
df.ord.hell.WV$PCoA1 = x
df.ord.hell.WV$PCoA2 = y
df.ord.hell.WV$Horizon = factor(df.ord.hell.WV$Horizon, levels = c("O", "M"))


#get values for axes
Var.Axis.1 = hell.pcoa.ord.WV$values$Relative_eig[1]
Var.Axis.2 = hell.pcoa.ord.WV$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

#rename factors for the key
levels(df.ord.hell.WV$Horizon) <- list(Organic = "O", Mineral="M")
levels(df.ord.hell.WV$Trt) <-list(`Less Fire` = "LF", `More Fire`="MF")

#Color by treatment
e = ggplot(df.ord.hell.WV)
e = e + geom_point(aes(x=PCoA1,y=PCoA2,color=factor(Trt),shape=Horizon), size=3)
e = e + theme_bw() + labs(title="Warner Valley", x=Axis1.label, y=Axis2.label)
e = e + scale_color_manual(name="Fire History", values=c('#0097d3', '#ffbb22')) 
e = e + scale_shape_manual(values = c(15,17)) + theme(legend.position = "none")
e

#color by burn severity
f = ggplot(df.ord.hell.WV)
f = f + geom_point(aes(x=PCoA1,y=PCoA2,color=SampleSev,shape=Horizon), size=3)
f = f + theme_bw() + labs(title="Warner Valley", x=Axis1.label, y=Axis2.label)
f = f + scale_color_manual(name="Burn Severity", values=c( '#f98e09', '#bc3754', '#57106e'))
f = f + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none")
f

#color by pH
j = ggplot(df.ord.hell.WV)
j = j + geom_point(aes(x=PCoA1,y=PCoA2,color=pH,shape=Horizon), size=3)
j = j + theme_bw() + labs(title="Warner Valley", x=Axis1.label, y=Axis2.label)
j = j + scale_color_viridis(limits=c(3, 9), option="D", direction=-1)
j = j + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none")
j


a + c + e
b+d+f

h+i+j

#####################Create dataframe of cohorts################################
####ONLY DO IF WANT TO ADJUST COHORT #s TO SOMETHING OTHER THAN 500!!!
#get number of plots 
temp <- data.frame(sample_data(hell.ps.3)) %>%
  distinct(Plot_Name, .keep_all=TRUE)
plots <- temp$Plot_Name
length(plots) #58 plots

rm(temp)

#how many cohorts (group of randomly selected plots) do you want to test
n <- 100

#create dataframe of the samples randomly selected for each cohort  
cohort_df <- data.frame(Plot_Name=character(), Rep=character(), Cohort=numeric())

for(i in 1:n) {
  set.seed(i)
  my_list = list(Plot_Name=plots, Rep = sample(c("A", "B", "C"), length(plots), replace = TRUE), Cohort=replicate(n=length(plots), i))
  cohort_df = rbind(cohort_df, my_list, stringsAsFactors=FALSE)
}
rm(i, my_list, plots, n)

#create columns to subset out phyloseq object
cohort_df$Sample_Core_Name <- paste(cohort_df$Plot_Name, cohort_df$Rep, sep="_")
cohort_df$O_Samples <- paste(cohort_df$Sample_Core_Name, "O", sep="") 
cohort_df$O_Samples[!(cohort_df$O_Samples %in% as.data.frame(sample_data(hell.ps.3))$Full_Sample_Name)] = "NA" #replaces O samples that don't exist in the dataset with "NA"
cohort_df$M_Samples <- paste(cohort_df$Sample_Core_Name, "M", sep="")
cohort_df$M_Samples[!(cohort_df$M_Samples %in% as.data.frame(sample_data(hell.ps.3))$Full_Sample_Name)] = "NA" #replaces M samples that don't exist in the dataset with "NA"
cohort_df_long <- gather(cohort_df, type, Sample, O_Samples:M_Samples, factor_key=TRUE) %>%
  filter(., Sample!="NA")

rm(cohort_df)



############500 Cohorts Trt and Trt+bs #################################
#3/13/24 edit treat burn severity as categorical; ordinal 
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
  
  ###SampleSev + Trt
  ###Trt + SampleSev
  #Butte Lake
  BL.bs.M = adonis2(dist.BL.M ~ SampleSev + Trt, SamDat.BL.M)
  BL.bs.M$term=row.names(BL.bs.M)
  BL.bs.M = BL.bs.M %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  BL.bs.O = adonis2(dist.BL.O ~ SampleSev + Trt, SamDat.BL.O)
  BL.bs.O$term=row.names(BL.bs.O)
  BL.bs.O = BL.bs.O %>%
    mutate(cohort=i, covars="bs+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.bs.M = adonis2(dist.H.M ~ SampleSev + Trt, SamDat.H.M)
  H.bs.M$term=row.names(H.bs.M)
  H.bs.M = H.bs.M %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  H.bs.O = adonis2(dist.H.O ~ SampleSev + Trt, SamDat.H.O)
  H.bs.O$term=row.names(H.bs.O)
  H.bs.O = H.bs.O %>%
    mutate(cohort=i, covars="bs+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.bs.M = adonis2(dist.WV.M ~ SampleSev + Trt, SamDat.WV.M)
  WV.bs.M$term = row.names(WV.bs.M)
  WV.bs.M = WV.bs.M %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  WV.bs.O = adonis2(dist.WV.O ~ SampleSev+Trt, SamDat.WV.O)
  WV.bs.O$term = row.names(WV.bs.O)
  WV.bs.O = WV.bs.O %>%
    mutate(cohort=i, covars="bs+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #bind together outputs of PERMANOVAs
  b = rbind(BL.M, BL.O, H.M, H.O, WV.M, WV.O, BL.bs.M, BL.bs.O, H.bs.M, H.bs.O, WV.bs.M, WV.bs.O)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_bsordinal_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(unit, horizon, covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1)) %>%
  filter(covars=="trt")

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

#write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_trtandtrtbs_bsordinal_summary.csv")

#plot histogram of p-values
trt.trtbs.raw <- read.csv("~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_bsordinal_rawresults.csv")

bs.p.values <- filter(trt.trtbs.raw, p!="NA") %>%
  filter(., covars=="bs+trt")

#bs histogram p-values
a <- ggplot(data=bs.p.values[bs.p.values$unit=="BL" & bs.p.values$term=="SampleSev", ], aes(x=p)) + geom_histogram() + facet_wrap(~horizon) + geom_vline(aes(xintercept=0.05)) + labs(title="Butte Lake")
b <- ggplot(data=bs.p.values[bs.p.values$unit=="H" & bs.p.values$term=="SampleSev", ], aes(x=p)) + geom_histogram() + facet_wrap(~horizon) + geom_vline(aes(xintercept=0.05)) + labs(title="Hole")
c <- ggplot(data=bs.p.values[bs.p.values$unit=="WV" & bs.p.values$term=="SampleSev", ], aes(x=p)) + geom_histogram() + facet_wrap(~horizon) + geom_vline(aes(xintercept=0.05)) + labs(title="Warner Valley")

a+b+c

#median p-values for bs
bs.p.values.sum <- bs.p.values %>%
  group_by(term, unit, horizon) %>%
  summarise(median.p=median(p))


############## 500 cohorts PERMANOVAs with horizon as strata ###################

final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort i
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  ps.BL <- subset_samples(ps.cohort, Unit=="BL")
  ps.H <- subset_samples(ps.cohort, Unit=="H")
  ps.WV <- subset_samples(ps.cohort, Unit=="WV")
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat.BL = data.frame(sample_data(ps.BL))
  SamDat.H = data.frame(sample_data(ps.H))
  SamDat.WV = data.frame(sample_data(ps.WV))
  
  #make dissimilarity matrix
  dist.BL = phyloseq::distance(ps.BL, method="bray")
  dist.H = phyloseq::distance(ps.H, method="bray")
  dist.WV = phyloseq::distance(ps.WV, method="bray")
  
  ###only Trt
  #Butte Lake
  BL = adonis2(dist.BL ~ Trt, strata=SamDat.BL$Horizon, data=SamDat.BL)
  BL$term=row.names(BL)
  BL = BL %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H = adonis2(dist.H ~ Trt, strata=SamDat.H$Horizon, data=SamDat.H)
  H$term=row.names(H)
  H = H %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV = adonis2(dist.WV ~ Trt, strata=SamDat.WV$Horizon, data=SamDat.WV)
  WV$term = row.names(WV)
  WV = WV %>%
    mutate(cohort=i, covars="trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###SampleSev + Trt
  ###Trt + SampleSev
  #Butte Lake
  BL.bs = adonis2(dist.BL ~ SampleSev + Trt, strata=SamDat.BL$Horizon, data=SamDat.BL)
  BL.bs$term=row.names(BL.bs)
  BL.bs = BL.bs %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Hole
  H.bs = adonis2(dist.H ~ SampleSev + Trt, strata=SamDat.H$Horizon, data=SamDat.H)
  H.bs$term=row.names(H.bs)
  H.bs = H.bs %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #Warner Valley
  WV.bs = adonis2(dist.WV ~ SampleSev + Trt, strata=SamDat.WV$Horizon, data=SamDat.WV)
  WV.bs$term = row.names(WV.bs)
  WV.bs = WV.bs %>%
    mutate(cohort=i, covars="bs+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #bind together outputs of PERMANOVAs
  b = rbind(BL, H, WV, BL.bs, H.bs, WV.bs)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

#write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_stratahorizon_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(unit, horizon, covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1)) %>%
  filter(covars=="trt")

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_stratahorizon_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))

#write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_stratahorizon_summary.csv")

############## 500 cohorts horizon and unit as fixed effects ################
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){
  #list of Full Name of Mineral samples for cohort i
  cohort_samples <- cohort_df_long[cohort_df_long$Cohort==i, ]
  
  
  #create 2 subsetted phyloseq objects for each cohort; one for mineral one for organic
  ps.cohort <- subset_samples(hell.ps.3, Full_Sample_Name %in% cohort_samples$Sample)
  
  # Get sample_data as dataframe for each phyloseq object
  SamDat = data.frame(sample_data(ps.cohort))

  #make dissimilarity matrix
  dist = phyloseq::distance(ps.cohort, method="bray")
  
  ###only Trt
  adon.i = adonis2(dist ~ Unit + Horizon + Trt, data=SamDat)
  adon.i$term=row.names(adon.i)
  adon.i = adon.i %>%
    mutate(cohort=i, covars="trt", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  ###SampleSev + Trt
  adon.bs = adonis2(dist ~ Unit + Horizon +SampleSev + Trt, data=SamDat)
  adon.bs$term=row.names(adon.bs)
  adon.bs = adon.bs %>%
    mutate(cohort=i, covars="bs+trt",p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #bind together outputs of PERMANOVAs
  b = rbind(adon.i, adon.bs)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
}

###Summarise Outputs

write.csv(final_df, "~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_horizonunitfixed_rawresults.csv")

cohor500_df_summary <- final_df %>%
  filter(term=="Total") %>%
  group_by(covars) %>%
  summarise(mean.n=mean(Df+1), min.n=min(Df+1), max.n=max(Df+1), sd.n=sd(Df+1)) %>%
  filter(covars=="trt")

#write.csv(cohor500_df_summary, "~/Desktop/LAVO22_GG2/500cohorts_horizonunitfixed_samplesummary.csv")

cohor500.sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(covars, term, p.sig) %>%
  summarise(n=n(), perc=n/500, mean.R2=mean(R2), sd.R2=sd(R2))
#write.csv(cohor500.sum, "~/Desktop/LAVO22_GG2/500cohorts_trtandbstrt_horizonunitfixed_summary.csv")

############### PERMANOVAs including Soil and Veg type as covariates #########
###############Subset phyloseq object and run permanovas######################
#all 6 PERMANOVAs separately
#empty dataframe to store coefficients
final_df = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())


for(i in 1:max(cohort_df_long$Cohort)){

  ############Empty adonis outputs 
  BL.all.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  BL.all.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  BL.allbs.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  BL.allbs.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  BL.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  BL.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  
  H.all.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  H.all.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  H.allbs.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  H.allbs.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  H.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  H.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  
  WV.all.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  WV.all.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  WV.allbs.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  WV.allbs.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  WV.M = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  WV.O = data.frame(Df=numeric(), SumOfSqs=numeric(), R2=numeric(), F=numeric(), term=character(), cohort=numeric(), covars=character(), horizon=character(), unit=character(), p=numeric())
  
  
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
  
  ###PERMANOVAs
  ##Butte Lake
  #variables of interest (soil and veg and trt)
  if(length(unique(factor(SamDat.BL.M$Veg_Type)))>1 & length(unique(factor(SamDat.BL.M$Soil_mu)))>1){
  BL.all.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + Trt, SamDat.BL.M)
  BL.all.M$term=row.names(BL.all.M)
  BL.all.M = BL.all.M %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  if(length(unique(factor(SamDat.BL.O$Veg_Type)))>1 & length(unique(factor(SamDat.BL.O$Soil_mu)))>1){
  BL.all.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + Trt,  SamDat.BL.O)
  BL.all.O$term=row.names(BL.all.O)
  BL.all.O = BL.all.O %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  #all variables of interest (soil, veg, burn severity, and trt)
  if(length(unique(factor(SamDat.BL.M$Veg_Type)))>1 & length(unique(factor(SamDat.BL.M$Soil_mu)))>1){
  BL.allbs.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.M)
  BL.allbs.M$term=row.names(BL.allbs.M)
  BL.allbs.M = BL.allbs.M %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  if(length(unique(factor(SamDat.BL.O$Veg_Type)))>1 & length(unique(factor(SamDat.BL.O$Soil_mu)))>1){
  BL.allbs.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.O)
  BL.allbs.O$term=row.names(BL.allbs.O)
  BL.allbs.O = BL.allbs.O %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  #only Trt
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

  
  
  ##Hole
  #variables of interest (soil and veg and trt)
  if(length(unique(factor(SamDat.H.M$Veg_Type)))>1 & length(unique(factor(SamDat.H.M$Soil_mu)))>1){
  H.all.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + Trt, SamDat.H.M)
  H.all.M$term=row.names(H.all.M)
  H.all.M = H.all.M %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #add in burn severity 
  H.allbs.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.M)
  H.allbs.M$term=row.names(H.allbs.M)
  H.allbs.M = H.allbs.M %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  if(length(unique(factor(SamDat.H.O$Veg_Type)))>1 & length(unique(factor(SamDat.H.O$Soil_mu)))>1){
  H.all.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + Trt,  SamDat.H.O)
  H.all.O$term=row.names(H.all.O)
  H.all.O = H.all.O %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #all variables of interest (soil, veg, burn severity, and trt)
  H.allbs.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.O)
  H.allbs.O$term=row.names(H.allbs.O)
  H.allbs.O = H.allbs.O %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  #only Trt
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

  
  ###Warner Valley
  #variables of interest (soil and veg and trt)
  if(length(unique(factor(SamDat.WV.M$Veg_Type)))>1 & length(unique(factor(SamDat.WV.M$Soil_mu)))>1){
  WV.all.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + Trt, SamDat.WV.M)
  WV.all.M$term=row.names(WV.all.M)
  WV.all.M = WV.all.M %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #add burn severity as a term
  WV.allbs.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.M)
  WV.allbs.M$term=row.names(WV.allbs.M)
  WV.allbs.M = WV.allbs.M %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  if(length(unique(factor(SamDat.WV.O$Veg_Type)))>1 & length(unique(factor(SamDat.WV.O$Soil_mu)))>1){
  WV.all.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + Trt,  SamDat.WV.O)
  WV.all.O$term=row.names(WV.all.O)
  WV.all.O = WV.all.O %>%
    mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  
  #all variables of interest (soil, veg, burn severity, and trt)
  WV.allbs.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.O)
  WV.allbs.O$term=row.names(WV.allbs.O)
  WV.allbs.O = WV.allbs.O %>%
    mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
    select(-`Pr(>F)`)
  }
  
  #only Trt
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
  
  ###combine all permanova outputs from the cohort into one table
  b = rbind(BL.all.M, BL.all.O, BL.allbs.M, BL.allbs.O, BL.M, BL.O, H.all.M, H.all.O, H.allbs.M, H.allbs.O, H.M, H.O, WV.all.M, WV.all.O, WV.allbs.M, WV.allbs.O, WV.M, WV.O)
  
  rm(BL.all.M, BL.all.O, BL.allbs.M, BL.allbs.O, BL.M, BL.O, H.all.M, H.all.O, H.allbs.M, H.allbs.O, H.M, H.O, WV.all.M, WV.all.O, WV.allbs.M, WV.allbs.O, WV.M, WV.O)
  
  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)

  #########Now that we have the permanovas made, let's store the coefficients
  
  #for loop to create permanovas and store coefficients in dataframe
  
  #for loop to make PERMANOVAs for a single cohort 
  #(1:4 because 4 burn severity metrics to compare to start)
  
  #make dataframe containing all of the outputs of the PERMANOVAs
  #can't start for loop from 0 so need to make an empty dataframe size 5x3 to 
  #include at the beginning 
  #temp = data.frame(replicate(6, c(0,0,0)))
  #names(temp) <- names(BL.all.M)
  
#  b = rbind(BL.all.M, BL.all.O, BL.allbs.M, BL.allbs.O, BL.M, BL.O, H.all.M, H.all.O, H.allbs.M, H.allbs.O, H.M, H.O, WV.all.M, WV.all.O, WV.allbs.M, WV.allbs.O, WV.M, WV.O)
  
  #add to the final dataframe
#  final_df <- rbind(final_df, b, stringsAsFactors=FALSE)
  
#  a = data.frame(cohort=numeric(), unit=character(), horizon=character(), terms.included=character(), term=character(), total_df=numeric(),R2=numeric(), F=numeric(), p_fdr=numeric())
  
#  for(x in 1:18) {
#    a[x, 1] = i #cohort number
#    a[x, 2] = ifelse(x<7, "BL", ifelse(x>12, "WV", "H"))
#    a[x, 2] = ifelse(x %% 2==0, 'O', 'M')
#    a[x, 3] = c[i]
#    a[x, 3] = b$row[x*3+1]
#    a[x, 4] = b$Df[x*3+2]
#    a[x, 5] = b$R2[x*3+1]
#    a[x, 6] = b$F[x*3+1]
#    a[x, 7] = b$`Pr(>F)`[x*3+1]
    
#  }
#  final_df <- rbind(final_df, a, stringsAsFactors=FALSE)
  
}

################# summarise outputs ############
final_sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/3, mean.R2=mean(R2), sd.R2=sd(R2))

#write.csv(final_sum, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/FireHistoryPERMANOVAs/3cohorts_allcovars_results.csv")

############### results treating all samples as independent #############
#subset phyloseq objects
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

###PERMANOVAs
##Butte Lake
#variables of interest (soil and veg and trt)
BL.all.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + Trt, SamDat.BL.M)
BL.all.M$term=row.names(BL.all.M)
BL.all.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + Trt,  SamDat.BL.O)
BL.all.O$term=row.names(BL.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
BL.allbs.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.M)
BL.allbs.M$term=row.names(BL.allbs.M)
BL.allbs.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.O)
BL.allbs.O$term=row.names(BL.allbs.O)
#only Trt
BL.M = adonis2(dist.BL.M ~ Trt, SamDat.BL.M)
BL.M$term=row.names(BL.M)
BL.O = adonis2(dist.BL.O ~ Trt, SamDat.BL.O)
BL.O$term=row.names(BL.O)

BL.all.M = BL.all.M %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.all.O = BL.all.O %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.allbs.M = BL.allbs.M %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.allbs.O = BL.allbs.O %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.M = BL.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.O = BL.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)



##Hole
#variables of interest (soil and veg and trt)
H.all.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + Trt, SamDat.H.M)
H.all.M$term=row.names(H.all.M)
H.all.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + Trt,  SamDat.H.O)
H.all.O$term=row.names(H.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
H.allbs.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.M)
H.allbs.M$term=row.names(H.allbs.M)
H.allbs.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.O)
H.allbs.O$term=row.names(H.allbs.O)
#only Trt
H.M = adonis2(dist.H.M ~ Trt, SamDat.H.M)
H.M$term=row.names(H.M)
H.O = adonis2(dist.H.O ~ Trt, SamDat.H.O)
H.O$term=row.names(H.O)

H.all.M = H.all.M %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.all.O = H.all.O %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.allbs.M = H.allbs.M %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.allbs.O = H.allbs.O %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.M = H.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.O = H.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Warner Valley
#variables of interest (soil and veg and trt)
WV.all.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + Trt, SamDat.WV.M)
WV.all.M$term=row.names(WV.all.M)
WV.all.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + Trt,  SamDat.WV.O)
WV.all.O$term=row.names(WV.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
WV.allbs.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.M)
WV.allbs.M$term=row.names(WV.allbs.M)
WV.allbs.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.O)
WV.allbs.O$term=row.names(WV.allbs.O)
#only Trt
WV.M = adonis2(dist.WV.M ~ Trt, SamDat.WV.M)
WV.M$term = row.names(WV.M)
WV.O = adonis2(dist.WV.O ~ Trt, SamDat.WV.O)
WV.O$term = row.names(WV.O)

WV.all.M = WV.all.M %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.all.O = WV.all.O %>%
  mutate(cohort=i, covars="soil+veg+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.allbs.M = WV.allbs.M %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.allbs.O = WV.allbs.O %>%
  mutate(cohort=i, covars="soil+veg+bs+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.M = WV.M %>%
  mutate(cohort=i, covars="trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.O = WV.O %>%
  mutate(cohort=i, covars="trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

adonis.all <- rbind(BL.all.M, BL.all.O, BL.allbs.M, BL.allbs.O, BL.M, BL.O, H.all.M, H.all.O, H.allbs.M, H.allbs.O, H.M, H.O, WV.all.M, WV.all.O, WV.allbs.M, WV.allbs.O, WV.M, WV.O)

#look at only trt effects
all.Trt.results <- filter(adonis.all, term=="Trt")
#save results
write.csv(adonis.all, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/FireHistoryPERMANOVAs/allsamples_allcovars_results.csv")

################# summarise outputs ############
final_sum <- final_df %>%
  filter(p!="NA") %>%
  mutate(p.range=ifelse(p>0.05, "p>0.05", ifelse(p<=0.001, "p<0.001", "p<0.05")), p.sig=ifelse(p<0.05, "Y", "N")) %>%
  group_by(unit, horizon, covars, term, p.sig) %>%
  summarise(n=n(), perc=n/3, mean.R2=mean(R2), sd.R2=sd(R2))

#write.csv(final_sum, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/FireHistoryPERMANOVAs/3cohorts_allcovars_results.csv")

################ PERMANOVAs from averaged OTUs ##############################
#load averaged OTU ps object
ps.avg <- readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3.avg")

#hellinger transformation for PERMANOVAs
hell.ps.avg = transform_sample_counts(ps.avg, function(x) (x / sum(x))^0.5 )

#make Soil type a factored variable instead of integer
sample_data(hell.ps.avg)$Soil_mu = as.factor(sample_data(hell.ps.avg)$Soil_mu)

###make PERMANOVAs
ps.BL.O <- subset_samples(hell.ps.avg, Unit=="BL" & Horizon=="O")
ps.H.O <- subset_samples(hell.ps.avg, Unit=="H" & Horizon=="O")
ps.WV.O <- subset_samples(hell.ps.avg, Unit=="WV"& Horizon=="O")
ps.BL.M <- subset_samples(hell.ps.avg, Unit=="BL" & Horizon=="M")
ps.H.M <- subset_samples(hell.ps.avg, Unit=="H" & Horizon=="M")
ps.WV.M <- subset_samples(hell.ps.avg, Unit=="WV"& Horizon=="M")



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

###PERMANOVAs
##Butte Lake
#variables of interest (soil and veg and trt)
BL.all.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + Trt, SamDat.BL.M)
BL.all.M$term=row.names(BL.all.M)
BL.all.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + Trt,  SamDat.BL.O)
BL.all.O$term=row.names(BL.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
#BL.allbs.M = adonis2(dist.BL.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.M)
#BL.allbs.M$term=row.names(BL.allbs.M)
#BL.allbs.O = adonis2(dist.BL.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.BL.O)
#BL.allbs.O$term=row.names(BL.allbs.O)

#only Trt
BL.M = adonis2(dist.BL.M ~ Trt, SamDat.BL.M)
BL.M$term=row.names(BL.M)
BL.O = adonis2(dist.BL.O ~ Trt, SamDat.BL.O)
BL.O$term=row.names(BL.O)

BL.all.M = BL.all.M %>%
  mutate(covars="soil+veg+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.all.O = BL.all.O %>%
  mutate(covars="soil+veg+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
#BL.allbs.M = BL.allbs.M %>%
#  mutate(covars="soil+veg+bs+trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
#BL.allbs.O = BL.allbs.O %>%
#  mutate(covars="soil+veg+bs+trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
BL.M = BL.M %>%
  mutate(covars="trt", horizon="M", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
BL.O = BL.O %>%
  mutate(covars="trt", horizon="O", unit="BL", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)



##Hole
#variables of interest (soil and veg and trt)
H.all.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + Trt, SamDat.H.M)
H.all.M$term=row.names(H.all.M)
H.all.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + Trt,  SamDat.H.O)
H.all.O$term=row.names(H.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
#H.allbs.M = adonis2(dist.H.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.M)
#H.allbs.M$term=row.names(H.allbs.M)
#H.allbs.O = adonis2(dist.H.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.H.O)
#H.allbs.O$term=row.names(H.allbs.O)
#only Trt
H.M = adonis2(dist.H.M ~ Trt, SamDat.H.M)
H.M$term=row.names(H.M)
H.O = adonis2(dist.H.O ~ Trt, SamDat.H.O)
H.O$term=row.names(H.O)

H.all.M = H.all.M %>%
  mutate(covars="soil+veg+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.all.O = H.all.O %>%
  mutate(covars="soil+veg+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
#H.allbs.M = H.allbs.M %>%
#  mutate(covars="soil+veg+bs+trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
#H.allbs.O = H.allbs.O %>%
#  mutate(covars="soil+veg+bs+trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
H.M = H.M %>%
  mutate(covars="trt", horizon="M", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
H.O = H.O %>%
  mutate(covars="trt", horizon="O", unit="H", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

###Warner Valley
#variables of interest (soil and veg and trt)
WV.all.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + Trt, SamDat.WV.M)
WV.all.M$term=row.names(WV.all.M)
WV.all.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + Trt,  SamDat.WV.O)
WV.all.O$term=row.names(WV.all.O)
#all variables of interest (soil, veg, burn severity, and trt)
#WV.allbs.M = adonis2(dist.WV.M ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.M)
#WV.allbs.M$term=row.names(WV.allbs.M)
#WV.allbs.O = adonis2(dist.WV.O ~ Soil_mu + Veg_Type + SampleSev + Trt, SamDat.WV.O)
#WV.allbs.O$term=row.names(WV.allbs.O)
#only Trt
WV.M = adonis2(dist.WV.M ~ Trt, SamDat.WV.M)
WV.M$term = row.names(WV.M)
WV.O = adonis2(dist.WV.O ~ Trt, SamDat.WV.O)
WV.O$term = row.names(WV.O)

WV.all.M = WV.all.M %>%
  mutate(covars="soil+veg+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.all.O = WV.all.O %>%
  mutate(covars="soil+veg+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
#WV.allbs.M = WV.allbs.M %>%
#  mutate(covars="soil+veg+bs+trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
#WV.allbs.O = WV.allbs.O %>%
#  mutate(covars="soil+veg+bs+trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
#  select(-`Pr(>F)`)
WV.M = WV.M %>%
  mutate(covars="trt", horizon="M", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)
WV.O = WV.O %>%
  mutate(covars="trt", horizon="O", unit="WV", p=`Pr(>F)`) %>%
  select(-`Pr(>F)`)

adonis.all <- rbind(BL.all.M, BL.all.O, BL.M, BL.O, H.all.M, H.all.O, H.M, H.O, WV.all.M, WV.all.O, WV.M, WV.O)


#save results
write.csv(adonis.all, "~/Desktop/LAVO22_GG2/avgsamples_allcovars_results.csv")




