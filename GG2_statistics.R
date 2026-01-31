# Isabella Muscettola
# 11/25/23
#goal: 1) dissimilartity PERMANOVA on GG2 assigned taxonomy (Bray-Curtis & UniFrac
#distances) using fire history, quadrat median burn severity, horizon, and unit 
#as explanatory variables 2) dissimilarity PERMANOVA on GG2 assigned taxonomy 
# investigating the efficacy of different metrics of burn severity (a) % substrate
# (% litter, % bare, & % rock), b) sampled soil level, c) weighted average of 
#quadrat burn severity, d) canopy mortality and e) remote sense dNBR burn sev

#load libraries
library(patchwork)
library(ggplot2)
library(phyloseq)
library(vegan)
library(dplyr)
library(tidyr)
library(viridis)


#load phyloseq objects
#ps = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.samples")
#ps.3 = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")
ps.3 = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")
hell.ps.3 = transform_sample_counts(ps.3, function(x) (x / sum(x))^0.5 )

#check that ps.3 has 312 samples
check = data.frame(sample_data(ps.3))
rm(check)

#####As of 11/25/2023 the phyloseq objects have been updated to include the new
#metadata that has cleaned up metrics of burn severities.

##################PERMANOVAs###########################################

#11/27 Originally, I was ordinating all of the communities together and 
#faceting by Unit. But since the vegetation and soils differ across the 3 units, 
#a more appropriate visualization would be to have each Unit subset 
######Bray-Curtis Dissimilarities
# Create ordination
#ps.hell = transform_sample_counts(ps.3, function(x) (x / sum(x))^0.5 )
#hell.pcoa.ord = ordinate(ps.hell, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
#x = data.frame(hell.pcoa.ord$vectors)$Axis.1
#y = data.frame(hell.pcoa.ord$vectors)$Axis.2
#df.ord.hell = data.frame(sample_data(ps.hell))
#df.ord.hell$PCoA1 = x
#df.ord.hell$PCoA2 = y


#get values for axes
#Var.Axis.1 = hell.pcoa.ord$values$Relative_eig[1]
#Var.Axis.2 = hell.pcoa.ord$values$Relative_eig[2]
#Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
#Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

#Color by burn severity (weighted average of the quadrat level)
#p = ggplot(df.ord.hell)
#p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Wt_Avg_QuadBS,shape=Horizon))
#p = p + theme_bw()
#p = p + facet_grid(~Unit)
#p = p + xlab(Axis1.label)
#p = p + ylab(Axis2.label) + scale_color_gradient(low ="yellow", high="purple")
#p

#color by treatment
#p = ggplot(df.ord.hell)
#p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon))
#p = p + theme_bw()
#p = p + facet_grid(~Unit)
#p = p + xlab(Axis1.label)
#p = p + ylab(Axis2.label) +
#p


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

#Color by treatment
a = ggplot(df.ord.hell.BL)
a = a + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon), size=2)
a = a + theme_bw() + labs(title="Butte Lake", x=Axis1.label, y=Axis2.label)
a = a + scale_color_manual(values=c('#ffbb22', '#0097d3')) 
a = a + scale_shape_manual(values = c(15,17)) + theme(legend.position = "none")
a

#Color by burn severity (sample severity)
b = ggplot(df.ord.hell.BL)
b = b + geom_point(aes(x=PCoA1,y=PCoA2,color=mtbs_sev,shape=Horizon), size=2)
b = b + theme_bw() + labs(title="Butte Lake", x=Axis1.label, y=Axis2.label)
b = b + scale_color_gradient(low ="yellow", high="purple")
b = b + scale_shape_manual(values=c(15,17)) + theme(legend.position = "none") 
b

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
c = c + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon), size=2)
c = c + theme_bw() + labs(title="Hole", x=Axis1.label, y=Axis2.label)
c = c + scale_color_manual(values=c('#ffbb22', '#0097d3')) 
c = c + scale_shape_manual(values = c(15,17))+ theme(legend.position = "none")
c

d = ggplot(df.ord.hell.H)
d = d + geom_point(aes(x=PCoA1,y=PCoA2,color=mtbs_sev,shape=Horizon), size=2)
d = d + theme_bw() + labs(title="Hole", x=Axis1.label, y=Axis2.label)
d = d + scale_color_gradient(low ="yellow", high="purple") + theme(legend.position = "none")
d = d + scale_shape_manual(values=c(15,17))
d

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

#Color by treatment
e = ggplot(df.ord.hell.WV)
e = e + geom_point(aes(x=PCoA1,y=PCoA2,color=factor(Trt),shape=Horizon), size=2)
e = e + theme_bw() + labs(title="Warner Valley", x=Axis1.label, y=Axis2.label)
e = e + scale_color_manual(values=c('#ffbb22', '#0097d3')) 
e = e + scale_shape_manual(values = c(15,17)) + theme(legend.position="none")
e

f = ggplot(df.ord.hell.WV)
f = f + geom_point(aes(x=PCoA1,y=PCoA2,color=mtbs_sev,shape=Horizon), size=2)
f = f + theme_bw() + labs(title="Warner Valley", x=Axis1.label, y=Axis2.label)
f = f + scale_color_gradient(name="Core Burn Severity", limits =c(1,4),  low ="yellow", high="purple") + theme(legend.position = "none")
f = f + scale_shape_manual(values=c(15,17))
f


a + c + e
b+d+f

?ggsave()
ggsave("Trt_Ordination.jpeg", plot = last_plot(), device = NULL, path = NULL,  scale = 1, width = 10, height = 7.2, units = c("in"),dpi = 600, limitsize = TRUE)

#################STATS############################################
#11/27/23; meeting with Stats consultant I may be overly permuting by treating 
#the replicates within a plot as independent observations. This is true if the
#bacterial #communities do indeed cluster based on plot number, which is not true
#when ordinating. This suggests that treating the replicates within a plot as
#independent may not be unreasonable and in my case is preferential considering
#that the almost all other variables (aside from canopy mortality and other plo
#-level classifications like vegetation type) are collected at the quadrat level
#and thus, I would be losing not only the variability of community composition by 
#pooling that data, but also the variation in all of my other metrics. I will 
#need to do some more digging once I have soil property data (pH, total C/N)
#to better understand whether there is clustering of my soil samples at a plot 
#level. 
#Additionally, since my sampling design was meant to ask the question of does 
#fire history predict changes in community composition within a unit where 
#vegetation and soil type have been controlled, I need to constrain my permutations
#to within a unit only. To do that I am using the "strata" variable within adonis2.

###Test for significance of Treatment where strata=Unit
# Get distance matrix
dist.hell = phyloseq::distance(hell.ps.3, method="bray")

#separated by unit
dist.hell.BL = phyloseq::distance(ps.BL, method="bray")
dist.hell.H = phyloseq::distance(ps.H, method="bray")
dist.hell.WV = phyloseq::distance(ps.WV, method="bray")


# Get sample_data as dataframe
SamDat = data.frame(sample_data(hell.ps.3))
SamDat.BL = data.frame(sample_data(ps.BL))
SamDat.H = data.frame(sample_data(ps.H))
SamDat.WV = data.frame(sample_data(ps.WV))

# test if Trt, horizon, unit, and median_bs are significant
hell.adonis = adonis2(dist.hell ~ Horizon + Trt, strata = SamDat$Unit, SamDat)
hell.adonis # p < 0.001 R2 = 0.048

#is unit or horizon better describe variance in community composition? 
hell.adonis.U = adonis2(dist.hell ~ Unit + Horizon + Trt, SamDat)
hell.adonis.U #Unit R2 =0.08063, Horizon R2=0.08673
hell.adonis.H = adonis2(dist.hell ~ Horizon + Unit + Trt, SamDat)
hell.adonis.H #Horizon R2=0.9380, Unit R2=0.07356

hell.adonis.BL = adonis2(dist.hell.BL ~ Horizon + Trt, SamDat.BL)
hell.adonis.H = adonis2(dist.hell.H ~ Horizon + Trt, SamDat.H)
hell.adonis.WV = adonis2(dist.hell.WV ~ Horizon + Trt, SamDat.WV)
hell.adonis.BL
hell.adonis.H
hell.adonis.WV

hell.adonis.BL = adonis2(dist.hell.BL ~ Horizon + Wt_Avg_QuadBS, SamDat.BL)
hell.adonis.H = adonis2(dist.hell.H ~ Horizon + Wt_Avg_QuadBS, SamDat.H)
hell.adonis.WV = adonis2(dist.hell.WV ~ Horizon + Wt_Avg_QuadBS, SamDat.WV)

#1/10/24; 6 PERMANOVAs to test whether more/fewer fires is significant predictor
#of community composition 



# OLD: when not including strata in the adonis function all 3 explanatory values
#included (Horizon, Unit & Treatment) are significant in explaining the variance 
#of community composition (p<0.001) R2 values vary: Horizon: R2=0.094, 
#Unit: R2 = 0.074, Trt: R2 = 0.043.


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
permutest(b.Trt) # here the results of this permutation test is significant; this means that we cannot reject the null hypothesis that our groups have the same dispersions; we need to be skeptical of the adonis....

#let's dig deeper into what is driving this significant difference in dispersion
b.Trt.BL = betadisper(dist.hell.BL, SamDat.BL$Trt)
b.Trt.H = betadisper(dist.hell.H, SamDat.H$Trt)
b.Trt.WV = betadisper(dist.hell.WV, SamDat.WV$Trt)
permutest(b.Trt.BL) #significant (p < 0.001)
permutest(b.Trt.H) #NOT significant
permutest(b.Trt.WV) #NOT significant

# test if Trt, horizon, unit, and median_bs are significant
hell.adonis.strata = adonis2(dist.hell ~ Horizon + Trt, strata = SamDat$Unit, SamDat)
hell.adonis.strata
# all 3 explanatory values included (Horizon, Unit & Treatment) are significant
# in explaining the variance of community composition (p<0.001) 
#R2 values vary: Horizon: R2=0.094, Unit: R2 = 0.074, Trt: R2 = 0.043.


#############BurnSeverity
ps.bs.drop = transform_sample_counts(subset_samples(ps.3, Wt_Avg_QuadBS > 0), function(x) (x / sum(x))^0.5 )
d = data.frame(sample_data(ps.bs.drop))

dist.hell.drop = phyloseq::distance(ps.bs.drop, method="bray")

adonis.unit = adonis2(dist.hell.drop ~ Horizon + Wt_Avg_QuadBS, data = d, strata=d$Unit)
adonis.unit #R2 0.103, p<0.001

###############Now let's look at different metrics of burn severity############

#############Smallest spatial scale
#######Core burn severity 
ps.hell2 = transform_sample_counts(subset_samples(ps.3, SampleSev>0), function(x) (x / sum(x))^0.5 )

#adonis
dist.hell2 = phyloseq::distance(ps.hell2, method="bray")
SamDat2 = data.frame(sample_data(ps.drop)) #282 observations

hell.adonis2 = adonis2(dist.hell2 ~ Horizon + Unit + Trt + SampleSev, SamDat2, strata = SamDat2$Unit)
hell.adonis2
#The burn severity at the sampled core level is a significant predictor of 
#community dissimilarity (p<0.001) R2 = 0.048.

###THIS NEEDS WORK: subset by unit because BL has unbalanced comparisons
#BL
ps.hell2.BL = subset_samples(ps.hell2, Unit == "BL")
dist.hell2.BL = phyloseq::distance(ps.hell2.BL, method="bray")
SamDat.BL = data.frame(sample_data(ps.hell2.BL))
hell.adonis.BL.SampleSev = adonis2(dist.hell2.BL ~ Horizon +  Trt + SampleSev, SamDat.BL)

hell.adonis.BL.SampleSev

#WV
ps.hell2.WV = subset_samples(ps.hell2, Unit == "WV")
dist.hell2.WV = phyloseq::distance(ps.hell2.WV, method="bray")
SamDat.WV = data.frame(sample_data(ps.hell2.WV))
hell.adonis.WV.SampleSev = adonis2(dist.hell2.WV ~ Horizon + Trt + SampleSev, SamDat.WV)
#Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : contrasts can be applied only to factors with 2 or more levels
hell.adonis.WV.SampleSev

#H
ps.hell2.H = subset_samples(ps.hell2, Unit == "H")
dist.hell2.H = phyloseq::distance(ps.hell2.H, method="bray")
SamDat.H = data.frame(sample_data(ps.hell2.H))
hell.adonis.H.SampleSev = adonis2(dist.hell2.H ~ Horizon + Trt + SampleSev, SamDat.H)
#Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : contrasts can be applied only to factors with 2 or more levels
hell.adonis.H.SampleSev




#####Quadrat weighted average of burn severity 
#######Core burn severity 
# 1 Quadrat did not have burn severity assigned (WV 01 C)
ps.drop = subset_samples(ps.3, Wt_Avg_QuadBS>0)
ps.hell2 = transform_sample_counts(ps.drop, function(x) (x / sum(x))^0.5 )

#adonis
dist.hell2 = phyloseq::distance(ps.hell2, method="bray")
SamDat2 = data.frame(sample_data(ps.drop))

hell.adonis3.1 = adonis2(dist.hell2 ~ Horizon + Trt, strata = SamDat2$Unit, SamDat2)
hell.adonis3.1 #p<0.001 but R2 4.8%

hell.adonis3 = adonis2(dist.hell2 ~ Horizon + Wt_Avg_QuadBS + Trt, strata = SamDat2$Unit, SamDat2)
hell.adonis3 # BS R2 = 10.3%; Trt R2 = 1.3% p < 0.001
#The average burn severity at the quadrat level is a significant predictor of 
#community dissimilarity (p<0.001) R2 = 0.042. 
#may not be fair to compare R2 of the 2 models bc size of dataframes are different

#Is this true for all 3 units? 
dist.hell2.BL = phyloseq::distance(subset_samples(ps.hell2, Unit =="BL"), method="bray")
hell.ad.BL2 = adonis2(dist.hell2.BL ~ Horizon + Wt_Avg_QuadBS + Trt, strata = SamDat.BL$Unit, SamDat.BL)
hell.ad.BL2 #R2 = 2.9%, p < 0.001

dist.hell2.H = phyloseq::distance(subset_samples(ps.hell2, Unit =="H"), method="bray")
hell.ad.H2 = adonis2(dist.hell2.H ~ Horizon + Wt_Avg_QuadBS + Trt, strata = SamDat.H$Unit, SamDat.H)
hell.ad.H2 # p < 0.05 R2 = 2.6%

dist.hell2.WV = phyloseq::distance(subset_samples(ps.hell2, Unit =="WV"), method="bray")
ps.drop.hell.WV = subset_samples(ps.hell2, Unit=="WV")
SamDat2.WV = data.frame(sample_data(ps.drop.hell.WV))
hell.ad.WV2 = adonis2(dist.hell2.WV ~ Horizon + Wt_Avg_QuadBS + Trt, SamDat2.WV)
hell.ad.WV2 # p<0.01 R2 = 2.9%



#####Substrate
##Litter + (Bare + Rock)
# NA Substrate in 1 Quadrat: BL PIJE 06 A
ps.hell.drop = subset_samples(ps.hell, LITT!="NA")
dist.hell2 = phyloseq::distance(ps.hell.drop, method="bray")
SamDat2 = data.frame(sample_data(ps.hell.drop))

hell.adonis4 = adonis2(dist.hell2 ~ Horizon + Trt + LITT + BARE_ROCK, strata = SamDat2$Unit, SamDat2)
hell.adonis4

## Switch order of LITT and BARE_ROCK
hell.adonis5 = adonis2(dist.hell2 ~ Horizon + Trt + BARE_ROCK + LITT, strata = SamDat2$Unit, SamDat2)
hell.adonis5


########Canopy mortality 
hell.adonis6 = adonis2(dist.hell ~ Horizon + Trt + Overstory_mort, strata = check$Unit, check)
hell.adonis6

############ FULL PERMANOVA ####################################

hell.adonis3 = adonis2(dist.hell2 ~ Horizon + Wt_Avg_QuadBS + Trt, strata = SamDat2$Unit, SamDat2)
hell.adonis3
