#installing BiocManager and phyloseq
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("phyloseq")


#interaction term: shows effect of treatment is different depending on the unit you are in #other option; subset data for each unit
#####BL: Horizon and Trt are significant (p<0.001), R2 horizon=0.097, R2 trt=0.14


# Make a very fast ordination plot

library(phyloseq)
library(tidyr)
library(ggplot2)
library(vegan)
library(dplyr)


# Import ps object

ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
ps

ps = subset_samples(ps,Unit != "UB")
ps

# Get relative abundances
ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5 )
otu_table(ps.hell)[1:6,1:6]


# Get relative abundances
ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5 )

# Create an ordination
ord.norm = ordinate(ps.norm,method="PCoA",distance="bray")
ord.hell = ordinate(ps.hell,method="PCoA",distance="bray")

# Collect the percent var explained by axes 1 and 2
Var.Axis.1 = ord.norm$values$Relative_eig[1]
Var.Axis.2 = ord.norm$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

# Plot the ordination
plot_ordination(physeq=ps.norm,ordination=ord.norm,color="Trt",shape="Unit")
plot_ordination(physeq=ps.hell,ordination=ord.hell,color="Trt",shape="Unit")


# Plotting ordination without phyloseq
x = data.frame(ord.norm$vectors)$Axis.1
y = data.frame(ord.norm$vectors)$Axis.2
df.ord.norm = data.frame(sample_data(ps.norm))
df.ord.norm$PCoA1 = x
df.ord.norm$PCoA2 = y

# Plotting ordination without phyloseq
x = data.frame(ord.hell$vectors)$Axis.1
y = data.frame(ord.hell$vectors)$Axis.2
df.ord.hell = data.frame(sample_data(ps.hell))
df.ord.hell$PCoA1 = x
df.ord.hell$PCoA2 = y


#saveRDS(df.ord.norm,"df.ord.norm")
# Bring in the ordination data
df.ord.norm = readRDS("df.ord.norm")

x = data.frame(ord.hell$vectors)$Axis.1
y = data.frame(ord.hell$vectors)$Axis.2
df.ord.hell = data.frame(sample_data(ps.hell))
df.ord.hell$PCoA1 = x
df.ord.hell$PCoA2 = y


# Plot faceted ordination for relabund data
p = ggplot(df.ord.norm)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# Let's do some stats

# PERMANOVA

# Get distance matrix
dist.norm = phyloseq::distance(ps.norm,method="bray")
dist.hell = phyloseq::distance(ps.hell,method="bray")

# Save dist.norm
#saveRDS(dist.norm,"dist.norm")
dist.norm = readRDS("dist.norm")

# Test if Trt is significant
norm.adonis = adonis2(dist.norm ~ Trt,df.ord.norm)
norm.adonis

# Test trt and horizon
norm.adonis = adonis2(dist.norm ~ Horizon*Trt,df.ord.norm)
norm.adonis
# BL: Horizon 0.097, Trt 0.149, p<0.001 for both > Bonferroni corr Trt p<0.003
# H: Horizon 0.084, p<0.001; Trt 0.034, p=0.006 > Bonferroni corr Trt p<0.018
# WV: Horiozn 0.16, p<0.001; Trt 0.02, p=0.015 > Bonferroni corr Trt p<0.045

##  more parameters
norm.adonis = adonis2(dist.norm ~ Unit + Horizon + Trt, df.ord.norm)
norm.adonis
# Unit, Horizon and Burn History were all significant predictors of bacterial community
# composition (PERMANOVA, p<0.001,R2unit=0.08)

df.ord.hell.nona <- filter(df.ord.hell, Slope!="NA")

## most parameters
hell.adonis = adonis2(dist.hell ~ Unit + Plot_Name + Horizon + Trt + Rep + Slope, df.ord.hell.nona)
hell.adonis
# Unit, Horizon and Burn History were all significant predictors of bacterial community
# composition (PERMANOVA, p<0.001,R2unit=0.08)

## interactions
norm.adonis = adonis2(dist.norm ~ Unit * Trt, df.ord.norm)
norm.adonis

###################

# Plot faceted ordination for Hellinger-transformed relative abundance data
p = ggplot(df.ord.hell)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon))
p = p + theme_bw()
p = p + facet_wrap(~Unit)
p


##########################
# Try plotting with BSev
df.ord.norm$Core_burnsev.factor = factor(df.ord.norm$Core_burnsev,order=TRUE,levels=c("1","2","3","4"))
p = ggplot(df.ord.norm)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Core_burnsev.factor,shape=Trt))
p = p + theme_bw()
p = p + facet_wrap(~Unit)
p



############################

# Create an ordination
ord.norm = ordinate(ps.norm,method="PCoA",distance="bray")
ord.hell = ordinate(ps.hell,method="PCoA",distance="bray")

# Plot the ordination
plot_ordination(physeq=ps.hell,ordination=ord.hell,color="Trt",shape="Unit")

# Plotting ordination without phyloseq
x = data.frame(ord.norm$vectors)$Axis.1
y = data.frame(ord.norm$vectors)$Axis.2

df.ord.norm = data.frame(sample_data(ps.norm))

df.ord.norm$x = x
df.ord.norm$y = y

p = ggplot(df.ord.norm)
p = p + geom_point(aes(x=x,y=y,color=Trt,shape=Unit))
p = p + theme_bw()
p