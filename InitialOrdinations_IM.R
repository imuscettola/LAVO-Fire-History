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


# Import ps object

ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
ps


# Get relative abundances
ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5 )
otu_table(ps.norm)[1:6,1:6]

# Create an ordination
ord.norm = ordinate(ps.norm,method="PCoA",distance="bray")
ord.hell = ordinate(ps.hell,method="PCoA",distance="bray")

# Plot the ordination
plot_ordination(physeq=ps.norm,ordination=ord.norm,color="Trt",shape="Unit")

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

# Create a UniFrac ordination of hellinger transformed relative abundances

