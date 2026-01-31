# goal: look at the Resequenced samples to see if BSA concentration during the 
# amplification process of sample prep led to noticable differences in community
# composition of the samples

library(phyloseq)
library(vegan)

#import phyloseq object
ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")

# 1. subset ps to just resequenced samples 
ps.dup = subset_samples(ps, Resequence=="Y")
# check output
with(sample_data(ps.dup), table(Resequence))
# good, we have 12 samples in this phyloseq object all of them are resequenced

# 2. transform otu table to hellinger transformed relative abundance
ps.hell.dup = transform_sample_counts(ps.dup, function(x) (x / sum(x))^0.5 )
otu_table(ps.hell.dup)[1:6,1:6]

# 3. ordinate 
ord.hell.dup = ordinate(ps.hell.dup,method="PCoA",distance="bray")

# 4. plot ordination to visually inspect whether pairwise samples are different 
# from each other
plot_ordination(physeq=ps.hell.dup,ordination=ord.hell.dup,color="Full_Sample_Name")


# key findings: 
# here it appears that the unburned samples have more spread between the pairs in 
# this ordination than the other samples (BL and WV). This indicates that 
# differences in community composition observed between samples that were 
# amplified using dilute BSA versus the correct concentration of BSA may not be 
# attributable to BSA concentration. In other words, even though BSA concentration 
# appears to have changed the amount of total sequences observed, the communities 
# are relatively equal.