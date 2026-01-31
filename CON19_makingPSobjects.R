library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)

setwd("~/Desktop/Community Assembly/2019")

# Import .biom file
ps.raw = import_biom("biom_files/feature-table-metaD-tax_json_CON19-FULL.biom" , parseFunction = parse_taxonomy_default)
ps.raw
sample_data(ps.raw)

#### Check taxonomy ####
head(tax_table(ps.raw))
# Yup, looks ugly.

##Fixing the tax_table object to "clean-up" the column names

x = data.frame(tax_table(ps.raw))
# Making a dummy variable to store the taxonomy data

colnames(x) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# Assigning the proper column names instead of SILVA ranks

x$Domain = gsub("d__", "", as.character(x$Domain))
x$Phylum = gsub("p__", "", as.character(x$Phylum))
x$Class = gsub("c__", "", as.character(x$Class))
x$Order = gsub("o__", "", as.character(x$Order))
x$Family = gsub("f__", "", as.character(x$Family))
x$Genus = gsub("g__", "", as.character(x$Genus))
x$Species = gsub("s__", "", as.character(x$Species))
# Substituting the characters we don't want with nothing in the taxonomy

x=tax_table(as.matrix(x,dimnames=list(row.names(x),colnames(x))))
# Turning it into a taxonomy table, while saving the rownames and column names
tax_table(ps.raw)=x
# Reassigning the taxonomy table in ps_xxx to the new modified one

head(tax_table(ps.raw))
# Check for success
# Looks good.

#### Shape up the full PS object, merge sequencing runs  ####

# We have duplicate sequencing for each sample. We need to merge each sample's duplicates.
# Use merge function in phylogeny to add the sample counts together (for this dataset, this is what we want.)
# You’ll need a column in the sample data that has the ID the same for all the paired samples.
# here it is "MySampleID"

ps = merge_samples(ps.raw, sample_data(ps.raw)$MySampleID)
warnings()
ps #The taxonomy is there...
sample_data(ps) #...But the sample metadata is not.

#rename row names
row.names(sample_data(ps)) = sample_data(ps)$MySampleID
sample_data(ps)

metadata = read.csv(file ="sample_metadata_connors_2019.csv", header=TRUE, stringsAsFactors=TRUE)
metadata
#Fix SampleID column name in ps object
sample_data(ps)$SampleID = sample_data(ps)$MySampleID
sample_data(ps)

#merge metadata df with sample_data to put the metadata in the correct order
df2 = merge(data.frame(sample_data(ps)),metadata,by="SampleID")
head(df2)
colnames(df2)
df3 = df2[,c(1, 19:35)]
colnames(df3)
colnames(df3) = c("SampleID", "Rep", "Blank", "TimesMixed", "ExtrBlank", "IncubBlank", "Initial",
                  "WaterBlank", "VortexControl", "Trt", "GroupOrAgit", "Group2",
                  "O2", "Moisture", "Meta", "MetaOrder", "MixingStudy", "ConnectivityStudy")
colnames(df3)
head(df3)

row.names(df3) = df3$SampleID
sample_data(ps) = sample_data(df3)
str(sample_data(ps))
sample_data(ps)
 
# get the stragglers into factors
sample_data(ps)$Rep = factor(sample_data(ps)$Rep)
sample_data(ps)$MetaOrder = factor(sample_data(ps)$MetaOrder)
# check that the correct number of levels for each factor, ie, that NAs are not factor levels
str(sample_data(ps)) 
sample_data(ps)

# Remove chloroplast and mitochondria from the dataset
ps
ps <- ps %>%
  subset_taxa(
    Family != "Mitochondria" &
      Class != "Chloroplast" &
      Order != "Chloroplast" &
      Domain != "Eukaryota"
  )
ps = prune_taxa(taxa_sums(ps) > 0, ps)
ps

#### make a ps object for just the mixed community "generous donor" sample ####
ps.mixed = subset_samples(ps, SampleID=="499")
nsamples(ps.mixed) #should be 1
sample_data(ps.mixed)
ps.mixed = prune_taxa(taxa_sums(ps.mixed) > 0, ps.mixed)
ps.mixed
saveRDS(ps.mixed,"ps.CON19.generousdonor")

#### make a ps object including all blanks ####
#filter out the mixed community sample, and resequenced samples
nsamples(ps)
ps1 = subset_samples(ps, SampleID!="499")
ps1 = subset_samples(ps1, SampleID!="3087")
ps1 = subset_samples(ps1, SampleID!="3093")
nsamples(ps1)
ps1 = prune_taxa(taxa_sums(ps1) > 0, ps1)
ps1
saveRDS(ps1,"ps.CON19.withallblanks")

#### cleaning up ps for the real work. ####

# filter out all blanks
ps2 = subset_samples(ps1, Blank!="Y")
nsamples(ps2) #should be 446 (12 trts x 32 + 32 initial + 32 agitation....minus 2 that were dropped during extraction)
#You'll save it a little later..

# remove columns that designate the blanks
sample_data(ps2)
sample_data(ps2) = sample_data(ps2)[,-c(3,5,6,8)]

# Put factors in order
head(sample_data(ps2))
sample_data(ps2)$TimesMixed <- ordered(sample_data(ps2)$TimesMixed, levels=c("Initial", "1X", "4X", "16X"))
sample_data(ps2)$Trt <- ordered(sample_data(ps2)$Trt, levels=c("Initial", "1X_O_U","1X_O_S","1X_A_U","1X_A_S",
                                                               "4X_O","4XV_O","4X_A","4XV_A",
                                                               "16X_O","16X_A","16XV_O","16XV_A",
                                                               "Meta_O_U", "Meta_O_S", "Meta_A_U", "Meta_A_S"))
sample_data(ps2)$GroupOrAgit <- ordered(sample_data(ps2)$GroupOrAgit, levels=c("Initial", "1", "4", "16", "Vortex"))
sample_data(ps2)$O2 <- ordered(sample_data(ps2)$O2, levels=c("Initial", "Oxic", "Anoxic"))
sample_data(ps2)$Moisture <- ordered(sample_data(ps2)$Moisture, levels=c("Initial", "Unsaturated", "Saturated"))
sample_data(ps2)$Meta <- ordered(sample_data(ps2)$Meta, levels=c("Initial", "N", "Meta"))


str(sample_data(ps2))

ps2
ps2 = prune_taxa(taxa_sums(ps2) > 0, ps2)
ps2

saveRDS(ps2,"ps.CON19")


##### Get observed OTUs for the full dataset (including blanks) for metadata.
#estimate_richness(ps1, measures="Observed")