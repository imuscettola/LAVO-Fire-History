#updating metadata
library(dplyr)
library(tidyverse)
library(ggplot2)

#03112024 add pH to metadata; add the NEG & POS controls back to the metadata
oldmetadata <- read.csv("~/Desktop/LAVO22_GG2/03112024_LAVO22_Metadata.csv")
oldoldmetadata <- read.csv("~/Desktop/LAVO22_GG2/02152024_LAVO22_Metadata.csv")
POSandNEG <- filter(oldoldmetadata, Unit=="NEG" | Unit=="POS") %>%
  mutate()

m <- colnames(oldmetadata)
o <- colnames(POSandNEG)
t <- c(m, o)
!unique(t)
unique(t)

newmetadata<-rbind(oldmetadata, POSandNEG)

#####################NEED TO WORK ON THIS.... ###############################

#03112024 add pH to metadata NOTE THIS METADATA IS MISSING NEG and POS controls! 
oldmetadata <- read.csv("~/Desktop/LAVO22_GG2/03072024_LAVO22_Metadata.csv")
pH <- read.csv("~/Desktop/LAVO22_GG2/pH/LAVO22_pH.csv")
pH <- select(pH, Sample, pH) %>%
  filter(., pH!="NA")
temp <- left_join(oldmetadata, pH, by=c("Full_Sample_Name" = "Sample"))

write.csv(temp, "~/Desktop/LAVO22_GG2/03112024_LAVO22_Metadata.csv")

# 03072024 metadata
oldmetadata <- read.csv("~/Desktop/LAVO22_GG2/03072024_LAVO22_Metadata.csv")
soilkey <- read.csv("~/Desktop/LAVO22_GG2/LAVO_SoilKey.csv")
sievemass <- read.csv("~/Desktop/LAVO22_GG2/LAVO22_SievedMass.csv")
sievemass <- sievemass %>%
  select(SampleID, Perc_Fine)

temp <- left_join(oldmetadata, soilkey, by="Soil_mu")
temp2 <- left_join(temp, sievemass, by=c("Full_Sample_Name" = "SampleID"))
write.csv(temp2, "~/Desktop/LAVO22_GG2/03072024_LAVO22_Metadata.csv")

# 02222024 metadata
mtbs <- read.csv("~/Documents/Mapping/mtbs/3_final/plots_mtbs.csv")
mtbs <- mtbs %>%
  select(1, 14:16) #mtbs_sev
oldmetadata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/OldMetadata/02212024_LAVO22_Metadata.csv")
oldmetadata <- select(oldmetadata, -35)

newmetadata <- merge(oldmetadata, mtbs, by.x="Plot_Name", by.y="Plot", all=TRUE) %>%
  mutate(., mtbs_sev=mtbs) %>%
  select(., -mtbs)

write.csv(newmetadata, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/02222024_LAVO22_Metadata.csv")

## Greengenes2 database - Making phyloseq object  
#code was compiled from the "0_16S_Seq_processing.Rmd" file created on 4/24/2023


library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)


# Import .biom file 
#ps.raw = import_biom("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/full-gg_OTU_table/feature-table-metaD-tax_json.biom" , parseFunction = parse_taxonomy_default)
ps.raw = import_biom("~/Desktop/LAVO22_GG2/full-gg_OTU_table/feature-table-metaD-tax_json.biom" , parseFunction = parse_taxonomy_default)
ps.raw
#

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

rm(x)

#### Shape up the full PS object, merge sequencing runs  ####

# We have duplicate sequencing for each sample. We need to merge each sample's duplicates.
# Use merge function in phylogeny to add the sample counts together (for this dataset, this is what we want.)
# You’ll need a column in the sample data that has the ID the same for all the paired samples.
# here it is "Submission_Sample_Name" (previously SampleID)

length(sample_data(ps.raw)$Submission_Sample_Name) / length(levels(as.factor(sample_data(ps.raw)$Submission_Sample_Name)))

ps = merge_samples(ps.raw, sample_data(ps.raw)$Submission_Sample_Name)
warnings()
?merge_samples
ps #The taxonomy is there...
samdat = sample_data(ps) #...But the sample metadata is not.

sample_data(ps)$Submission_Sample_Name

#rename row names
#row.names(sample_data(ps)) = sample_data(ps)$Submission_Sample_Name ##error: non-unique values when setting row.names; duplicate 'row.names' are not allowed.... could be because I have some reruns of samples. try merging on numeric ascending list "MySampleID"
#sample_data(ps)
#sample_data(ps)$Submission_Sample_Name

samdatdftemp = data.frame(SampleID = sample_names(ps),Submission_Sample_Name = sample_names(ps))
row.names(samdatdftemp) = samdatdftemp$SampleID
samdatdftemp = sample_data(samdatdftemp)
sample_data(ps) = samdatdftemp


#metadata = read.csv(file ="~/Box/MSWhitman/LAVO_FireHistory/Data/R/03072024_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)
metadata = read.csv(file ="~/Desktop/LAVO22_GG2/03112024_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)
head(metadata)

metadata$Run = as.factor(substring(metadata$Run_Submission_Name,1,1))
metadata = metadata %>%
  filter(Run=="1")

#merge metadata df with sample_data to put the metadata in the correct order
df2 = merge(data.frame(sample_data(ps)),metadata,by="Submission_Sample_Name")
head(df2)

row.names(df2) = df2$SampleID
sample_data(ps) = sample_data(df2)



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

# Save ps object ps.raw
#saveRDS(ps, "~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.raw")
#saveRDS(ps, "~/Desktop/LAVO22_GG2/phyloseqobjects/ps.raw")
#ps = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.raw")
ps = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.raw")
ps

# Remove blanks, positive control, and duplicates
ps.samples = subset_samples(ps, Pos_Cntrl == "N")
ps.samples = subset_samples(ps.samples, Blank == "N")
ps.samples = subset_samples(ps.samples, !endsWith(Submission_Sample_Name, "2"))
ps.samples

# Save ps object ps.samples
#saveRDS(ps.samples,"~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.samples")
#saveRDS(ps.samples,"~/Desktop/LAVO22_GG2/phyloseqobjects/ps.samples")
#ps.samples = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.samples")
ps.samples = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.samples")

# Remove Unburn Unit so just the Burned Units Remain
ps.3 = subset_samples(ps.samples,Unit != "UB")
#add Plot_horizon to sample_data
sample_data(ps.3)$Plot_Horizon = as.factor(paste(sample_data(ps.3)$Plot_Name, sample_data(ps.3)$Horizon, sep="_"))
ps.3

# Save ps object ps.3units
#saveRDS(ps.3, "~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3units")
#saveRDS(ps.3, "~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")

ps.3 = readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")
ps.3

ps.3.SamDat = as.data.frame(sample_data(ps.3))
ps.3.SamDat

# save df.3 as csv: ERROR : invalid class "sample_data" object: Sample Data must have non-zero dimensions.
#write.csv(ps.3.SamDat, file ="~/Desktop/LAVO22_GG2/phyloseqobjects/ps3SamDat.csv")

######### Average the relative abundances of the OTUs within a plot  #########
#note: sample relative abundance = (sample otu count) / (sample total count)
#      plot relative abundance = average(rep A rel. abund, rep B rel. abund, rep C rel. abund)

ps.3.norm <- transform_sample_counts(ps.3, function(x) x / sum(x))
ps.3.norm.avg <- merge_samples(ps.3.norm, group='Plot_Horizon', fun=mean)

#check otu table output
otu_table(ps.3.norm.avg)[1:3, 1:5] #how is it possible that there are relative abundances > 1.0?

#fix sample_data
samdatdftemp = data.frame(Plot_Horizon = sample_names(ps.3.norm.avg))
row.names(samdatdftemp) = samdatdftemp$Plot_Horizon
samdatdftemp = sample_data(samdatdftemp)
sample_data(ps.3.norm.avg) = samdatdftemp

#metadata_raw = read.csv(file ="~/Box/MSWhitman/LAVO_FireHistory/Data/R/03072024_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)
metadata_raw = read.csv(file ="~/Desktop/LAVO22_GG2/03112024_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)
metadata = metadata_raw %>%
  mutate(Plot_Horizon=paste(Plot_Name, Horizon, sep="_")) %>%
  select(1, 9:11, 17:21, 35, 36) %>%
  distinct(.keep_all = TRUE) %>%
  filter(Unit!="UB") %>%
  mutate(Plot_Horizon=paste(Plot_Name, Horizon, sep="_"))
head(metadata)

#merge metadata df with sample_data to put the metadata in the correct order
df2 = merge(data.frame(sample_data(ps.3.norm.avg)),metadata,by="Plot_Horizon")
head(df2)

row.names(df2) = df2$Plot_Horizon
sample_data(ps.3.norm.avg) = sample_data(df2)

rm(metadata_raw, df2, samdatdftemp)

#Save ps object 
#saveRDS(ps.3.norm.avg, "~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3.norm.avg")
#saveRDS(ps.3.norm.avg, "~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3.norm.avg")

##Avg relative abundances by merging otu counts within a plot
ps.3.avg <- merge_samples(ps.3, group='Plot_Horizon', fun=mean)

otu_table(ps.3.avg)[1:3, 1:10]

#fix metadata
samdatdftemp = data.frame(Plot_Horizon = sample_names(ps.3.avg))
row.names(samdatdftemp) = samdatdftemp$Plot_Horizon
samdatdftemp = sample_data(samdatdftemp)
sample_data(ps.3.avg) = samdatdftemp

df2 = merge(data.frame(sample_data(ps.3.avg)),metadata,by="Plot_Horizon")
head(df2)

row.names(df2) = df2$Plot_Horizon
sample_data(ps.3.avg) = sample_data(df2)

temp <- as.data.frame(sample_data(ps.3.avg))

#save ps object
#saveRDS(ps.3.avg, "~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.3.avg")
saveRDS(ps.3.avg, "~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3.avg")


