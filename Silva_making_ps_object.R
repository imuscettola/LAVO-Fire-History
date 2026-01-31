## Silva database - Making phyloseq object  
#code was compiled from the "0_16S_Seq_processing.Rmd" file created on 4/24/2023


library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)


# Import .biom file 
ps.raw = import_biom("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/full-silva_OTU_table/feature-table-metaD-tax_json.biom" , parseFunction = parse_taxonomy_default)
ps.raw

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


metadata = read.csv(file ="~/Box/MSWhitman/LAVO_FireHistory/Data/R/11222023_LAVO22_Metadata.csv", header=TRUE, stringsAsFactors=TRUE)
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
#saveRDS(ps, "~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.raw")
#ps

# Remove blanks and positive control
ps.samples = subset_samples(ps, Pos_Cntrl == "N")
ps.samples = subset_samples(ps.samples, Blank == "N")
ps.samples = subset_samples(ps.samples, !endsWith(Submission_Sample_Name, "2"))
ps.samples

# Save ps object
#saveRDS(ps.samples,"~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
#ps.samples = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
#ps.samples


# Remove Unburn Unit so just the Burned Units Remain
ps.3 = subset_samples(ps.samples,Unit != "UB")
ps.3

# Save ps object ps.3units
#saveRDS(ps.3, "~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3units")
#ps.3 = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3units")
#ps.3

# why are there 312 samples remaining here while there were 322 samples remaining
# from the Greengenes2 taxonomy?
df.3 = sample_data(ps.3)
df.3

# save df.3 as csv
#write.csv(df.3, file ="~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3.csv")

#comparing the df.3 of the silva and gg2 databases, despite gg2 output suggesting
#that it has 322 samples, the sample_data table shows only 312 samples. Not 
#sure what's going on here, but will proceed regardless.

