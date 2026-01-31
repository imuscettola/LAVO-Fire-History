#11/22/23 
#objectives:1)  to visualize the OTU relative abundances for each sample with 
#taxonomy assigned by the Greengenes2 database. 2) to compare the taxonomy assigned
#by the Silva v. Greengenes2 database. 3) Evaluate which taxonomy assigning classifier
#is best for subsequent data analyses. 

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

#import phyloseq object
ps = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.samples")
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5)
df = sample_data(ps)

ps.raw = readRDS("~/Box/MSWhitman/LAVO_FireHistory/Data/Sequencing/Greengenes2_Taxonomy/ps.raw")
ps.hell.raw = transform_sample_counts(ps.raw, function(x) (x / sum(x))^0.5)

##############Creating dataframes; only do once################################
#thanks to code from: https://github.com/joey711/phyloseq/issues/1521
#############Create Phylum dataframe
phylum.df = ps.hell.raw %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt() %>%
  select(Phylum, Sample, Abundance) %>%
  filter(Abundance > 0) %>%
  mutate(SampleID = Sample) %>%
  mutate(Hellinger_relAbundance = Abundance) %>%
  select(SampleID, Phylum, Hellinger_relAbundance)

#pull metadata
df.full = sample_data(ps.hell.raw)
#merge phylum table with metadata
phylum.df = left_join(phylum.df, df.full, by ="SampleID")
phylum.df = select(phylum.df, 1:3, 9:26, 28:31)
#save csv
write.csv(phylum.df, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_Phylum.csv")


#######Create OTU table 
write.csv(ps.hell.raw %>% 
            psmelt() %>%
            arrange(OTU) %>%
            select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
            spread(Sample, Abundance),
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_HellingerAll.csv")

write.csv(ps.raw %>% 
            psmelt() %>%
            arrange(OTU) %>%
            select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
            spread(Sample, Abundance),
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_All.csv")

write.csv(ps.raw %>% 
            psmelt() %>%
            arrange(OTU) %>%
            select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance), 
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_All_long.csv")

###########Create Genus dataframe
genus.df = ps.hell.raw %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt() %>%
  select(Genus, Sample, Abundance) %>%
  filter(Abundance > 0) %>%
  mutate(SampleID = Sample) %>%
  mutate(Hellinger_relAbundance = Abundance) %>%
  select(SampleID, Genus, Hellinger_relAbundance)

#pull metadata
df.full = sample_data(ps.hell.raw)
#merge phylum table with metadata
genus.df = left_join(genus.df, df.full, by ="SampleID")
genus.df = select(genus.df, 1:3, 9:26, 28:31)
#save csv
write.csv(genus.df, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_Genus.csv")

#cleanup
rm(df.full)
rm(genus.df.long)
rm(phylum.df)




################### Start from here #####################################

ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5) 
phylum.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_Phylum.csv")
genus.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_Genus.csv")
raw.otu.table = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/GG2/Raw_OTU_All_long.csv")

#######################Summarize OTU Table###################################
#Goal is to compare the taxonomy assigning using either Silva or Greengenes database
raw.genus.summary = raw.otu.table %>%
  group_by(Genus) %>%
  summarise(sum(Abundance)) #number of NA Genuses 3491430 # of empty values: 413549
raw.species.summary = raw.otu.table %>%
  group_by(Species) %>%
  summarise(sum(Abundance)) # NA: 7381751, empty: 1080563

unique.genus = raw.otu.table %>%
  count(Genus)



########################## OTU TABLES ######################################
# transform total counts to relative abundance using the hellinger transformation
otu.hell = otu_table(ps.hell)
tibble(otu.hell)

#plot bar graph of OTUs by Unit
plot_bar(ps.hell, fill = "Phylum") + facet_wrap(~Unit, scales = "free_x", nrow =2)
#takes SO LONG moving forward might be worthwhile to subset the phyloseq object
#by Unit
# I want to make dataframes so that I can plot using ggplot instead of using
# plot_bar in phyloseq...


##################OTU Phylum Barplots###########################################

#######Graph Phylum for BL Unit facet by plot
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="BL" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="BL" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


############## Warner Valley
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="WV" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="WV" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p

############# Hole
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="H" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="H" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


############Unburned
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="UB" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="UB" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


########################Ordinations colored by Plot##########################
#########NOTE: 
## This produces the same results as the Silva_exploring_OTUs code bc it pulls 
#from OTU matrix, does not use any assigned taxonomy in the ordinations. 
hell.pcoa.ord = ordinate(ps.hell, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(hell.pcoa.ord$vectors)$Axis.1
y = data.frame(hell.pcoa.ord$vectors)$Axis.2
df.ord.hell = data.frame(sample_data(ps.hell))
df.ord.hell$PCoA1 = x
df.ord.hell$PCoA2 = y


#get values for axes
Var.Axis.1 = hell.pcoa.ord$values$Relative_eig[1]
Var.Axis.2 = hell.pcoa.ord$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")


###############Butte Lake
#BL Organic LF 
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "BL" & df.ord.hell$Trt=="LF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

#BL Organic MF 
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "BL" & df.ord.hell$Trt=="MF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# BL Mineral LF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "BL" & df.ord.hell$Trt=="LF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

#BL Mineral MF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "BL" & df.ord.hell$Trt=="MF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

############Warner Valley 
#WV Organic LF 
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "WV" & df.ord.hell$Trt=="LF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

#WV Organic MF 
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "WV" & df.ord.hell$Trt=="MF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# WV Mineral LF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "WV" & df.ord.hell$Trt=="LF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

#WV Mineral MF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "WV" & df.ord.hell$Trt=="MF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p


###########Hole
# H Organic
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "H", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# H Mineral LF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "H" & df.ord.hell$Trt=="LF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# H Mineral MF
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "H" & df.ord.hell$Trt=="MF", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

###########Unburned
# UB Organic
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="O" & df.ord.hell$Unit == "UB", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

# UB Mineral
p = ggplot(df.ord.hell[df.ord.hell$Horizon=="M" & df.ord.hell$Unit == "UB", ])
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Plot_Name))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p



##################Now lets look at the Resequenced Samples#####################
resequenced = phylum.df %>%
  filter(Resequence =="Y")
resequenced = mutate(resequenced, BSA = ifelse(endsWith(resequenced$SampleID, "2"), "correct", "concentrated"))



#####################OTU tables of resequenced samples########################
p = ggplot() + theme_bw()
p = p + geom_col(data=resequenced, aes(x=BSA, y=Hellinger_relAbundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Full_Sample_Name))
p

##################Blanks########################################################
#blanks alone
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Blank =="Y", ], aes(x=SampleID, y=Hellinger_relAbundance, fill=Phylum), position = "fill") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p

#blanks with samples from the same DNA Extract Round
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Extract_Rd == "1", ], aes(x=SampleID, y=Hellinger_relAbundance, fill=Phylum), position = "fill") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p

#let's visualize the blanks in an ordination to more easily judge whether they are
#different from real samples
#NOTE: These ordinations are the same as the Silva assigned taxonomy ps objects
blanks.ord = ordinate(ps.hell, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(blanks.ord$vectors)$Axis.1
y = data.frame(blanks.ord$vectors)$Axis.2
df.ord.hell = data.frame(sample_data(ps.hell))
df.ord.hell$PCoA1 = x
df.ord.hell$PCoA2 = y


#get values for axes
Var.Axis.1 = blanks.ord$values$Relative_eig[1]
Var.Axis.2 = blanks.ord$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")


df.ord.hell$Extract_Rd = as.character(df.ord.hell$Extract_Rd)
p = ggplot(df.ord.hell)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Blank, shape=Extract_Rd))
p = p + theme_bw()
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p