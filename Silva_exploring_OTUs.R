#11/22/23 
#objective: to visualize the OTU relative abundances for each sample with taxonomy assigned by the Silva database. 
#long term goal: to compare Silva database taxonomy classifications to Greengenes2. 

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

#import phyloseq object
ps = readRDS("C:/Users/iemus/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.initial.samples")
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5)
df = sample_data(ps)

ps.raw = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.raw")
ps.hell.raw = transform_sample_counts(ps.raw, function(x) (x / sum(x))^0.5)

##############Creating dataframes; only do once################################
#thanks to code from: https://github.com/joey711/phyloseq/issues/1521
#############Create Phylum dataframe
write.csv(ps.hell %>%
            tax_glom(taxrank = "Phylum") %>% 
            psmelt() %>%
            select(Phylum, Sample, Abundance) %>% 
            spread(Sample, Abundance), 
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv")


#read the newly created dataframe
phylum.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv", header = TRUE, sep = ",")
phylum.df = phylum.df[,-1]

df <- data.frame(t(phylum.df[-1]))
colnames(df) <- phylum.df[, 1]
df$SampleID = row.names(df)
df = pivot_longer(df, cols=1:41, names_to="Phylum", values_to = "Hellinger Relative Abundance")
#save edited phylum dataframe
write.csv(df, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv")

#adjust the column where sample is separated by "." instead of "-"
phylum.df.long = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv")
phylum.df.long = phylum.df.long %>%
  select(2:4)

#pull metadata
df.full = sample_data(ps.hell)
#merge phylum table with metadata
phylum.df.long = left_join(phylum.df.long, df.full, by ="SampleID")
phylum.df.long = select(phylum.df.long, 1:3, 9:26, 28:31)
#save csv
write.csv(phylum.df.long, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv")

#cleanup
rm(df)
rm(df.full)
rm(phylum.df.long)



#######Create OTU table 
write.csv(ps.hell.raw %>% 
            psmelt() %>%
            arrange(OTU) %>%
            select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
            spread(Sample, Abundance),
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_HellingerAll.csv")

write.csv(ps.raw %>% 
            psmelt() %>%
            arrange(OTU) %>%
            select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance) %>%
            spread(Sample, Abundance),
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_All.csv")

raw.otu.table = ps.raw %>% 
  psmelt() %>%
  arrange(OTU) %>%
  select(OTU, Domain, Phylum, Class, Order, Family, Genus, Species, Sample, Abundance)
write.csv(raw.otu.table, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_All_long.csv")

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
genus.df.long = left_join(genus.df, df.full, by ="SampleID")
genus.df.long = select(genus.df.long, 1:3, 9:26, 28:31)
#save csv
write.csv(genus.df.long, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Genus.csv")

#cleanup
rm(df.full)
rm(genus.df.long)
rm(phylum.df)




################### Start from here #####################################

ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5) 
phylum.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Samples_OTU_Phyla.csv")
genus.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Genus.csv")
raw.otu.table = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_All_long.csv")

#######################Summarize OTU Table###################################
#Goal is to compare the taxonomy assigning using either Silva or Greengenes database
raw.genus.summary = raw.otu.table %>%
  group_by(Genus) %>%
  summarise(sum(Abundance)) #number of NA Genuses 3505780 # of uncultured 992141
raw.species.summary = raw.otu.table %>%
  group_by(Species) %>%
  summarise(sum(Abundance)) # NA 8963699




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
p = p + geom_col(data=phylum.df[phylum.df$Unit=="BL" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="BL" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


############## Warner Valley
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="WV" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="WV" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p

############# Hole
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="H" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="H" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


############Unburned
#Organic
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="UB" & phylum.df$Horizon == "O", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p
#Mineral
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Unit=="UB" & phylum.df$Horizon == "M", ], aes(x=Rep, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Plot_Name))
p


########################Ordinations colored by Plot##########################
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

#########################Genus Taxa table#####################################
write.table(ps %>% tax_glom(taxrank = "Genus") %>% 
              transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt() %>% 
              select(Genus, Sample, Abundance) %>% spread(Sample, Abundance), 
            file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/ps.relative_abundance.genus.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


