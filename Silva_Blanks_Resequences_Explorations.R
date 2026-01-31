# load libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)

#load data; phyloseq object
ps.raw = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.raw")

#transform to Hellinger transformed relative abundance
ps.hell = transform_sample_counts(ps.raw, function(x) (x / sum(x))^0.5 )

#cleanup some storage space
rm(ps.raw)

#############Only do this step once: Create Phylum dataframe####################
write.csv(ps.hell %>%
            tax_glom(taxrank = "Phylum") %>% 
            psmelt() %>%
            select(Phylum, Sample, Abundance) %>% 
            spread(Sample, Abundance),
          file = "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv")

#read the newly created dataframe
phylum.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv", header = TRUE, sep = ",")
phylum.df = phylum.df[,-1]

df <- data.frame(t(phylum.df[-1]))
colnames(df) <- phylum.df[, 1]
df$SampleID = row.names(df)
df = pivot_longer(df, cols=1:41, names_to="Phylum", values_to = "Hellinger Relative Abundance")
#save edited phylum dataframe
write.csv(df, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv")

#adjust the column where sample is separated by "." instead of "-"
phylum.df.long = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv")
phylum.df.long = phylum.df.long %>%
  select(2:4)

#pull metadata
df.full = sample_data(ps.hell)
#merge phylum table with metadata
phylum.df.long = left_join(phylum.df.long, df.full, by ="SampleID")
phylum.df.long = phylum.df.long %>%
  select(1:3, 9:26, 28:31) 
#save csv
write.csv(phylum.df.long, "~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv")

#cleanup
rm(df)
rm(df.full)
rm(phylum.df.long)


##################Start from here############################################
#load data frames
phylum.df = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/Silva/Raw_OTU_Phylum.csv")
resequenced = phylum.df %>%
  filter(Resequence =="Y")%>%
  mutate(BSA = ifelse(endsWith(resequenced$SampleID, "2"), "correct", "concentrated"))



#####################OTU tables of resequenced samples########################
p = ggplot() + theme_bw()
p = p + geom_col(data=resequenced, aes(x=BSA, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") 
p = p + facet_wrap(vars(Full_Sample_Name))
p

##################Blanks########################################################
#blanks alone
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Blank =="Y", ], aes(x=SampleID, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p

#blanks with samples from the same DNA Extract Round
p = ggplot() + theme_bw()
p = p + geom_col(data=phylum.df[phylum.df$Extract_Rd == "1", ], aes(x=SampleID, y=Hellinger.Relative.Abundance, fill=Phylum), position = "fill") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
p

#let's visualize the blanks in an ordination to more easily judge whether they are
#different from real samples
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



