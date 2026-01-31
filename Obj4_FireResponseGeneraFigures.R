#relative abundance boxplot of fire responsive genera

#load libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)

#load ps object
ps.3 <- readRDS("~/Desktop/LAVO22_GG2/phyloseqobjects/ps.3units")

#clean up taxa names so OTU does not start with a number 
ps.3 <- clean_taxa_names(ps.3)
ps.3.norm <- transform_sample_counts(ps.3, function(x) x/sum(x))

#load vector of OTUs to extract from phyloseq object
toextract <- read.csv("~/Desktop/LAVO22_GG2/corncob/toextract.csv")


#vector of OTUs to extract
taxaofint <- toextract$OTU

#extract OTUs from phyloseq object
ps.extract <- prune_taxa(taxaofint, ps.3.norm)
#check taxa table
tax_table(ps.extract)

#convert phyloseq object to dataframe
mdf <- psmelt(ps.extract)

#subset samples for OTUs where abundance was significantly different
BL_firehist <- filter(toextract, unit=="Butte Lake" & comparison=="fire history")
H_firehist <- filter(toextract, unit=="Hole" & comparison =="fire history")
WV_firehist <- filter(toextract, unit=="Warner Valley" & comparison=="fire history")

BL_bs <- filter(toextract, unit=="Butte Lake" & comparison=="burn severity")
H_bs <- filter(toextract, unit=="Hole" & comparison =="burn severity")
WV_bs <- filter(toextract, unit=="Warner Valley" & comparison=="burn severity")

BL <- filter(toextract, unit=="Butte Lake")
vec <- BL$OTU

BL_OTUs <- mdf %>%
  filter(Unit=="BL") %>%
  subset(OTU %in% BL$OTU)

H_OTUs <- mdf %>%
  filter(Unit=="H") %>%
  subset(OTU %in% toextract[toextract$unit=="Hole", ]$OTU)

WV_OTUs <- mdf %>%
  filter(Unit=="WV") %>%
  subset(OTU %in% toextract[toextract$unit=="Warner Valley", ]$OTU)

unique(BL_OTUs$OTU)
unique(H_OTUs$OTU)
unique(WV_OTUs$OTU)

finalmdf <- rbind(BL_OTUs, H_OTUs, WV_OTUs)

mlmdf <- finalmdf %>%
  subset(OTU %in% toextract[toextract$comparison=="fire history", ]$OTU)

bsmdf <- finalmdf %>%
  subset(OTU %in% toextract[toextract$comparison=="burn severity", ]$OTU)

write.csv(mlmdf, "~/Desktop/LAVO22_GG2/corncob/responsivegenera_firehistory.csv")
write.csv(bsmdf, "~/Desktop/LAVO22_GG2/corncob/responsivegenera_bs.csv")

####### VISUALIZING ##############
#plot taxa of interest
bsmdf$SampleSev <- as.factor(bsmdf$SampleSev)
fill_labs <- c("Butte Lake", "Hole", "Warner Valley")

p = ggplot(bsmdf, aes(fill=factor(Unit, labels=c("Butte Lake", "Hole", "Warner Valley"))))
p = p + geom_boxplot(position="dodge", aes(x=SampleSev,y=log(Abundance)))
p = p + facet_grid(Genus~OTU) + theme_bw()
p

#burn severity taxa of interest position dodge, unit
Mass_bs <- ggplot(bsmdf[bsmdf$Genus=="Massilia", ], aes(fill=factor(Unit, labels=c("Butte Lake", "Hole", "Warner Valley"))))
Mass_bs <- Mass_bs + geom_boxplot(position="dodge", aes(x=SampleSev,y=log(Abundance)))
Mass_bs <- Mass_bs + facet_grid(Genus~OTU) + theme_bw()
Mass_bs <- Mass_bs + labs(x="Burn Severity", y='log(Relative Abundance)', fill="Unit") 
Mass_bs <- Mass_bs + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Mass_bs

Novi_bs <- ggplot(bsmdf[startsWith(bsmdf$Genus, "Novi"), ], aes(fill=factor(Unit, labels=c("Butte Lake", "Hole"))))
Novi_bs <- Novi_bs + geom_boxplot(position="dodge", aes(x=SampleSev,y=log(Abundance)))
Novi_bs <- Novi_bs + facet_grid(~OTU) + theme_bw()
Novi_bs <- Novi_bs + labs(x="Burn Severity", y='log(Relative Abundance)', fill="Unit") 
Novi_bs <- Novi_bs + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Novi_bs

Mode_trt <- ggplot(bsmdf[bsmdf$Genus=="Modestobacter", ], aes(fill=factor(Unit, labels=c("Butte Lake", "Hole", "Warner Valley"))))
Mode_bs <- Mode_bs + geom_boxplot(position="dodge", aes(x=SampleSev,y=log(Abundance)))
Mode_bs <- Mode_bs + facet_grid(Genus~OTU) + theme_bw()
Mode_bs <- Mode_bs + labs(x="Burn Severity", y='log(Relative Abundance)', fill="Unit") 
Mode_bs <- Mode_bs + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Mode_bs

#fire history taxa of interest position dodge, unit
Mass_trt <- ggplot(mlmdf[mlmdf$Genus=="Massilia", ], aes(fill=factor(Unit, labels=c("Butte Lake", "Hole", "Warner Valley"))))
Mass_trt <- Mass_trt + geom_boxplot(position="dodge", aes(x=Trt, y=log(Abundance)))
Mass_trt <- Mass_trt + facet_grid(Genus~OTU) + theme_bw()
Mass_trt <- Mass_trt + labs(x="Fire History", y='log(Relative Abundance)', fill="Unit") 
#Mass_trt <- Mass_trt + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Mass_trt

Novi_trt <- ggplot(mlmdf[startsWith(mlmdf$Genus, "Novi"), ], aes(fill=factor(Unit, labels=c("Butte Lake", "Warner Valley")))) 
Novi_trt <- Novi_trt + geom_boxplot(position="dodge", aes(x=Trt,y=log(Abundance)))
Novi_trt <- Novi_trt + facet_grid(~OTU) + theme_bw()
Novi_trt <- Novi_trt + labs(x="Fire History", y='log(Relative Abundance)', fill="Unit") 
# <- Novi_trt + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Novi_trt

Mode_trt <- ggplot(mlmdf[mlmdf$Genus=="Modestobacter", ], aes(fill=factor(Unit, labels=c("Butte Lake", "Hole", "Warner Valley"))))
Mode_trt <- Mode_trt + geom_boxplot(position="dodge", aes(x=Trt,y=log(Abundance)))
Mode_trt <- Mode_trt + facet_grid(Genus~OTU) + theme_bw()
Mode_trt <- Mode_trt + labs(x="Fire History", y='log(Relative Abundance)', fill="Unit") 
#Mode_trt <- Mode_trt + scale_x_discrete(labels=c('unburned', 'low', 'moderate', 'severe'))
Mode_trt


