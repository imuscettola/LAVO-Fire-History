#Objective: to better understand how the variables in my dataset are correlated
#or related to each other

#Experimental Design
#Designed this experiment of three units (aka discrete areas of the park) with
#treatments of more v. fewer fires where the soil & vegetation types of the 
#plots were selected to be equal.
#(the way I have approached the analyses, it makes sense to treat the units
#separately but still talk about them all together in the discussion. This is
#further solidified by Unit being the strongest correlate with community 
#composition) with the strata 

library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)

#upload data
#metadata = read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/11252023_LAVO22_Metadata.csv")
metadata <- read.csv("~/Box/MSWhitman/LAVO_FireHistory/Data/R/02222024_LAVO22_Metadata.csv")
#minimac
#metadata = read.csv("~/Desktop/LAVO22_GG2/11252023_LAVO22_Metadata.csv")
metadata.temp = metadata %>%
  filter(str_detect(Run_Submission_Name, "^2")) %>%
  distinct(Full_Sample_Name, .keep_all = TRUE) %>%
  filter(Blank=="N") %>%
  filter(Pos_Cntrl=="N") %>%
  select(Plot_Name, Full_Sample_Name:Core, Soil_mu:Slope, SampleSev:mtbs_sev)

veg.burn.data = read.csv("~/Desktop/LAVO22_GG2/LAVO_Soil_Data_vegburn.csv")

veg.burn.data.temp = veg.burn.data %>%
  unite(Sample_Core_Name, Plot, Quadrat, sep="_", remove = FALSE)

### From tall to wide Vegetation and Seedling Data to make easier to import into 
#metadata later

##pivot vegetation data from long to wide for future use
#veg.wide = veg.burn.data.temp %>%
#  select(Sample_Core_Name, Vegetation, VegFrac) %>%
#  drop_na(VegFrac) %>%
#  pivot_wider(names_from = Vegetation, values_from=VegFrac)
##save dataframe to LAVO_FireHistory > Data > R 
#write.csv(veg.wide, "~/Desktop/LAVO22_GG2/LAVO22_VegWide.csv", row.names=FALSE)

#veg.wide = read.csv("~/Desktop/LAVO22_GG2/LAVO22_VegWide.csv")

##pivot seedling data from long to wide for future use
#seedling.wide <- veg.burn.data.temp %>%
#  select(Sample_Core_Name, Seedlings, SeedCount) %>%
#  drop_na(SeedCount) %>%
#  pivot_wider(names_from=Seedlings, values_from = SeedCount)
##save dataframe to LAVO_FireHistory > Data > R
#write.csv(seedling.wide, "~/Desktop/LAVO22_GG2/LAVO22_SeedlingsWide.csv", row.names = FALSE)

#seedling.wide = read.csv("~/Desktop/LAVO22_GG2/LAVO22_SeedlingsWide.csv")


###add O_depthcm to metadata
O_depth = veg.burn.data.temp %>%
  select(Sample_Core_Name, O_Depthcm) %>%
  drop_na(O_Depthcm)

#check number of rows in O_depth matches the number of samples in dataset 
metadata.temp.sum <- metadata.temp %>%
  count(Horizon)
View(metadata.temp.sum) #185 mineral soil samples, yes it matches!

#add O depth to metadata
data = left_join(metadata.temp, O_depth, by="Sample_Core_Name") 
#write.csv(data, "~/Desktop/LAVO22_GG2/LAVO22_SampleMetadata.csv", row.names = FALSE)


############### Table summarizing dependent variables #####################

###count O & M horizons
data.O.M.count = data %>%
  group_by(Unit, Trt, Horizon) %>%
  summarize(n = n())
#write.csv(data.O.M.count, "~/Desktop/LAVO22_GG2/LAVO22_HorizonSummary.csv")

data.O.M.count = data %>%
  mutate(M=ifelse(Horizon=="M", 1, 0)) %>%
  mutate(O=ifelse(Horizon=="O", 1, 0)) %>%
  group_by(Unit, Trt) %>%
  summarize(O.n=sum(O), M.n=sum(M))

###count Veg types
#make dataframe of plot level info only
#data.plot <- data %>%
#  distinct(Plot_Name, .keep_all = TRUE) %>%
#  select(Plot_Name:Trt, Veg_Type:Slope)

data.plot <- data %>%
  distinct(Plot_Name, .keep_all=TRUE) %>%
  select(Plot_Name, Unit:Trt, Soil_mu, Overstory_mort:Slope, Veg_Type:mtbs_sev)

#write.csv(data.plot, "~/Desktop/LAVO22_GG2/LAVO22_PlotMetadata.csv", row.names=FALSE)

data.Veg.count <- data.plot %>%
  group_by(Unit, Trt, Veg_Type) %>%
  summarise(n = n())

#write.csv(data.Veg.count, "~/Desktop/LAVO22_GG2/LAVO22_VegTypeSummary.csv")

data.Veg.count$Veg_Type=as.factor(data.Veg.count$Veg_Type)
data.Veg.count.wide <- data.Veg.count %>%
  group_by(Unit, Trt) %>%
  pivot_wider(names_from = "Veg_Type", values_from = "n")

write.csv(data.Veg.count, "~/Desktop/plane/LAVO22_VegTypeSummary.csv")
write.csv(data.Veg.count.wide, "~/Desktop/plane/LAVO22_VegTypeSummaryWide.csv")

### summarize O depth by Unit, Trt, Veg_Type
data.O.depth <- data %>%
  group_by(Unit, Trt, Veg_Type) %>%
  summarise(O.mean =mean(O_Depthcm), O.sd=sd(O_Depthcm))

data.O.dep <- data %>%
  group_by(Unit, Trt) %>%
  summarise(O.mean=mean(O_Depthcm), O.sd=sd(O_Depthcm))
View(data.O.dep)

write.csv(data.O.dep, "~/Desktop/LAVO22_GG2/LAVO22_ODepthSummary.csv")

### summarize sample burn severity
data.sample.bs <- data %>%
  drop_na(SampleSev) %>%
  group_by(Unit, Trt) %>%
  summarise(sample.mean.bs = mean(SampleSev), sample.bs.sd=sd(SampleSev))

write.csv(data.sample.bs, "~/Desktop/LAVO22_GG2/LAVO22_SampleBSSummary.csv")

### summarize weighted mean quadrat burn severity
data.quad.bs <- data %>%
  drop_na(Wt_Avg_QuadBS) %>%
  group_by(Unit, Trt) %>%
  summarise(quad.mean.bs = mean(Wt_Avg_QuadBS), quad.bs.sd=sd(Wt_Avg_QuadBS))

write.csv(data.quad.bs, "~/Desktop/LAVO22_GG2/LAVO22_QuadBSSummary.csv")

### summarize Soil types by Unit and Treatment
data.soil <- data.plot %>%
  group_by(Unit, Trt, Soil_mu) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Soil_mu, values_from = n)

### by Unit, Trt, Veg type
data.soil.veg <- data.plot %>%
  group_by(Unit, Trt, Veg_Type, Soil_mu) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = Soil_mu, values_from=n)

#merge sample level 
data.summary <- merge(data.O.dep, data.sample.bs, by.x=c("Unit", "Trt"), by.y=c("Unit", "Trt"))
data.summary <- merge(data.summary, data.quad.bs, by.x=c("Unit", "Trt"), by.y=c("Unit", "Trt"))
data.summary <- merge(data.summary, data.O.M.count, by.x=c("Unit", "Trt"), by.y=c("Unit", "Trt"))
data.summary <- merge(data.summary, data.Veg.count.wide, by.x=c("Unit", "Trt"), by.y=c("Unit", "Trt"))
data.summary <- merge(data.summary, data.soil, by.x=c("Unit", "Trt"), by.y=c("Unit", "Trt"))
data.summary <- data.summary %>%
  mutate(Sample.n=O.n+M.n)

write.csv(data.summary, "~/Desktop/LAVO22_GG2/LAVO22_MetadataSummary.csv")


###### all summary data got saved to the Box > MSWhitman > LAVO_FireHistory >
# Data > R > TableSummary Folder 



################## burn sev correlations ##################################




#figures to show counts of each dependent variable by type
temp = metadata %>%
  select(Plot_Name, Veg_Type, Soil_mu, Unit, Trt)

meta_summary = metadata %>%
  group_by(Plot_Name) %>%
  summarise(samples_n = n()) %>%
  filter(Plot_Name!="NA")

meta_summary = left_join(meta_summary, temp, by="Plot_Name")

meta_summary$Soil_mu = as.character(meta_summary$Soil_mu)

all = ggplot(meta_summary, aes(x=Soil_mu, Veg_Type)) + geom_count() + facet_wrap(~Unit)
all

BL = ggplot(meta_summary[meta_summary$Unit=="BL", ], aes(x=Soil_mu, Veg_Type)) + geom_count() + facet_wrap(~Trt)
BL


veg_soil_summary = meta_summary %>%
  group_by(Unit) %>%
  unite(veg_soil, Veg_Type, Soil_mu, sep="_") 
#%>%
  summarise(veg_soil_n = count(veg_soil))


##summarize dependent variables












