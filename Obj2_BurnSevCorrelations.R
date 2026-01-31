#Objective: 1) to determine whether burn severity is correlated with treatment for 
#a) samplesev, b) quadrat severity, c) canopy mortality, d) mtbs dnbr. 2) to determine 
#which burn severity metrics are correlated with one another a) samplesev, b) quadrat 
#severity, c) canopy mortality, d) mtbs dnbr, and e) % exposed mineral soil

###load libraries and data ####
library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)
library(lmerTest)
citation("ggplot2")
 metadata <- read.csv("~/Desktop/03112024_LAVO22_Metadata.csv")

###tidy dataframe ####
#remove run 2, blanks, positive controls, and duplicates
data <- metadata %>%
  filter(Pos_Cntrl=="N" & Blank=="N" & !endsWith(Submission_Sample_Name, "2")) %>%
  filter(startsWith(Run_Submission_Name, "1")) %>%
  select(Plot_Name, Full_Sample_Name:Core, Soil_mu:Slope, SampleSev:mtbs_sev) %>%
  filter(Unit!="UB") %>%
  mutate(Plot_Horizon=paste(Plot_Name, Horizon, sep="_"))

#average burn severity metrics for replicates at the plot
plothorizondata <- data %>%
  group_by(Plot_Horizon) %>%
  summarise(AvgSampleSev=mean(SampleSev), AvgQuadBS=mean(Wt_Avg_QuadBS), AvgBAREROCK=mean(BARE_ROCK)) 

temp <- data %>%
  select(Plot_Name, Unit:Horizon, Overstory_mort, Veg_Type:Plot_Horizon) %>%
  distinct(.keep_all = TRUE)

#dataframe for analyses
plothorizondata <- left_join(plothorizondata, temp, by = "Plot_Horizon") %>%
  mutate(unit=ifelse(Unit=="BL", "Butte Lake", ifelse(Unit=="H", "Hole", 'Warner Valley'))) %>%
  mutate(Treatment=ifelse(Trt=="LF", "less fire", "more fire"))
plothorizondata$Overstory_mort <- factor(plothorizondata$Overstory_mort, ordered=TRUE, levels=c(0, 1, 2, 3, 4))
plothorizondata$mtbs_sev <- factor(plothorizondata$mtbs_sev, ordered=TRUE, levels=c(1, 2, 3, 4))
plothorizondata$FireHistory <- ifelse(plothorizondata$Trt=="LF", "Less Fire", "More Fire")

rm(temp, data, metadata)

unitsummary <- plothorizondata %>%
  group_by(Unit) %>%
  summarise(meansample = mean(AvgSampleSev), sdsample=sd(AvgSampleSev), meanquad =mean(AvgQuadBS), sdquad=sd(AvgQuadBS), meanmtbs=mean(mtbs_sev), sdmtbs=sd(mtbs_sev), meanOverstory = mean(Overstory_mort), sdOverstory=sd(Overstory_mort))



### Objective 1: Is burn severity correlated with treatment? ####

#visualize Average Sample severity by treatment facet by unit
a <- ggplot(plothorizondata, aes(x=FireHistory, y=AvgSampleSev, fill=FireHistory)) 
a <- a + geom_boxplot() + facet_wrap(~unit) + scale_fill_manual(values=c('#0097d3','#ffbb22')) 
a <- a + labs(y="Mean Core Burn Severity", x="Fire History") + theme_bw()
a

#### a) samplesev ####

#check the distribution across units
ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgSampleSev, fill=Treatment)) + geom_histogram(binwidth = 1, position=position_dodge(width=1)) + facet_wrap(~unit) + labs(title="Distribution of Plot Average Sample Severity")


###Two-samples Wilcoxon test
BL.samplesev<-wilcox.test(AvgSampleSev~Treatment, data=plothorizondata[plothorizondata$Unit=="BL", ], exact=FALSE)
BL.samplesev$p.value ###BL p-value = 2.20e-09
WV.samplesev<-wilcox.test(AvgSampleSev~Treatment, data=plothorizondata[plothorizondata$Unit=="WV", ], exact=FALSE)
WV.samplesev$p.value #####WV p-value=0.857--> burn sev not significantly different
H.samplesev<-wilcox.test(AvgSampleSev ~ Treatment, data=plothorizondata[plothorizondata$Unit=="H", ], exact=FALSE)
H.samplesev$p.value #####0.636


## b) quadrat sev ####

#distributions across the units
ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgQuadBS, fill=Treatment)) + geom_histogram(binwidth = 1, position=position_dodge(width=1)) + facet_wrap(~unit) + labs(title="Distribution of Plot Average Quadrat Severity")


###Two-samples Wilcoxon test
BL.quadsev<-wilcox.test(AvgQuadBS~Treatment, data=plothorizondata[plothorizondata$Unit=="BL", ], exact=FALSE)
BL.quadsev$p.value ###BL p-value = 2.84e-09
WV.quadsev<-wilcox.test(AvgQuadBS~Treatment, data=plothorizondata[plothorizondata$Unit=="WV", ], exact=FALSE)
WV.quadsev$p.value #####WV p-value=0.740--> burn sev not significantly different
H.quadsev<-wilcox.test(AvgQuadBS ~ Treatment, data=plothorizondata[plothorizondata$Unit=="H", ], exact=FALSE)
H.quadsev$p.value #####0.728

## c) canopy mortality  ####

ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=Overstory_mort, fill=Treatment)) + geom_histogram(binwidth = 1, position=position_dodge(width=1)) + facet_wrap(~unit) + labs(title="Distribution of Canopy Mortality")


###Two-samples Wilcoxon test
BL.canopy<-wilcox.test(Overstory_mort~Treatment, data=plothorizondata[plothorizondata$Unit=="BL", ], exact=FALSE)
BL.canopy$p.value ###BL p-value = 5.763e-06 ***
WV.canopy<-wilcox.test(Overstory_mort~Treatment, data=plothorizondata[plothorizondata$Unit=="WV", ], exact=FALSE)
WV.canopy$p.value #####WV p-value=0.0357 *
H.canopy<-wilcox.test(Overstory_mort ~ Treatment, data=plothorizondata[plothorizondata$Unit=="H", ], exact=FALSE)
H.canopy$p.value #####0.875

## d) mtbs thresholded bs ####

ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=mtbs_sev, fill=Treatment)) + geom_histogram(binwidth = 1, position=position_dodge(width=1)) + facet_wrap(~unit) + labs(title="Distribution of MTBS Severity")


###Two-samples Wilcoxon test
BL.mtbs<-wilcox.test(mtbs_sev~Treatment, data=plothorizondata[plothorizondata$Unit=="BL", ], exact=FALSE)
BL.mtbs$p.value ###BL p-value = 0.000378 ***
WV.mtbs<-wilcox.test(mtbs_sev~Treatment, data=plothorizondata[plothorizondata$Unit=="WV", ], exact=FALSE)
WV.mtbs$p.value #####WV p-value=0.6713
H.mtbs<-wilcox.test(mtbs_sev ~ Treatment, data=plothorizondata[plothorizondata$Unit=="H", ], exact=FALSE)
H.mtbs$p.value #####0.0127 *


###Objective 2 correlation of burn severities
#a) samplesev, b) quadrat severity, c) canopy mortality, d) mtbs dnbr, and e) % exposed mineral soil

a <- ggplot(data=data[data$Horizon=="M", ], aes(x=SampleSev, y=Wt_Avg_QuadBS)) + geom_point() +geom_smooth(method="lm") + labs(x="Core Substrate Burn Severity", y="Average Quadrat Substrate Burn Severity")
a # p<2.2e-16; adj R2=0.936; y=0.301 + 0.89x
b <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgSampleSev, y=Overstory_mort)) + geom_point() +geom_smooth(method="lm") + labs(x="Core Substrate Burn Severity", y="Canopy Mortality")
b #p=3.72e-07; adj R2=0.361 y=-0.832+1.17x
c <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgSampleSev, y=mtbs_sev)) + geom_point() +geom_smooth(method="lm") + labs(x="Core Substrate Burn Severity", y="MTBS thresholded dNBR")
c #p=3.40e-05 adj R2=0.253; y=0.495+0.525x
d <- ggplot(data=data[data$Horizon=="M", ], aes(x=SampleSev, y=BARE_ROCK)) + geom_point() + labs(x="Core Substrate Burn Severity", y="% Exposed Mineral Soil")
d #p=0.225; adj r2=0.003; y=0.246+0.037x
e <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgQuadBS, y=Overstory_mort)) + geom_point() +geom_smooth(method="lm")+ labs(x="Average Quadrat Substrate Burn Severity", y="Canopy Mortality")
e #p=1.125e-08; r2=0.4402; y=-1.032+1.24x
f <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=AvgQuadBS, y=mtbs_sev)) + geom_point() +geom_smooth(method="lm")+ labs(x="Average Quadrat Substrate Burn Severity", y="MTBS thresholded dNBR")
f #p=1.127e-06; r2=0.341; y=0.336+0.582x
g <- ggplot(data=data[data$Horizon=="M", ], aes(x=Wt_Avg_QuadBS, y=BARE_ROCK)) + geom_point() + labs(x="Average Quadrat Substrate Burn Severity", y="% Exposed Mineral Soil")
g #p=0.358; r2=-0.00009; y=0.267+0.029x
h <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=Overstory_mort, y=mtbs_sev)) + geom_point() +geom_smooth(method="lm") + labs(x="Canopy Mortality", y="MTBS thresholded dNBR")
h #p=4.86e-07; r2=0.355; y=1.141+0.321x
i <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=Overstory_mort, y=AvgBAREROCK)) + geom_point() + labs(x="Canopy Mortality", y="% Exposed Mineral Soil")
i #p=0.415; r2=-0.0059; y=0.39-0.024x
j <- ggplot(data=plothorizondata[plothorizondata$Horizon=="M", ], aes(x=mtbs_sev, y=AvgBAREROCK)) + geom_point()  + labs(x="MTBS thresholded dNBR", y="% Exposed Mineral Soil")
j #p=0.118; r2=0.026; y=0.186+0.085x


#regressions
summary(lm(AvgQuadBS ~ AvgSampleSev, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(Overstory_mort ~ AvgSampleSev, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(mtbs_sev ~ AvgSampleSev, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(BARE_ROCK ~ SampleSev, data=data[data$Horizon=="M", ]))
summary(lm(Overstory_mort ~ AvgQuadBS, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(mtbs_sev ~ AvgQuadBS, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(BARE_ROCK ~ Wt_Avg_QuadBS, data=data[data$Horizon=="M", ]))
summary(lm(mtbs_sev ~ Overstory_mort, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(AvgBAREROCK ~ Overstory_mort, data=plothorizondata[plothorizondata$Horizon=="M", ]))
summary(lm(AvgBAREROCK ~ mtbs_sev, data=plothorizondata[plothorizondata$Horizon=="M", ]))

#Conclusions: sample and quadrat substrate severity are highly correlated. % exposed mineral soil is not correlated with any other burn severity metric.
#distribution of mtbs severity is skewed lower severity than the other burn severity metrics. 

a+b+c+d+e+f+g+h+i +j



