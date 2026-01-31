# Exploring NMDS ordinations 

#load libraries
library(ggplot2)
library(phyloseq)
library(plyr)
library(dplyr)
library(vegan)

#load phyloseq object 
ps = readRDS("~/Desktop/LAVO22_16S_Library/30042023_MiniMac/ps.3units")
sample_data(ps)

# hellinger transformation of relative abundances
ps.hell = transform_sample_counts(ps, function(x) (x / sum(x))^0.5 )

#distribution of sample total otus
#d = df %>% 
#  group_by(Sample_ID,TotalSeqs)%>%
#  summarize()%>%
#  arrange(-TotalSeqs)
#d$Sample_ID = ordered(d$Sample_ID,levels=d$Sample_ID)
#p = ggplot(d)
#p = p + geom_bar(aes(x=Sample_ID, y=TotalSeqs), stat="identity")
#p = p + geom_hline(yintercept=3000)
#p

?ordinate

hell.pcoa.ord = ordinate(ps.hell2, method="PCoA", distance="bray")

# Plotting ordination without phyloseq
x = data.frame(hell.pcoa.ord$vectors)$Axis.1
y = data.frame(hell.pcoa.ord$vectors)$Axis.2
df.ord.hell = data.frame(sample_data(ps.hell2))
df.ord.hell$PCoA1 = x
df.ord.hell$PCoA2 = y


#get values for axes
Var.Axis.1 = hell.pcoa.ord$values$Relative_eig[1]
Var.Axis.2 = hell.pcoa.ord$values$Relative_eig[2]
Axis1.label = paste("PCoA Axis 1 (",round(Var.Axis.1*100,1)," %)")
Axis2.label = paste("PCoA Axis 2 (",round(Var.Axis.2*100,1)," %)")

df.ord.hell$median_bs = as.character(df.ord.hell$median_bs)

p = ggplot(df.ord.hell)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=median_bs,shape=Horizon))
p = p + theme_bw()
p = p + facet_grid(~Unit)
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p


p = plot_ordination(hell.pcoa.ord, color="Trt", shape="Horizon", axes=c(1,2))
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Trt,shape=Horizon))
p = p + theme_bw()
p = p + facet_grid(~Unit) 
p = p + xlab(Axis1.label)
p = p + ylab(Axis2.label)
p

p = plot_ordination(hell.pcoa.ord , color = "Trt", shape = "Horizon",  axes=c(1,2))
p = p  + geom_point(size=3) #+ scale_colour_gradient(low="red",high="blue")
p = p + guides(colour = guide_legend(), shape = guide_legend(""))
p = p + theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,strip.text.x = element_text(size=14, face="bold"),
    axis.title = element_text(size=22, face="bold"),
    axis.text = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size = 22),
    strip.background = element_rect(colour="white", fill="white"))
p = p + facet_wrap(~Land_Class)
p

#
df.ord.norm$Core_burnsev.factor = factor(df.ord.norm$Core_burnsev,order=TRUE,levels=c("1","2","3","4"))
p = ggplot(df.ord.norm)
p = p + geom_point(aes(x=PCoA1,y=PCoA2,color=Core_burnsev.factor,shape=Trt))
p = p + theme_bw()
p = p + facet_wrap(~Unit)
p

