library(tidyverse)
library(ggpubr)
###################
###################
# Figure 5A pyramidal SSTR3 cilia
###################
###################
SSTR3exp <-read.csv(file.choose())   

cols <- c("pos" = "chartreuse3", 
          "neg"= "grey51")

ggplot(SSTR3exp, aes(fill=SSTR3, y=Fraction, x=Layer)) + 
  geom_bar(position="stack", stat="identity",width = 0.8) +
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = cols)+
  theme_pubr() +
  theme(axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour='black')) +
  ylab("Fraction of cells")

###################
###################
#  Figure 5D interneuron SSTR3 cilia
###################
###################

IN_AVG <-read.csv(file.choose())   

cols <- c("SSTR3+" = "chartreuse", 
          "SSTR3-"= "midnightblue")

IN_AVG$Subtype <- factor(IN_AVG$Subtype, levels = c("EXC", "GAD67","ChAT","PV","SOM"))

ggplot(IN_AVG, aes(fill=SSTR3, y=Fraction, x=Subtype)) + 
  scale_y_continuous(expand = c(0,0))+
  geom_bar(position="stack", stat="identity",width = 0.8) +
  scale_fill_manual(values = cols)+
  theme_pubr() +
  theme(axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour='black')) +
  ylab("Fraction of cells")


###################
###################
# Fig 5F Percent cells with SSTR3+ cilia in culture, excitatory vs inhibitory
###################
###################
bar_plain <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean, aes(color = {{ fillcol }}),
                    fill = cols1, width = 0.75, alpha = 1), 
       scale_y_continuous(expand = c(0,0)), 
       coord_cartesian(ylim = c(0,100)), 
       theme_pubr(),
       theme(legend.position = "none",
             axis.title.x = element_blank(),
             axis.text.x = element_text(size=10, colour="black"),
             axis.text.y = element_text(size=10, colour="black"),
             axis.ticks = element_line(colour="black",size=0.5))
  )
}

CILIA <-read.csv(file.choose(), header=TRUE)  

cols1 <- c("SSTR3+" = "black", "SSTR3+GAD+" =  'black')

p1<-CILIA %>% ggplot(aes(SSTR3,Percent)) + bar_plain(fillcol = Treatment) + 
  scale_colour_manual(values = cols1) + ylab("Percent") 

p1
