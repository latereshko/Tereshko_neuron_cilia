library(ggpubr)
library(tidyverse)
library(patchwork)
library(dplyr)
library(FSA)
library(ggbeeswarm)
###################
#plotting functions
###################
###################
bees_bars <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean,
                    aes(color = {{ fillcol }}), fill = NA, 
                    width = 0.75, alpha = 1), 
       geom_quasirandom(aes(colour = {{ fillcol }}),
                        shape = 16, size=0.8, width = 0.15, alpha = 1), 
       stat_summary(geom = "errorbar",
                    fun.data = mean_se, width = 0.5), 
       scale_y_continuous(expand = c(0,0)), 
       #facet_grid(cols = vars(Dissociation),as.table = FALSE, switch = NULL), 
       theme_pubr(),
       theme(legend.position = "none", 
             axis.title.x = element_blank(),
             axis.text.x = element_text(size=10, colour="black"),
             axis.text.y = element_text(size=10, colour="black"),
             axis.ticks = element_line(colour="black",size=0.5))
  )
}
###################
###################
###################
# 48H dendrites and nodes
###################
###################
dends_48 <-read.csv(file.choose(), header=TRUE)  

cols6 <- c("CTL" = "grey51", "shARL13b_2" = 'deepskyblue3', "shIFT88_CEP164" ='darkgoldenrod')

p9<-dends_48 %>% ggplot(aes(Treatment,Lengths)) + bees_bars(fillcol = Treatment)  + scale_colour_manual(values = cols6) +
  coord_cartesian(ylim = c(0,3500),clip = "off")+ylab("Total dendritic length (uM)")

p10<-dends_48 %>% ggplot(aes(Treatment,Nodes)) + bees_bars(fillcol = Treatment)  + scale_colour_manual(values = cols6) +
  coord_cartesian(ylim = c(0,60),clip = "off")+ylab("No. dendritic nodes")
p9+p10

#look at data compared to normal
ggqqplot(dends_48$Lengths)
ggdensity(dends_48$Lengths)

#Dunn Kruskal-Wallis multiple comparison
dunnTest(Nodes ~ Treatment,
         data=dends_48,
         method="bh") 
#Dunn Kruskal-Wallis multiple comparison
dunnTest(Lengths ~ Treatment,
         data=dends_48,
         method="bh") 
########################

###################
###################
# KD efficiency SHIFT88&CEP164 normalized total intensity
###################
###################
KD_IFT <-read.csv(file.choose(), header=TRUE)  
KD_IFT <-KD_IFT %>% filter(Marker=="IFT88", !Treatment %in% c("shIFT88"))

cols5 <- c("CTL" = "grey51", "shIFT88_CEP164" =  'darkgoldenrod')

p7<-KD_IFT %>% ggplot(aes(x = Treatment, y = NormInt)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols5) + 
  coord_cartesian(ylim = c(0,1.95), clip = "off") + ylab("Normalized intensity")

p8<-KD_IFT %>% ggplot(aes(x = Treatment, y = Length)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols5) + 
  coord_cartesian(ylim = c(0,10), clip = "off") + ylab("Length (um)")
p7+p8

#test
wilcox.test(NormInt ~ Treatment, data = KD_IFT)
wilcox.test(Length ~ Treatment, data = KD_IFT)
