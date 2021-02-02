library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(patchwork)
library(dplyr)
library(FSA)
###################
#plotting functions
###################
###################
bees_bars <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean,
                    aes(fill = {{ fillcol }}), 
                    width = 0.75, alpha = 1), 
       geom_quasirandom(aes(colour = {{ fillcol }}),
                        shape = 16, size=0.8, width = 0.15, alpha = 1), 
       stat_summary(geom = "errorbar",
                    fun.data = mean_se, width = 0.5), 
       scale_y_continuous(expand = c(0,0)), 
       scale_fill_manual(values = cols),
       scale_colour_manual(values = dots),
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
###################
#Figure 3A Absolute value amplitudes of mEPSCs
###################
###################
mEPSC_AVG<-read.csv(file.choose(), header=TRUE)

cols <- c("CTL" = "grey51","shARL13b_1"= "midnightblue","shARL13b_2" = 'deepskyblue3')
dots <- c("CTL" = "grey80", "shARL13b_1" = 'blue2',"shARL13b_2" = 'deepskyblue')

mEPSC_AVG %>% ggplot(aes(x = Treatment, y = AVG_AmpY)) + bees_bars(fillcol = Treatment) + 
  coord_cartesian(ylim = c(0,20)) + ylab("Abs. amplitude (pA)") 

#look at data 
ggqqplot(mEPSC_AVG,"AVG_AmpY",facet.by = "Treatment")
ggdensity(mEPSC_AVG,"AVG_AmpY",color = "Treatment",palette = cols)

#Dunn Kruskal-Wallis multiple comparison
dunnTest(AVG_AmpY ~ Treatment,
         data=mEPSC_AVG,
         method="bh") 


###################
###################
# Figure 3B_C mini CDFs
###################
###################
imi <-read.csv(file.choose(), header=TRUE)

first22_imi <- imi %>% group_by(Cell,Treatment,Dissociation) %>% slice(1:22)
CHECKminiperCell <-first22_imi %>% group_by(Cell,Treatment, Dissociation) %>% summarise(minicount=n())

first22_imi %>% 
  ggplot(aes(Abs_AmpY, colour = Treatment)) + 
  stat_ecdf(size=0.75) + 
  stat_ecdf(geom = "line", size = 0.1) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,1.05))+
  theme_classic() + 
  xlab("IMI (mS) ") + 
  ylab("Cumulative Probability") + 
  scale_color_manual(values=c('grey51','midnightblue','deepskyblue3'))  

first22_imi %>% 
  ggplot(aes(IMI_Manual, colour = Treatment)) + 
  stat_ecdf(size=0.75) + 
  stat_ecdf(geom = "line", size = 0.1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() + 
  xlab("IMI (mS) ") + 
  ylab("Cumulative Probability") + 
  scale_color_manual(values=c('grey51','midnightblue','deepskyblue3'))    

kruskal.test(Abs_AmpY ~ Treatment, data = first22_imi)
dunnTest(Abs_AmpY ~ Treatment,
         data=first22_imi,
         method="bonferroni") 

kruskal.test(IMI_Manual ~ Treatment, data = first22_imi)
dunnTest(IMI_Manual ~ Treatment,
         data=first22_imi,
         method="bonferroni") 

###################
###################
# Figure 1D,E Passive properties
###################
###################
AVG<-read.csv(file.choose(), header=TRUE)

AVG %>% ggplot(aes(x = Treatment, y = AVG_Vm)) + bees_bars(fillcol = Treatment)  +  
  coord_cartesian(ylim = c(0,-75)) + ylab("Vm")

AVG %>% ggplot(aes(x = Treatment, y = AVG_Rin)) + bees_bars(fillcol = Treatment)  +
  coord_cartesian(ylim = c(0,590)) + ylab("Rin") 


#look at data 
ggqqplot(AVG,"AVG_Vm",facet.by = "Treatment")
ggdensity(AVG,"AVG_Vm",color = "Treatment",palette = cols)

ggqqplot(AVG,"AVG_Rin",facet.by = "Treatment")
ggdensity(AVG,"AVG_Rin",color = "Treatment",palette = cols)

#test
#Dunn Kruskal-Wallis multiple comparison
dunnTest(AVG_Rin ~ Treatment,
         data=AVG,
         method="bh")
#Dunn Kruskal-Wallis multiple comparison
dunnTest(AVG_Vm ~ Treatment,
         data=AVG,
         method="bh") 
