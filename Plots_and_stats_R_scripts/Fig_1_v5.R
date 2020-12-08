library(ggpubr)
library(tidyverse)
library(patchwork)
library(dplyr)
library(FSA)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggbeeswarm)
library(kableExtra)
###################
#plotting functions
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
# Fig 1A Percent cells with cilia, excitatory vs inhibitory
CILIA <-read.csv(file.choose(), header=TRUE)  

cols1 <- c("Excitatory" = "black", "Inhibitory" =  'black')

p1<-CILIA %>% ggplot(aes(Treatment,Percent)) + bar_plain(fillcol = Treatment) + 
  scale_colour_manual(values = cols1) + ylab("Percent") 

p1
###################
###################
# Fig 1B Lengths of cilia excitatory vs inhibitory
###################
###################
CILIA_length <-read.csv(file.choose(), header=TRUE)  

cols2 <- c("N" = "blue", "Y" =  'red')

p2<-CILIA_length %>% ggplot(aes(x = GAD., y = Length)) + bees_bars(fillcol = GAD.)+
  coord_cartesian(ylim = c(0,12), clip = "off") + scale_colour_manual(values = cols2) + scale_x_discrete(labels=c("exc","inh"))

p2

#test
wilcox.test(Length ~ GAD., data = CILIA_length)

###################
###################
#Fig1 D+E KD efficiency shARL13b normalized total intensity
###################
###################

KD <-read.csv(file.choose(), header=TRUE)  
KD <-KD %>% filter(Channel=="ARL13b")

cols <- c("CTL" = "grey51", "shARL13b_1" = 'midnightblue',"shARL13b_2" = 'deepskyblue3')

p3<-KD %>% ggplot(aes(x = Treatment, y = NormInt)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols) + 
  coord_cartesian(ylim = c(0,1.95), clip = "off") + ylab("Normalized intensity")

p4<-KD %>% ggplot(aes(x = Treatment, y = Length)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols) + 
  coord_cartesian(ylim = c(0,10), clip = "off") + ylab("Length (um)")
p3+p4

#test
#Dunn Kruskal-Wallis multiple comparison
dunnTest(NormInt ~ Treatment,
         data=KD,
         method="bh") 

#Dunn Kruskal-Wallis multiple comparison
dunnTest(Length ~ Treatment,
         data=KD,
         method="bh") 

KDint.lm <- KD %>% lmer(data = ., formula = log(NormInt) ~ Treatment + (1 | Dissociation))

KDint.lm.emm <- KDint.lm %>% emmeans("trt.vs.ctrl" ~ Treatment )
KDint.lm.emm
KDint.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>%kable_minimal()

### LENGTHS

KDlength.lm <- KD %>% lmer(data = ., formula= log(Length) ~ Treatment + (1 | Dissociation))

KDlength.lm.emm <- KDlength.lm %>% emmeans("trt.vs.ctrl" ~ Treatment )
KDlength.lm.emm
KDlength.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>%kable_minimal()

plot(KDint.lm)
shapiro.test(resid(KDint.lm))
qqnorm(resid(KDint.lm))
qqline(resid(KDint.lm)) 
###################

###################
###################
#gross morpho
###################
###################
dends <-read.csv(file.choose(), header=TRUE)  

cols <- c("CTL" = "grey51", "shARL13b_1" =  'midnightblue')

p5<-dends %>% ggplot(aes(Treatment,Total.dendritic.length)) + bees_bars(fillcol = Treatment)  + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,3500),clip = "off")+ylab("Total dendritic length (uM)")

p6<-dends %>% ggplot(aes(Treatment,Dendritic.Nodes)) + bees_bars(fillcol = Treatment)  + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,45),clip = "off")+ylab("No. dendritic nodes")
p5+p6

#test
wilcox.test(Total.dendritic.length ~ Treatment, data = dends)

wilcox.test(Dendritic.Nodes ~ Treatment, data = dends)

