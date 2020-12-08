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
bees_bars_dose <- function(fillcol) {
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
bees_bars_time <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean,
                    aes(color = {{ fillcol }}), fill = NA, 
                    width = 0.75, alpha = 1), 
       geom_quasirandom(aes(colour = {{ fillcol }}),
                        shape = 16, size=0.8, width = 0.15, alpha = 1), 
       stat_summary(geom = "errorbar",
                    fun.data = mean_se, width = 0.5), 
       scale_y_continuous(expand = c(0,0)), 
       facet_grid(cols = vars(Time),as.table = FALSE, switch = NULL), 
       theme_pubr(),
       theme(legend.position = "none", 
             axis.title.x = element_blank(),
             axis.text.x = element_text(size=10, colour="black"),
             axis.text.y = element_text(size=10, colour="black"),
             axis.ticks = element_line(colour="black",size=0.5))
  )
}
###################
#Figure 6B stats for Times of merck intensity
###################
###################
SSTRDrugs <- read_csv(file.choose()) 

SSTRDrugs_TimeR <- SSTRDrugs %>% filter(!Time %in% c("24H"),
                                        Label %in% c("CTL","L_2","MK_1"),
                                        Overlay %in% c("GRB_R"))

overlaylabels <- c("Shank3", "VGLUT1")
names(overlaylabels) <- c("GRB_R", "GRB_B")

cols <- c("CTL" = "grey51", 
          "L_05"= "orange1","L_1"= "darkorange","L_2"= "orangered",
          "MK_0125"= "skyblue3","MK_05"= "deepskyblue2", "MK_1"= "deepskyblue3","MK_2"= "deepskyblue4")

SSTRDrugs_TimeR %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_time(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,3),clip = "off")+xlab("Treatment ")+ ylab("Avg. Total Intensity Shank3")

#look at data 
ggqqplot(SSTRDrugs_TimeR$NormAvgTOT)
ggdensity(SSTRDrugs_TimeR$NormAvgTOT)
shapiro.test(SSTRDrugs_TimeR$NormAvgTOT)
bartlett.test(NormAvgTOT ~ Label, data = SSTRDrugs_TimeR)

# linear model w. lmer Time Shank3
timeR.lm <- SSTRDrugs_TimeR %>% 
  lmer(data = ., formula = (log(NormAvgTOT)) ~ Label*Time + (1 | Dissociation) + (1 | SLIDE))

timeR.lm.emm <- timeR.lm %>% 
  emmeans("trt.vs.ctrl" ~ Label | factor(Time))
timeR.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(timeR.lm)
shapiro.test(resid(timeR.lm))
qqnorm(resid(timeR.lm))
qqline(resid(timeR.lm)) 


############################
############################
#Figure 6C merck density Time
############################
############################
Dense <- read_csv(file.choose()) 

Dense_time <-Dense %>% filter(!Time %in% c("24H"),
               Label %in% c("CTL","L_2","MK_1"))

Dense_time %>% ggplot(aes(Label,Density)) + bees_bars_time(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,0.5),clip = "off")+xlab("Treatment")

#look at data 
ggqqplot(Dense_time$Density)
ggdensity(Dense_time$Density)
shapiro.test(Dense_time$Density)
bartlett.test(Density ~ Label, data = Dense_time)

# linear model w. lmer time denselinear model w. lmer
dense_time.lm <- Dense_time %>% 
  lmer(data = ., formula = (log(Density)) ~ Label*Time  (1 | Dissociation))

dense_time.lm.emm <- dense_time.lm %>% 
  emmeans("trt.vs.ctrl" ~ Label | factor(Time))
dense_time.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(dense_time.lm)
shapiro.test(resid(dense_time.lm))
qqnorm(resid(dense_time.lm))
qqline(resid(dense_time.lm)) 
