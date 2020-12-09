library(ggpubr)
library(tidyverse)
library(patchwork)
library(dplyr)
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
       theme_pubr(),
       theme(legend.position = "none", 
             axis.title.x = element_blank(),
             axis.text.x = element_text(size=10, colour="black"),
             axis.text.y = element_text(size=10, colour="black"),
             axis.ticks = element_line(colour="black",size=0.5))
  )
}
###################
######################################
###################
# Figure 4A spontaneous FR
###################
###################
spont_FR <-read.csv(file.choose(), header=TRUE)

spont_FR_SEP <-spont_FR %>% separate(Cell, into = c("Cell", "Recording"), sep = "_00")

spont_AVGS <- spont_FR_SEP %>% group_by(Cell, Treatment, Dissociation) %>% summarize(AVGSPIKE=mean(Spikes))
spont_AVGS$Time <- rep(c(20), times = c(52))
spont_AVGS$RATE <- spont_AVGS$AVGSPIKE/spont_AVGS$Time

cols10 <- c("CTL" = "grey51", "shARL13b_1" =  'midnightblue')

spont_AVGS %>% ggplot(aes(x = Treatment, y = RATE)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols10) + 
  coord_cartesian(ylim = c(0,0.8), clip = "off") + ylab("Avg. Firing Rate (Hz)") 

spont_AVGS <- spont_AVGS %>% filter(!RATE == 0)

#tests
wilcox.test(RATE ~ Treatment, data = spont_AVGS)

###################
###################
# Figure 4B FI curve instantaneous firing
###################
###################

FI <-read.csv(file.choose(), header=TRUE)  

avg_FI <- FI %>%
  group_by(Treatment,STEP) %>%
  summarise(AVGperSTEP = mean(AVG_FR), se=se(AVG_FR))

ggplot(data=avg_FI, 
       aes(x=STEP, y=AVGperSTEP, ymin=(AVGperSTEP+se), ymax=(AVGperSTEP-se), 
           fill=Treatment)) + 
  geom_line(aes(colour=Treatment), size =1) + 
  geom_ribbon(alpha=0.15)+
  scale_color_manual(values=c('grey51','midnightblue','deepskyblue3')) +
  scale_fill_manual(values=c('grey51','midnightblue','deepskyblue3')) +
  scale_y_continuous(expand = c(0,0)) +
  theme_pubr() +
  labs(x="Current (pA)", y = "Frequency (Hz)") + 
  theme(legend.position = "none") 

#ANOVA repeated measures

lmeModel = lmer(AVG_FR ~ Treatment*(as.factor(STEP)) + (1 | Dissociation), data=FI)

anova(lmeModel)
lmeModel.emm <- lmeModel %>% emmeans("trt.vs.ctrl" ~ Treatment | factor(STEP) )
lmeModel.emm #emmeans
lmeModel.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>%kable_minimal()
