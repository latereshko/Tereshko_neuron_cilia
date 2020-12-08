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
###################
###################
###################
###################
#Figure 2C,F exc/inh synapse intensity ()
normint_exc <-read.csv(file.choose(), header=TRUE)  
normint_exc48 <-read.csv(file.choose(), header=TRUE)
normint_inh <-read.csv(file.choose(), header=TRUE)

normint_exc_B <- normint_exc %>% filter(Overlay %in% c("GRB_B"))
normint_inh_B <- normint_inh %>% filter(Overlay %in% c("GRB_B"))
normint_exc48_B <- normint_exc48 %>% filter(Overlay %in% c("GRB_B"))

cols8 <- c("CTL" = "grey51","shARL13b_1"= "midnightblue")
cols9 <- c("CTL" = "grey51", "shIFT88_CEP164" = "darkgoldenrod","shARL13b_2"= "deepskyblue3")

p1<-normint_exc_B %>% ggplot(aes(x = Treatment, y = NormAvgTOT)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols8) +
  ylab("Avg. Total Intensity") + coord_cartesian(ylim = c(0,2.5), clip = "off") 

p2<-normint_exc48_B %>% ggplot(aes(x = Treatment, y = NormAvgTOT)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols9) +
  ylab("Avg. Total Intensity") + coord_cartesian(ylim = c(0,2.5), clip = "off") 

p3<-normint_inh_B %>% ggplot(aes(x = Treatment, y = NormAvgTOT)) + bees_bars(fillcol = Treatment) + scale_colour_manual(values = cols8) +
  ylab("Avg. Total Intensity") + coord_cartesian(ylim = c(0,2.5), clip = "off")

p1+p2+p3

#look at data compared to normal
ggqqplot(normint_exc_B$NormAvgTOT)
ggdensity(normint_exc_B$NormAvgTOT)
shapiro.test(normint_exc_B$NormAvgTOT)

ggqqplot(normint_exc48_B$NormAvgTOT)
ggdensity(normint_exc48_B$NormAvgTOT)
shapiro.test(normint_exc48_B$NormAvgTOT)

ggqqplot(normint_inh_B$NormAvgTOT)
ggdensity(normint_inh_B$NormAvgTOT)
shapiro.test(normint_inh_B$NormAvgTOT)

#linear model 1
lm <- normint_exc_B %>% lmer(data = ., formula = NormAvgTOT~ Treatment + (1 | Dissociation))
lm.emm <- lm %>% emmeans("trt.vs.ctrl" ~ Treatment )
lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(lm)
shapiro.test(resid(lm))
qqnorm(resid(lm))
qqline(resid(lm)) 

#linear model 2
lm <- normint_exc48_B %>% lmer(data = ., formula = log(NormAvgTOT) ~ Treatment + (1 | Dissociation))
lm.emm <- lm %>% emmeans("trt.vs.ctrl" ~ Treatment )
lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(lm)
shapiro.test(resid(lm))
qqnorm(resid(lm))
qqline(resid(lm)) 


#linear model 3
lm <- normint_inh_B %>% lmer(data = ., formula = log(NormAvgTOT) ~ Treatment + (1 | Dissociation))
lm.emm <- lm %>% emmeans("trt.vs.ctrl" ~ Treatment )
lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(lm)
shapiro.test(resid(lm))
qqnorm(resid(lm))
qqline(resid(lm)) 
