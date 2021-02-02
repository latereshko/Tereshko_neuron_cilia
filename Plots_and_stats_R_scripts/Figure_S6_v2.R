library(tidyverse)
library(ggpubr)
library(ggbeeswarm)
library(patchwork)
library(dplyr)
library(FSA)
library(lme4)
library(lmerTest)
library(emmeans)
library(kableExtra)
###################
#plotting functions
###################
bees_bars_dose <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean,
                    aes(fill = {{ fillcol }}),
                    width = 0.75, alpha = 1), 
       geom_quasirandom(aes(colour = {{ fillcol }}),
                        shape = 16, size=0.8, width = 0.15, alpha = 1), 
       stat_summary(geom = "errorbar",
                    fun.data = mean_se, width = 0.5), 
       scale_y_continuous(expand = c(0,0)), 
       #facet_grid(cols = vars(Dissociation),as.table = FALSE, switch = NULL), 
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
bees_bars_time <- function(fillcol) {
  list(stat_summary(geom = "bar", fun = mean,
                    aes(fill = {{ fillcol }}),
                    width = 0.75, alpha = 1), 
       geom_quasirandom(aes(colour = {{ fillcol }}),
                        shape = 16, size=0.8, width = 0.15, alpha = 1), 
       stat_summary(geom = "errorbar",
                    fun.data = mean_se, width = 0.5), 
       scale_y_continuous(expand = c(0,0)), 
       facet_grid(cols = vars(Time),as.table = FALSE, switch = NULL), 
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
cols <- c("CTL" = "grey51", 
         "L_05"= "orange1","L_1"= "darkorange","L_2"= "orangered",
         "MK_0125"= "skyblue3","MK_05"= "deepskyblue2", "MK_1"= "deepskyblue3","MK_2"= "deepskyblue4")

dots <- c("CTL" = "grey80", 
          "L_05"= "chocolate2","L_1"= "chocolate3","L_2"= "firebrick",
          "MK_0125"= "skyblue","MK_05"= "skyblue1", "MK_1"= "deepskyblue","MK_2"= "dodgerblue")
############################
#Figure S4A merck time VGLUT1 intensity
############################
SSTRDrugs <- read_csv(file.choose()) 

SSTRDrugs_TimeB <- SSTRDrugs %>% filter(!Time %in% c("24H"),
                                        Label %in% c("CTL","L_2","MK_1"),
                                        Overlay %in% c("GRB_B"))

SSTRDrugs_TimeB %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_time(fillcol = Label) + 
  coord_cartesian(ylim = c(0,3))+xlab("Treatment")+ ylab("Avg. Total Intensity VGLUT1")

#look at data 
ggqqplot(SSTRDrugs_TimeB,"NormAvgTOT",facet.by = "Label")
ggdensity(SSTRDrugs_TimeB,"NormAvgTOT",color = "Label",palette = cols)

############################
# linear model Time VGLUT1
timeB.lm <- SSTRDrugs_TimeB %>% 
  lmer(data = ., formula = (log(NormAvgTOT)) ~ Label*Time + (1 | Dissociation))
timeB.lm.emm <- timeB.lm %>%  emmeans("trt.vs.ctrl" ~ Label | factor(Time))
timeB.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(timeB.lm)
shapiro.test(resid(timeB.lm))
qqnorm(resid(timeB.lm))
qqline(resid(timeB.lm)) 

############################
############################
#Figure S4B Merck doses 24h
############################
############################
SSTRDrugs_doseR <- SSTRDrugs %>% filter(Time %in% c("24H"),
                                        !Label %in% c("MK_2"),
                                        Overlay %in% c("GRB_R"))

SSTRDrugs_doseB <- SSTRDrugs %>% filter(Time %in% c("24H"),
                                        !Label %in% c("MK_2"),
                                        Overlay %in% c("GRB_B"))

overlaylabels <- c("Shank3", "VGLUT1")
names(overlaylabels) <- c("GRB_R", "GRB_B")

SSTRDrugs_doseR %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_dose(fillcol = Label) + 
  coord_cartesian(ylim = c(0,3))+xlab("Treatment ")+ ylab("Avg. Total Intensity Shank3")

SSTRDrugs_doseB %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_dose(fillcol = Label) + 
  coord_cartesian(ylim = c(0,3))+xlab("Treatment ")+ ylab("Avg. Total Intensity VGLUT1")


#look at data 
ggdensity(SSTRDrugs_doseR,"NormAvgTOT",color = "Label",palette = cols)
ggdensity(SSTRDrugs_doseB,"NormAvgTOT",color = "Label",palette = cols)

# linear model dose Shank3
merck_R.lm <- SSTRDrugs_doseR %>% 
  lmer(data = ., formula = log(NormAvgTOT) ~ Label + ( 1 | Dissociation))

merck_R.lm.emm <- merck_R.lm %>% 
  emmeans("trt.vs.ctrl" ~ Label)
merck_R.lm.emm
merck_R.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(merck_R.lm)
shapiro.test(resid(merck_R.lm))
qqnorm(resid(merck_R.lm))
qqline(resid(merck_R.lm)) 

############################
# linear model dose VGLUT1
merck_B.lm <- SSTRDrugs_doseB %>% 
  lmer(data = ., formula = log(NormAvgTOT) ~ Label +(1 | SLIDE) + ( 1 | Dissociation))

merck_B.lm.emm <- merck_B.lm %>% 
  emmeans("trt.vs.ctrl" ~ Label)
merck_B.lm.emm
merck_B.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(merck_B.lm)
shapiro.test(resid(merck_B.lm))
qqnorm(resid(merck_B.lm))
qqline(resid(merck_B.lm)) 

############################
############################
# Figure S4C merck density doses 24h
############################
############################
Dense <- read_csv(file.choose()) 

Dense_dose <- Dense %>% filter(Time %in% c("24H"),
                               !Label %in% c("MK_2"))

Dense_dose %>% ggplot(aes(Label,Density)) + bees_bars_time(fillcol = Label) + 
  coord_cartesian(ylim = c(0,0.5))+xlab("Treatment")


#look at data 
ggqqplot(Dense_dose,"Density",facet.by = "Label")
ggdensity(Dense_dose,"Density",color = "Label",palette = cols)

# linear model w. lmer Dose density
dense_dose.lm <- Dense_dose %>% 
  lmer(data = ., formula = log(Density) ~ Label +(1 | SLIDE) + ( 1 | Dissociation))

dense_dose.lm.emm <- dense_dose.lm %>% 
  emmeans("trt.vs.ctrl" ~ Label)
dense_dose.lm.emm
dense_dose.lm.emm$contrasts %>% 
  rbind(adjust = "dunnett") %>% 
  kbl() %>% kable_minimal()

plot(dense_dose.lm)
shapiro.test(resid(dense_dose.lm))
qqnorm(resid(dense_dose.lm))
qqline(resid(dense_dose.lm)) 


############################
############################
#Figure S4 D merck cilia lengths
############################
############################

merck_cilia <-read.csv(file.choose(), header=TRUE)  
merck_cilia <- merck_cilia %>% filter(!Dose %in% c("L_1"))
# merck_cilia <- merck_cilia %>% filter(!Time %in% c("18H"))

cols7 <- c("grey51",'chocolate1', 'purple')

merck_cilia %>% ggplot(aes(Treatment,Length)) + bees_bars_dose(fillcol = Dose) + 
  coord_cartesian(ylim = c(0,12))

dunnTest(Length ~ Treatment,
         data=merck_cilia,
         method="bh") 

############################
############################
# Figure S4E viability
############################
############################
viability_AVG <-read.csv(file.choose())   

ggplot(viability_AVG, aes(fill=Avg, y=Fraction, x=Label)) + 
  geom_bar(position="stack", stat="identity",width = 0.8) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("dead" = "chartreuse3","live"= "grey51")) +
  theme_pubr() +
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour='black')) +
  ylab("Fraction of cells")

