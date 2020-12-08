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
############################
#Figure S4A merck time VGLUT1 intensity
############################
SSTRDrugs <- read_csv(file.choose()) 

SSTRDrugs_TimeB <- SSTRDrugs %>% filter(!Time %in% c("24H"),
                                        Label %in% c("CTL","L_2","MK_1"),
                                        Overlay %in% c("GRB_B"))

SSTRDrugs_TimeB %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_time(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,3),clip = "off")+xlab("Treatment")+ ylab("Avg. Total Intensity VGLUT1")

#look at data 
ggqqplot(SSTRDrugs_TimeB$NormAvgTOT)
ggdensity(SSTRDrugs_TimeB$NormAvgTOT)
shapiro.test(SSTRDrugs_TimeB$NormAvgTOT)
bartlett.test(NormAvgTOT ~ Label, data = SSTRDrugs_TimeB)

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

cols <- c("CTL" = "grey51", 
          "L_05"= "orange1","L_1"= "darkorange","L_2"= "orangered",
          "MK_0125"= "skyblue3","MK_05"= "deepskyblue2", "MK_1"= "deepskyblue3","MK_2"= "deepskyblue4")

SSTRDrugs_doseR %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_dose(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,3),clip = "off")+xlab("Treatment ")+ ylab("Avg. Total Intensity Shank3")

SSTRDrugs_doseB %>% ggplot(aes(Label,NormAvgTOT)) + bees_bars_dose(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,3),clip = "off")+xlab("Treatment ")+ ylab("Avg. Total Intensity VGLUT1")

#look at data 
ggqqplot(SSTRDrugs_doseR$NormAvgTOT)
ggdensity(SSTRDrugs_doseR$NormAvgTOT)
shapiro.test(SSTRDrugs_doseR$NormAvgTOT)
bartlett.test(NormAvgTOT ~ Label, data = SSTRDrugs_doseR)

#look at data 
ggqqplot(SSTRDrugs_doseB$NormAvgTOT)
ggdensity(SSTRDrugs_doseB$NormAvgTOT)
shapiro.test(SSTRDrugs_doseB$NormAvgTOT)
bartlett.test(NormAvgTOT ~ Label, data = SSTRDrugs_doseB)


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


Dense_dose %>% ggplot(aes(Label,Density)) + bees_bars_time(fillcol = Label) + scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0,0.5),clip = "off")+xlab("Treatment")


#look at data compared to normal
ggqqplot(Dense_dose$Density)
ggdensity(Dense_dose$Density)
shapiro.test(Dense_dose$Density)

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
# merck_cilia <- merck_cilia %>% filter(!Dose %in% c("L_1"))
# merck_cilia <- merck_cilia %>% filter(!Time %in% c("18H"))

cols7 <- c("grey51",'chocolate1', 'purple')

merck_cilia %>% ggplot(aes(Treatment,Length)) + bees_bars(fillcol = Dose) + scale_colour_manual(values = cols7) +
  coord_cartesian(ylim = c(0,12),clip = "off")

dunnTest(Length ~ Treatment,
         data=merck_cilia,
         method="bh") 

############################
############################
# Figure S4E viability
############################
############################
viability_AVG <-read.csv(file.choose())   

cols <- c("dead" = "chartreuse3", 
          "live"= "grey51")

ggplot(viability_AVG, aes(fill=Avg, y=Fraction, x=Label)) + 
  geom_bar(position="stack", stat="identity",width = 0.8) +
  scale_fill_manual(values = cols)+
  theme_pubr() +
  theme(#legend.position = "none",
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour='black')) +
  ylab("Fraction of cells")

