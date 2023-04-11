rm(list=ls())

library(tidyverse)
library(dplyr)
library(ggplot2)
library(effects)
library(glmmTMB)
library(emmeans)
library(lme4)
library(viridis)
library(viridisLite)
library(MetBrewer)
library(vegan)
library(MetBrewer)
library(ggpubr)
library(parameters)
library(MASS)
library(DHARMa)
library(performance)
library(stats)
library(car)
library(scales)

######################################################
### OBJECTIVE 2: SEED TRAPS WITH AND WITHOUT LURES ###
######################################################

### read and clean data
setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_2")
data <- read.csv("fakedata_obj2.csv")
head(data)

data_2 <- data %>%
  pivot_longer(cols = c("piperaceae","solanaceae","urticaceae","araceae","moraceae","clusiaceae", "myrtaceae", "melastomataceae", "unknown"),
               names_to = "plant_family",
               values_to = "count")

### total count of seeds 
data$site <- as.factor(data$site)
viridis(70, option = "inferno") # check color options 
# "#000004FF" "#280B54FF" "#65156EFF" "#9F2A63FF" "#D44842FF" "#F57D15FF" "#FAC127FF"
# "#FCFFA4FF"

figure_1_allseeds <- ggplot(data, aes(x = treatment, y = total)) +
  theme_pubclean(base_size = 15) +  
  scale_color_manual(values = c("#F57D15FF", "#9F2A63FF", "#280B54FF")) + 
  geom_jitter(aes(x = treatment, y = total, color = treatment), width = 0.1, size = 3) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures", "basal" = "Basal")) +
  stat_summary(fun.data = mean_se, color = "black", size = 1) +
  ylab ("Total number of seeds") +
  xlab ("") + 
  theme(legend.position = "none")

figure_1_allseeds
ggsave(file="figure_1_allseeds_FAKE.jpg", 
       plot= figure_1_allseeds,
       width=15,height=10,units="cm",dpi=300)

### seed community 

supp.labs <- c("Control", "Basal", "Chemical lures")
names(supp.labs) <- c("control", "basal", "lures")

figure_2_seedcommunity <-ggplot(data_2, aes(x = site, fill = plant_family, y = count)) +
  theme_pubclean(base_size = 15) +
  theme(panel.border = element_blank(), strip.background = element_rect(colour="NA", fill= "NA"))+
  geom_bar(stat = "identity", color = "NA") +
  scale_fill_viridis_d(option = "magma", direction = -1, name= "Plant family", labels = c("piperaceae" ="Piperaceae",
                                                                            "solanaceae"= "Solanaceae",
                                                                            "urticaceae"= "Urticacea",
                                                                            "araceae"= "Araceae",
                                                                            "moraceae"= "Moraceae",
                                                                            "clusiaceae"="Clusiaceae",
                                                                            "myrtaceae"= "Myrtaceae",
                                                                            "melastomataceae"= "Melastomataceae",
                                                                            "unknown"="Unknown")) + 
  scale_x_discrete(labels = c("1" = "Site 1", "2" = "Site 2", "3" = "Site 3", "4" = "Site 4")) +
  ylab ("Total number of seeds") +
  xlab (" ") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "right") + 
  facet_wrap(~treatment, labeller = labeller(treatment = supp.labs))

figure_2_seedcommunity 
ggsave(file="figure_2_seedcommunity _FAKE.jpg", 
       plot= figure_2_seedcommunity,
       width=18,height=18,units="cm",dpi=300)

head(data_2)

fruitbatsp_glmm <- glmmTMB(count ~ treatment*plant_family + (1|site) + (1|collection_time), data = data_2, family = poisson)
plot(allEffects(fruitbatsp_glmm))
tab_model(fruitbatsp_glmm)
Anova(fruitbatsp_glmm)


figure_2_allseeds <- ggplot(data_2, aes(x = plant_family, y = count, color = plant_family)) +
  theme_pubclean(base_size = 15) +  
  theme(panel.border = element_blank(), strip.background = element_rect(colour="NA", fill= "NA"))+
  scale_color_viridis_d(option = "magma", direction = -1) + 
  geom_jitter(width = 0.1, size = 3) +
  scale_x_discrete(labels = c("piperaceae" ="Piperaceae", "solanaceae"= "Solanaceae", "urticaceae"= "Urticacea","araceae"= "Araceae","moraceae"= "Moraceae", "clusiaceae"="Clusiaceae", "myrtaceae"= "Myrtaceae", "melastomataceae"= "Melastomataceae", "unknown"="Unknown")) +
  ylab ("Total number of seeds") +                                                                      
  xlab ("") +
  stat_summary(fun.data = mean_se, color = "black", size = 1) +
  facet_wrap(~treatment, labeller = labeller(treatment = supp.labs)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(legend.position = "none")

figure_2_allseeds
ggsave(file="figure_2_allseeds_FAKE.jpg", 
       plot= figure_2_allseeds,
       width=18,height=18,units="cm",dpi=300)
