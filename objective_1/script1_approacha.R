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

#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### Bats captured using mist nets 
setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_1") # set directory to objective 1 folder
data <- read.csv("data_nodates.csv")
head(data)
data <- data  %>% filter(!bat_species %in% c("bats", "fruit_bats"))
### Bar graphs comparing bat communities control versus chemical lures 
supp.labs <- c("Control", "Chemical lures")
names(supp.labs) <- c("control", "lures")

figure_1 <-ggplot(data, aes(x = site, fill = bat_species, y = values)) +
  theme_pubclean(base_size = 15) +
  theme(panel.border = element_blank(), strip.background = element_rect(colour="NA", fill="NA"))+
  geom_bar(stat = "identity", color = "NA") +
  scale_fill_viridis_d(option = "inferno", name= "Bat ID", labels = c("artibeus_spp" = "Artibeus spp.",
                                                      "carollia_spp" = "Carollia spp.",
                                                      "desmodus_rotundus" = "Desmodus rotundus",
                                                      "ectophylla.alba" = "Ectophylla alba",
                                                      "insectivorous_bats" = "Insectivorous bats",
                                                      "nectarivorous_bats" = "Nectarivorous bats",
                                                      "sturnira_spp" = "Sturnira spp.",
                                                      "uroderma_spp" = "Uroderma spp.")) + 
  scale_x_discrete(labels = c("site 1" = "A",
                              "site 2" = "B",
                              "site 3" = "C",
                              "site 4" = "D",
                              "site 5" = "E",
                              "site 6" = "F")) + 
  ylab ("Total number of bats captured per site") +
  xlab (" ") + 
  theme(legend.position = "right") + 
  facet_wrap(~treatment,
             labeller = labeller(treatment = supp.labs))

figure_1
ggsave(file="figure_1.jpg", 
       plot= figure_1,
       width=16,height=14,units="cm",dpi=300)

# all fruit bats
data <- read.csv("data_nodates.csv")
head(data)
fruitbats_data <- data %>% filter(bat_species %in% c("fruit_bats"))

# models 
fruitbats_glmm <- glmmTMB(values ~ treatment + (1|site), data = fruitbats_data, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = fruitbats_data, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = fruitbats_data, family = nbinom2(link = "log"))

# poisson
check_overdispersion(fruitbats_glmm) # Overdispersion detected.
check_zeroinflation(fruitbats_glmm) # Model is underfitting zeros (probable zero-inflation).
check_distribution(fruitbats_glmm) # neg. binomial (zero-infl.) 84%
plot(allEffects(fruitbats_glmm))

# negative binomial 1
check_overdispersion(fruitbats_glmm_nb1) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_nb1) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(fruitbats_glmm_nb1))

# negative binomial 2
check_overdispersion(fruitbats_glmm_nb2) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_nb2) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(fruitbats_glmm_nb2))

# compare all models
AICtab(fruitbats_glmm,
       fruitbats_glmm_nb1,
       fruitbats_glmm_nb2,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(fruitbats_glmm_nb1) ### tab best model: # negative binomial 1

# calculate effect size
glmm_emmeans <-emmeans(fruitbats_glmm_nb1,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 
foldchange_fruitbats <- 14.4/12.6

# compute predictions for graph using the function predict
fruitbats_data$prediction_fruitbats <- predict(fruitbats_glmm_nb1, fruitbats_data, type="response")
glmm_graph_allbats <- ggplot(fruitbats_data, aes(x = treatment, y = prediction_fruitbats)) +
  theme_classic(base_size = 15) + 
  geom_line(data = fruitbats_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = fruitbats_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of fruit bat species") +
  xlab ("")
glmm_graph_allbats
ggsave(file="fruitbats_a.jpg", 
       plot= glmm_graph_allbats,
       width=15,height=15,units="cm",dpi=300)

# carollia spp.
data <- read.csv("data_nodates.csv")
data
carollia_data <- data %>% filter(bat_species %in% c("carollia_spp"))

# models 
carollia_glmm <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = carollia_data, family = nbinom2(link = "log")) # Model convergence problem
carollia_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = carollia_data, family = nbinom1(link = "log")) # Model convergence problem
carollia_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = carollia_data, family = nbinom2(link = "log")) # Model convergence problem

# poisson
check_overdispersion(carollia_glmm) # No overdispersion detected.
check_zeroinflation(carollia_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
check_distribution(carollia_glmm) # neg. binomial (zero-infl.)         56%
plot(allEffects(carollia_glmm))

# negative binomial 1
check_overdispersion(carollia_glmm_nb1) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb1) # Model is overfitting zeros.
plot(allEffects(carollia_glmm_nb1))

# negative binomial 2
check_overdispersion(carollia_glmm_nb2) # No overdispersion detected.
check_zeroinflation(carollia_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(carollia_glmm_nb2))

# compare all models
AICtab(carollia_glmm,
       carollia_glmm_nb1,
       carollia_glmm_nb2,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(carollia_glmm) ### tab best model: carollia_glmm

# calculate effect size
glmm_emmeans <-emmeans(carollia_glmm,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 

# compute predictions for graph using the function predict
carollia_data$prediction_carollia <- predict(carollia_glmm, carollia_data, type="response")
glmm_graph_carollia <- ggplot(carollia_data, aes(x = treatment, y = prediction_carollia)) +
  theme_classic(base_size = 15) + 
  geom_line(data = carollia_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = carollia_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of Carollia spp. bats") +
  xlab ("")
glmm_graph_carollia
ggsave(file="carollia_a.jpg", 
       plot= glmm_graph_carollia,
       width=15,height=15,units="cm",dpi=300)

# uroderma_spp
data <- read.csv("data_nodates.csv")
data
uroderma_data <- data %>% filter(bat_species %in% c("uroderma_spp"))

# models 
uroderma_glmm <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = poisson)
uroderma_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = uroderma_data, family = poisson)
uroderma_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = nbinom1(link = "log")) # Model convergence problem
uroderma_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = uroderma_data, family = nbinom1(link = "log")) # Model convergence problem
uroderma_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = uroderma_data, family = nbinom2(link = "log")) # Model convergence problem

# poisson
check_overdispersion(uroderma_glmm) # No overdispersion detected.
check_zeroinflation(uroderma_glmm) # Model is overfitting zeros.
check_distribution(uroderma_glmm) # neg. binomial (zero-infl.) 56%
plot(allEffects(uroderma_glmm))

# poisson, zi
check_overdispersion(uroderma_glmm_zi) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_zi) # Model is overfitting zeros.
plot(allEffects(uroderma_glmm))

# negative binomial 1
check_overdispersion(uroderma_glmm_nb1) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_nb1) # Model is overfitting zeros.
plot(allEffects(uroderma_glmm_nb1))

# negative binomial 2
check_overdispersion(uroderma_glmm_nb2) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_nb2) # Model is overfitting zeros.
plot(allEffects(uroderma_glmm_nb2))

# compare all models
AICtab(uroderma_glmm,
       uroderma_glmm_zi, 
       uroderma_glmm_nb1,
       uroderma_glmm_nb2,
       uroderma_glmm_nb1_zi,
       uroderma_glmm_nb2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(uroderma_glmm) ### tab best model: uroderma_glmm

# calculate effect size
glmm_emmeans <-emmeans(uroderma_glmm,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 

# compute predictions for graph using the function predict
uroderma_data$prediction_uroderma <- predict(uroderma_glmm, uroderma_data, type="response")
glmm_graph_uroderma <- ggplot(uroderma_data, aes(x = treatment, y = prediction_uroderma)) +
  theme_classic(base_size = 15) + 
  geom_line(data = uroderma_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = uroderma_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of Uroderma spp. bats") +
  xlab ("")
glmm_graph_uroderma
ggsave(file="uroderma_a.jpg", 
       plot= glmm_graph_uroderma,
       width=15,height=15,units="cm",dpi=300)

# ectophylla alba
data <- read.csv("data_nodates.csv")
data
ectophylla_data <- data %>% filter(bat_species %in% c("ectophylla.alba"))

# models 
ectophylla_glmm <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = poisson)
ectophylla_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = ectophylla_data, family = poisson)
ectophylla_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = nbinom1(link = "log")) 
ectophylla_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = ectophylla_data, family = nbinom2(link = "log"))
ectophylla_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = ectophylla_data, family = nbinom1(link = "log")) 
ectophylla_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = ectophylla_data, family = nbinom2(link = "log")) # Model convergence problem

# poisson
check_overdispersion(ectophylla_glmm) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm) # Model is overfitting zeros.
check_distribution(ectophylla_glmm) # neg. binomial (zero-infl.) 53%
plot(allEffects(ectophylla_glmm))

# poisson, zi
check_overdispersion(ectophylla_glmm_zi) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm_zi) # Model is overfitting zeros.
plot(allEffects(ectophylla_glmm))

# negative binomial 1
check_overdispersion(ectophylla_glmm_nb1) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm_nb1) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(ectophylla_glmm_nb1))

# negative binomial 2
check_overdispersion(ectophylla_glmm_nb2) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(ectophylla_glmm_nb2))

# compare all models
AICtab(ectophylla_glmm,
       ectophylla_glmm_zi, 
       ectophylla_glmm_nb1,
       ectophylla_glmm_nb2,
       ectophylla_glmm_nb1_zi,
       ectophylla_glmm_nb2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(ectophylla_glmm_nb1) ### tab best model: ectophylla_glmm_nb1

# calculate effect size
glmm_emmeans <-emmeans(ectophylla_glmm_nb1,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 

# compute predictions for graph using the function predict
ectophylla_data$prediction_ectophylla <- predict(ectophylla_glmm_nb1, ectophylla_data, type="response")
glmm_graph_ectophylla <- ggplot(ectophylla_data, aes(x = treatment, y = prediction_ectophylla)) +
  theme_classic(base_size = 15) + 
  geom_line(data = ectophylla_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = ectophylla_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of Ectophylla alba bats") +
  xlab ("")
glmm_graph_ectophylla
ggsave(file="ectophylla_a.jpg", 
       plot= glmm_graph_ectophylla,
       width=15,height=15,units="cm",dpi=300)

# sturnira
data <- read.csv("data_nodates.csv")
data
sturnira_data <- data %>% filter(bat_species %in% c("sturnira_spp"))

# models 
sturnira_glmm <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = poisson)
sturnira_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = poisson)
sturnira_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom1(link = "log")) 
sturnira_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = nbinom1(link = "log")) 
sturnira_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = nbinom2(link = "log")) # Model convergence problem

# poisson
check_overdispersion(sturnira_glmm) # No overdispersion detected.
check_zeroinflation(sturnira_glmm) # Model is overfitting zeros.
check_distribution(sturnira_glmm) # neg. binomial (zero-infl.) 53%
plot(allEffects(sturnira_glmm))

# poisson, zi
check_overdispersion(sturnira_glmm_zi) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_zi) # Model is overfitting zeros.
plot(allEffects(sturnira_glmm))

# negative binomial 1
check_overdispersion(sturnira_glmm_nb1) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_nb1) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm_nb1))

# negative binomial 2
check_overdispersion(sturnira_glmm_nb2) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm_nb2))

# compare all models
AICtab(sturnira_glmm,
       sturnira_glmm_zi, 
       sturnira_glmm_nb1,
       sturnira_glmm_nb2,
       sturnira_glmm_nb1_zi,
       sturnira_glmm_nb2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(sturnira_glmm) ### tab best model: sturnira_glmm_nb1

# calculate effect size
glmm_emmeans <-emmeans(sturnira_glmm,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 

# compute predictions for graph using the function predict
sturnira_data$prediction_sturnira <- predict(sturnira_glmm, sturnira_data, type="response")
glmm_graph_sturnira <- ggplot(sturnira_data, aes(x = treatment, y = prediction_sturnira)) +
  theme_classic(base_size = 15) + 
  geom_line(data = sturnira_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = sturnira_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of Sturnira spp. bats") +
  xlab ("")
glmm_graph_sturnira
ggsave(file="sturnira_a.jpg", 
       plot= glmm_graph_sturnira,
       width=15,height=15,units="cm",dpi=300)

# sturnira
data <- read.csv("data_nodates.csv")
data
sturnira_data <- data %>% filter(bat_species %in% c("sturnira_spp"))

# models 
sturnira_glmm <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = poisson)
sturnira_glmm_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = poisson)
sturnira_glmm_nb1 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom1(link = "log")) 
sturnira_glmm_nb2 <- glmmTMB(values ~ treatment + (1|site), data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_nb1_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = nbinom1(link = "log")) 
sturnira_glmm_nb2_zi <- glmmTMB(values ~ treatment + (1|site), zi = ~treatment, data = sturnira_data, family = nbinom2(link = "log")) # Model convergence problem

# poisson
check_overdispersion(sturnira_glmm) # No overdispersion detected.
check_zeroinflation(sturnira_glmm) # Model is overfitting zeros.
check_distribution(sturnira_glmm) # neg. binomial (zero-infl.) 53%
plot(allEffects(sturnira_glmm))

# poisson, zi
check_overdispersion(sturnira_glmm_zi) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_zi) # Model is overfitting zeros.
plot(allEffects(sturnira_glmm))

# negative binomial 1
check_overdispersion(sturnira_glmm_nb1) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_nb1) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm_nb1))

# negative binomial 2
check_overdispersion(sturnira_glmm_nb2) # No overdispersion detected.
check_zeroinflation(sturnira_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm_nb2))

# compare all models
AICtab(sturnira_glmm,
       sturnira_glmm_zi, 
       sturnira_glmm_nb1,
       sturnira_glmm_nb2,
       sturnira_glmm_nb1_zi,
       sturnira_glmm_nb2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(sturnira_glmm) ### tab best model: sturnira_glmm_nb1

# calculate effect size
glmm_emmeans <-emmeans(sturnira_glmm,~treatment, type="response")
glmm_emmeans # fold change with emmeans response value 

# compute predictions for graph using the function predict
sturnira_data$prediction_sturnira <- predict(sturnira_glmm, sturnira_data, type="response")
glmm_graph_sturnira <- ggplot(sturnira_data, aes(x = treatment, y = prediction_sturnira)) +
  theme_classic(base_size = 15) + 
  geom_line(data = sturnira_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = sturnira_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of Sturnira spp. bats") +
  xlab ("")
glmm_graph_sturnira
ggsave(file="sturnira_a.jpg", 
       plot= glmm_graph_sturnira,
       width=15,height=15,units="cm",dpi=300)


