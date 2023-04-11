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

#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### Bats captured using mist nets 
setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_1") # set directory to objective 1 folder
data <- read.csv("data.csv")
head(data)

### All fruit bats 
head(data)

# models 
fruitbats_glmm <- glmmTMB(fruit_bats ~ treatment + (1|site), data = data, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(fruit_bats ~ treatment + (1|site), data = data, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(fruit_bats ~ treatment + (1|site), data = data, family = nbinom2(link = "log"))
fruitbats_glmm_zi <- glmmTMB(fruit_bats ~ treatment + (1|site), zi = ~site, data = data, family = poisson)
fruitbats_glmm_nb1_zi <- glmmTMB(fruit_bats ~ treatment + (1|site), zi = ~site, data = data, family = nbinom1(link = "log"))
fruitbats_glmm_nb2_zi <- glmmTMB(fruit_bats ~ treatment + (1|site), zi = ~site, data = data, family = nbinom2(link = "log"))
fruitbats_glmm_hurdle_poisson_zi  <- glmmTMB(fruit_bats ~ treatment + (1|site), zi = ~site, data = data, family = truncated_poisson)
fruitbats_glmm_hurdle_bn2_zi <- glmmTMB(fruit_bats ~ treatment + (1|site), zi = ~site, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(fruitbats_glmm) # neg. binomial (zero-infl.) 56%

AICtab(fruitbats_glmm,
       fruitbats_glmm_nb1,
       fruitbats_glmm_nb2,
       fruitbats_glmm_zi,
       fruitbats_glmm_nb1_zi,
       fruitbats_glmm_nb2_zi,
       fruitbats_glmm_hurdle_poisson_zi,
       fruitbats_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(fruitbats_glmm_nb2) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(fruitbats_glmm_nb2 )
check_overdispersion(fruitbats_glmm_nb2) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_nb2) # Model is overfitting zeros.
plot(allEffects(fruitbats_glmm_nb2 ))
check_model(fruitbats_glmm_nb2 )

# calculate effect size
glmm_emmeans <-emmeans(fruitbats_glmm_nb2,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
data$prediction_fruitbats <- predict(fruitbats_glmm_nb2, data, type="response")
glmm_graph_allbats <- ggplot(data, aes(x = treatment, y = prediction_fruitbats)) +
  theme_pubclean(base_size = 15) + 
  geom_jitter(data = data, aes(x = treatment, y = fruit_bats, color = site), width = 0.1, size = 4, alpha = 0.7) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis_d(option="D", name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total number of fruit bat species/night/site") +
  xlab ("") + 
  theme(legend.position = "bottom")
glmm_graph_allbats
ggsave(file="fruitbats_b.jpg", 
       plot= glmm_graph_allbats,
       width=15,height=15,units="cm",dpi=300)

# carollia
data

# models 
carollia_glmm <- glmmTMB(carollia_spp ~ treatment +(1|date) + (1|site), data = data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(carollia_spp ~ treatment + (1|date) +(1|site), data = data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(carollia_spp ~ treatment + (1|date) + (1|site), data = data, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(carollia_spp ~ treatment +(1|date)+ (1|site), zi = ~site, data = data, family = poisson)
carollia_glmm_nb1_zi <- glmmTMB(carollia_spp ~ treatment + (1|date)+ (1|site), zi = ~site, data = data, family = nbinom1(link = "log"))
carollia_glmm_nb2_zi <- glmmTMB(carollia_spp ~ treatment + (1|date) +(1|site), zi = ~site, data = data, family = nbinom2(link = "log"))
carollia_glmm_hurdle_poisson_zi  <- glmmTMB(carollia_spp ~ treatment +(1|date)+ (1|site), zi = ~site, data = data, family = truncated_poisson)
carollia_glmm_hurdle_bn2_zi <- glmmTMB(carollia_spp ~ treatment +(1|date)+ (1|site), zi = ~site, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(carollia_glmm) # poisson (zero-infl.)         38%

AICtab(carollia_glmm,
       carollia_glmm_nb1,
       carollia_glmm_nb2,
       carollia_glmm_zi,
       carollia_glmm_nb1_zi,
       carollia_glmm_nb2_zi,
       carollia_glmm_hurdle_poisson_zi,
       carollia_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(carollia_glmm) ### tab best model

# Check for overdispesion and zero inflation 
summary(carollia_glmm)
check_overdispersion(carollia_glmm) # No overdispersion detected.
check_zeroinflation(carollia_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(carollia_glmm))
check_model(carollia_glmm)

# calculate effect size
glmm_emmeans <-emmeans(carollia_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
data$prediction_carollia <- predict(carollia_glmm, data, type="response")
glmm_graph_carollia <- ggplot(data, aes(x = treatment, y = prediction_carollia)) +
  theme_pubclean(base_size = 15) + 
  geom_jitter(data = data, aes(x = treatment, y = carollia_spp, color = site), width = 0.1, size = 4, alpha = 0.7) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Carollia spp.") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,7,by=2))
glmm_graph_carollia
ggsave(file="carollia_b.jpg", 
       plot= glmm_graph_carollia,
       width=15,height=15,units="cm",dpi=300)

# uroderma
data

# models 
uroderma_glmm <- glmmTMB(uroderma_spp ~ treatment + (1|date) +(1|site), data = data, family = poisson)
uroderma_glmm_nb1 <- glmmTMB(uroderma_spp ~ treatment + (1|date)+ (1|site), data = data, family = nbinom1(link = "log"))
uroderma_glmm_nb2 <- glmmTMB(uroderma_spp ~ treatment + (1|date) +(1|site), data = data, family = nbinom2(link = "log"))
uroderma_glmm_zi <- glmmTMB(uroderma_spp ~ treatment + (1|date) +(1|site), zi = ~1, data = data, family = poisson)
uroderma_glmm_nb1_zi <- glmmTMB(uroderma_spp ~ treatment + (1|date) +(1|site), zi = ~1, data = data, family = nbinom1(link = "log"))
uroderma_glmm_nb2_zi <- glmmTMB(uroderma_spp ~ treatment + (1|date) +(1|site), zi = ~1, data = data, family = nbinom2(link = "log"))
uroderma_glmm_hurdle_poisson_zi  <- glmmTMB(uroderma_spp ~ treatment +(1|date)+ (1|site), zi = ~1, data = data, family = truncated_poisson)
uroderma_glmm_hurdle_bn2_zi <- glmmTMB(uroderma_spp ~ treatment +(1|date)+ (1|site), zi = ~1, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(uroderma_glmm) # beta-binomial 44%

AICtab(uroderma_glmm,
       uroderma_glmm_nb1,
       uroderma_glmm_nb2,
       uroderma_glmm_zi,
       uroderma_glmm_nb1_zi,
       uroderma_glmm_nb2_zi,
       uroderma_glmm_hurdle_poisson_zi,
       uroderma_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(uroderma_glmm_nb2) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(uroderma_glmm_nb2)
check_overdispersion(uroderma_glmm_nb2) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(uroderma_glmm_nb2))
check_model(uroderma_glmm_nb2)

# calculate effect size
glmm_emmeans <-emmeans(uroderma_glmm_nb2,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
data$prediction_uroderma <- predict(uroderma_glmm_nb2, data, type="response")
glmm_graph_uroderma <- ggplot(data, aes(x = treatment, y = prediction_uroderma)) +
  theme_pubclean(base_size = 15) + 
  geom_jitter(data = data, aes(x = treatment, y = uroderma_spp, color = site), width=0.1, size = 4, alpha = 0.7) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Uroderma spp.") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,25,by=5))
glmm_graph_uroderma
ggsave(file="uroderma_a.jpg", 
       plot= glmm_graph_uroderma,
       width=15,height=15,units="cm",dpi=300)

# ectophylla alba
head(data)

# models 
ectophylla_glmm <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), data = data, family = poisson)
ectophylla_glmm_nb1 <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), data = data, family = nbinom1(link = "log"))
ectophylla_glmm_nb2 <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), data = data, family = nbinom2(link = "log"))
ectophylla_glmm_zi <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = poisson)
ectophylla_glmm_nb1_zi <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom1(link = "log"))
ectophylla_glmm_nb2_zi <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom2(link = "log"))
ectophylla_glmm_hurdle_poisson_zi  <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = truncated_poisson)
ectophylla_glmm_hurdle_bn2_zi <- glmmTMB(ectophylla.alba ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(ectophylla_glmm) # beta-binomial 47%

AICtab(ectophylla_glmm,
       ectophylla_glmm_nb1,
       ectophylla_glmm_nb2,
       ectophylla_glmm_zi,
       ectophylla_glmm_nb1_zi,
       ectophylla_glmm_nb2_zi,
       ectophylla_glmm_hurdle_poisson_zi,
       ectophylla_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(ectophylla_glmm_nb2) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(ectophylla_glmm_nb2)
check_overdispersion(ectophylla_glmm_nb2) # No overdispersion detected.
check_zeroinflation(ectophylla_glmm_nb2) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(ectophylla_glmm_nb2))
check_model(ectophylla_glmm_nb2)

# calculate effect size
glmm_emmeans <-emmeans(ectophylla_glmm_nb2,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
data$prediction_ectophylla <- predict(ectophylla_glmm_nb2, data, type="response")
glmm_graph_ectophylla <- ggplot(data, aes(x = treatment, y = prediction_ectophylla)) +
  theme_pubclean(base_size = 15) + 
  geom_jitter(data = data, aes(x = treatment, y = ectophylla.alba, color = site), width = 0.1, size = 4, alpha = 0.7) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Ectophylla alba") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,60,by=15))
glmm_graph_ectophylla
ggsave(file="ectophylla_b.jpg", 
       plot= glmm_graph_ectophylla,
       width=15,height=15,units="cm",dpi=300)

# artibeus
head(data) ##### cannot run model because it's just one data

# models 
artibeus_glmm <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), data = data, family = poisson)
artibeus_glmm_nb1 <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), data = data, family = nbinom1(link = "log"))
artibeus_glmm_nb2 <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), data = data, family = nbinom2(link = "log"))
artibeus_glmm_zi <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = poisson)
artibeus_glmm_nb1_zi <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom1(link = "log"))
artibeus_glmm_nb2_zi <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom2(link = "log"))
artibeus_glmm_hurdle_poisson_zi  <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = truncated_poisson)
artibeus_glmm_hurdle_bn2_zi <- glmmTMB(artibeus_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(artibeus_glmm) # neg. binomial (zero-infl.) 53%

AICtab(artibeus_glmm,
       artibeus_glmm_nb1,
       artibeus_glmm_nb2,
       artibeus_glmm_zi,
       artibeus_glmm_nb1_zi,
       artibeus_glmm_nb2_zi,
       artibeus_glmm_hurdle_poisson_zi,
       artibeus_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(artibeus_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(artibeus_glmm_zi)
check_overdispersion(artibeus_glmm) # No overdispersion detected.
check_zeroinflation(artibeus_glmm) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(artibeus_glmm))
check_model(artibeus_glmm)

# calculate effect size
glmm_emmeans <-emmeans(artibeus_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
artibeus_data$prediction_artibeus <- predict(artibeus_glmm, artibeus_data, type="response")
glmm_graph_artibeus <- ggplot(artibeus_data, aes(x = treatment, y = prediction_artibeus)) +
  theme_classic(base_size = 15) + 
  geom_line(data = artibeus_data, aes(x = treatment, y = values, group = site), position = position_dodge(0.06), alpha = 0.5, color = "light gray") +
  geom_jitter(data = artibeus_data, aes(x = treatment, y = values, color = site), position = position_dodge(0.06), size = 4) +
  stat_summary(fun.data = mean_se, color = "gray", size = 0.8) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Artibeus sp.") +
  xlab ("") +
  scale_y_continuous(breaks=seq(0,1,by=1))

glmm_graph_artibeus
ggsave(file="artibeus_a.jpg", 
       plot= glmm_graph_artibeus,
       width=15,height=15,units="cm",dpi=300)

# sturnira
head(data)

# models 
sturnira_glmm <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), data = data, family = poisson)
sturnira_glmm_nb1 <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), data = data, family = nbinom1(link = "log"))
sturnira_glmm_nb2 <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), data = data, family = nbinom2(link = "log"))
sturnira_glmm_zi <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = poisson)
sturnira_glmm_nb1_zi <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom1(link = "log"))
sturnira_glmm_nb2_zi <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = nbinom2(link = "log"))
sturnira_glmm_hurdle_poisson_zi  <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family = truncated_poisson)
sturnira_glmm_hurdle_bn2_zi <- glmmTMB(sturnira_spp ~ treatment + (1|site) + (1|date), zi = ~1, data = data, family=truncated_nbinom2)

# compare all models
check_distribution(sturnira_glmm) # neg. binomial (zero-infl59%

AICtab(sturnira_glmm,
       sturnira_glmm_nb1,
       sturnira_glmm_nb2,
       sturnira_glmm_zi,
       sturnira_glmm_nb1_zi,
       sturnira_glmm_nb2_zi,
       sturnira_glmm_hurdle_poisson_zi,
       sturnira_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)

tab_model(sturnira_glmm) ### tab best model: # negative binomial 1

# Check for overdispesion and zero inflation 
summary(sturnira_glmm)
check_overdispersion(sturnira_glmm) # No overdispersion detected.
check_zeroinflation(sturnira_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(sturnira_glmm))
check_model(sturnira_glmm)

# calculate effect size
glmm_emmeans <-emmeans(sturnira_glmm,~treatment, type="response")
glmm_emmeans

# compute predictions for graph using the function predict
data$prediction_sturnira <- predict(sturnira_glmm, data, type="response")
glmm_graph_sturnira <- ggplot(data, aes(x = treatment, y = prediction_sturnira)) +
  theme_pubclean(base_size = 15) + 
  geom_jitter(data = data, aes(x = treatment, y = sturnira_spp, color = site), width = 0.1, size = 4, alpha = 0.7) +
  stat_summary(fun.data = mean_se, color = "black", size = 0.9) +
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Sturnira spp.") +
  xlab ("") + 
  scale_y_continuous(breaks=seq(0,5,by=1))
glmm_graph_sturnira
ggsave(file="sturnira_b.jpg", 
       plot= glmm_graph_sturnira,
       width=15,height=15,units="cm",dpi=300)

sp_bats <- ggarrange(glmm_graph_carollia,
                     glmm_graph_ectophylla,
                     glmm_graph_uroderma,
                     glmm_graph_sturnira,
                     #glmm_graph_artibeus,
                     ncol = 2,
                     nrow = 2,
                     align = "hv",
                     common.legend = TRUE,
                     legend="none")
sp_bats

figure_2 <- ggarrange(glmm_graph_allbats,
                      sp_bats,
                      ncol = 2,
                      nrow = 1)

figure_2
ggsave(file="figure_2_b.jpg", 
       plot= figure_2,
       width=30,height=18,units="cm",dpi=300)
