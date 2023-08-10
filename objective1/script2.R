#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### script 2: BAT SPP. GLMMs ###

rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(effects)
library(glmmTMB)
library(sjPlot)
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
library(bbmle) 

### Load data: bats captured using mist nets 
setwd("~/Desktop/lures/objective1") # set directory to objective 1 folder
data <- read.csv("data.csv")
head(data)

### Group fruit bats by date and site
total_bats <- data %>%
  group_by(date, treatment, site) %>%
  summarise(total_bats = sum(across(matches("fruit_bats"))))
head(total_bats)
# GLMM for count data
fruitbats_glmm <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), data =  total_bats, family = nbinom2(link = "log"))
fruitbats_glmm_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = poisson)
fruitbats_glmm_nb1_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = nbinom1(link = "log"))
fruitbats_glmm_nb2_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = nbinom2(link = "log"))
fruitbats_glmm_hurdle_poisson_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = truncated_poisson)
fruitbats_glmm_hurdle_bn2_zi <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family=truncated_nbinom2)
fruitbats_glmm_zi_site <- glmmTMB(total_bats ~ treatment  + (1|site) + (1|date), zi = ~1, data =  total_bats, family = poisson)
fruitbats_glmm_null <- glmmTMB(total_bats ~ 1  + (1|site) + (1|date), data =  total_bats, family = poisson)
# Model selection
AIC_total_bats <- AICtab(fruitbats_glmm_null,
                         fruitbats_glmm,
                         fruitbats_glmm_nb1,
                         fruitbats_glmm_nb2,
                         fruitbats_glmm_zi,
                         fruitbats_glmm_nb1_zi,
                         fruitbats_glmm_nb2_zi,
                         fruitbats_glmm_hurdle_poisson_zi,
                         fruitbats_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_total_bats
write.csv(AIC_total_bats, "AIC_total_bats.csv")
tab_model(fruitbats_glmm_hurdle_bn2_zi) ### tab best model: # negative binomial 2
parameters(fruitbats_glmm_hurdle_bn2_zi)
r2(fruitbats_glmm_hurdle_bn2_zi)

glmm_emmeans <-emmeans(fruitbats_glmm_hurdle_bn2_zi, ~ treatment, type="response")
glmm_emmeans
glmm_emmeans <- as.data.frame(glmm_emmeans)

# Check for overdispesion and zero inflation 
summary(fruitbats_glmm_hurdle_bn2_zi)
check_overdispersion(fruitbats_glmm_hurdle_bn2_zi) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_hurdle_bn2_zi) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(fruitbats_glmm_hurdle_bn2_zi)) # looks like no interaction 

total_bats$prediction <- predict(fruitbats_glmm_hurdle_bn2_zi, total_bats, type="response")
fruitbatstotal <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = total_bats, aes(x = treatment, y = total_bats, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("Total count of fruit bats") +
  xlab (" ") +
  theme(legend.position = "bottom") 
fruitbatstotal

ggsave(file="Figure4.jpg", 
       plot= fruitbatstotal,
       width=7,height=10,units="cm",dpi=500)

### Pivot longer 
data_long <- data %>%
  pivot_longer(cols = -c(date, site, treatment, hours, nets),
               names_to = "bat_species", values_to = "abundance") 
data_long <- data_long  %>% filter(!bat_species %in% c("bats", "fruit_bats", "cperspicillata", "csowelli", "ccastanea")) # eliminate redundant data
head(data_long)
histogram(data_long$abundance) # lots of zeros! 

data_long <- data_long %>% filter(!bat_species %in% c("nectarivorous_bats",
                                                      "insectivorous_bats",
                                                      "desmodus_rotundus",
                                                      "artibeus_spp")) # delete bats that are not fruit bats or fruit bats that have lots of zeros
head(data_long)
histogram(data_long$abundance) 

### All fruit bats, model selection 
# Interactive GLMMs
fruitbats_glmm <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), data = data_long, family = poisson)
fruitbats_glmm_nb1 <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), data = data_long, family = nbinom1(link = "log"))
fruitbats_glmm_nb2 <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), data = data_long, family = nbinom2(link = "log"))
fruitbats_glmm_zi <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family = poisson)
fruitbats_glmm_nb1_zi <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family = nbinom1(link = "log"))
fruitbats_glmm_nb2_zi <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family = nbinom2(link = "log"))
fruitbats_glmm_hurdle_poisson_zi <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family = truncated_poisson)
fruitbats_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family=truncated_nbinom2)
fruitbats_glmm_zi_site <- glmmTMB(abundance ~ treatment*bat_species + (1|site) + (1|date), zi = ~1, data = data_long, family = poisson)
fruitbats_glmm_null <- glmmTMB(abundance ~ 1 + (1|site) + (1|date), data = data_long, family = poisson)
# Model selection
AIC_interactive <- AICtab(fruitbats_glmm_null,
       fruitbats_glmm,
       fruitbats_glmm_nb1,
       fruitbats_glmm_nb2,
       fruitbats_glmm_zi,
       fruitbats_glmm_nb1_zi,
       fruitbats_glmm_nb2_zi,
       fruitbats_glmm_hurdle_poisson_zi,
       fruitbats_glmm_hurdle_bn2_zi,
       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_interactive, "AIC_interactive.csv")

tab_model(fruitbats_glmm_nb2) ### tab best model: # negative binomial 2
anova_int <- Anova(fruitbats_glmm_nb2) # interaction is not significant 
write.csv(anova_int, "anova.csv",)
glmm_emmeans <-emmeans(fruitbats_glmm_nb2, ~ treatment*bat_species, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# Check for overdispesion and zero inflation 
summary(fruitbats_glmm_nb2)
check_overdispersion(fruitbats_glmm_nb2) # No overdispersion detected.
check_zeroinflation(fruitbats_glmm_nb2) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(fruitbats_glmm)) # looks like no interaction 

data_long$prediction <- predict(fruitbats_glmm_nb2, data_long, type="response")
glmm_graph_allbats_interaction <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_long, aes(x = treatment, y = abundance, color = site),width = 0.1, size = 2, alpha=0.6) +
  #stat_summary(fun.data = mean_ci, color = "black", size = 0.8) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab ("") +
  xlab ("") +
  facet_wrap(~bat_species, scales = "free_y")
glmm_graph_allbats_interaction

# The interaction between treatment and bat species is not statistically significant (P = 0.6056)

### Subset the dataset by bat species 

# carollia
data_long
carollia_data <- data_long %>% filter(bat_species %in% c("carollia_spp"))

# models 
carollia_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = carollia_data, family = poisson)
carollia_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = poisson)
carollia_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = nbinom1(link = "log"))
carollia_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = nbinom2(link = "log"))
carollia_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = carollia_data, family = truncated_poisson)
carollia_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = carollia_data, family=truncated_nbinom2)
carollia_glmm_null <- glmmTMB(abundance ~ 1, data = carollia_data, family = poisson)

# compare all models
check_distribution(carollia_glmm) # neg. binomial (zero-infl.) 56%

AIC_carollia <- AICtab(carollia_glmm_null,
                       carollia_glmm,
                       carollia_glmm_nb1,
                       #carollia_glmm_nb2,
                       carollia_glmm_zi,
                       #carollia_glmm_nb1_zi,
                       carollia_glmm_nb2_zi,
                       carollia_glmm_hurdle_poisson_zi,
                       carollia_glmm_hurdle_bn2_zi,
                       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_carollia, "AIC_carollia.csv")

tab_model(carollia_glmm) ### tab best model
parameters(carollia_glmm)
r2(carollia_glmm)
# Check for overdispesion and zero inflation 
summary(carollia_glmm)
check_overdispersion(carollia_glmm) # No overdispersion detected.
check_zeroinflation(carollia_glmm) # Model seems ok, ratio of observed and predicted zeros is within the tolerance range.
plot(allEffects(carollia_glmm))

# calculate effect size
glmm_emmeans <- emmeans(carollia_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# control 0.29, lures 0.60
0.60 - 0.29 # 0.31
(0.31*100)/0.29 # 106.89
0.60/0.29

# compute predictions for graph using the function predict
carollia_text <- parse(text = "italic('Carollia') ~ 'spp.'")
glmm_graph_carollia <-  ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = carollia_data, aes(x = treatment, y = abundance, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab (carollia_text) +
  xlab ("") + 
  geom_text(x = 1.5, y = 5, label = "*", size = 10)
glmm_graph_carollia
ggsave(file="carollia.jpg", 
       plot= glmm_graph_carollia,
       width=15,height=15,units="cm",dpi=300)

# Uroderma
uroderma_data <- data_long %>% filter(bat_species %in% c("uroderma_spp"))

# models 
uroderma_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = uroderma_data, family = poisson)
uroderma_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = poisson)
uroderma_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = nbinom1(link = "log"))
uroderma_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = nbinom2(link = "log"))
uroderma_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = uroderma_data, family = truncated_poisson)
uroderma_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = uroderma_data, family=truncated_nbinom2)
uroderma_glmm_null <- glmmTMB(abundance ~ 1, data = uroderma_data, family = poisson)
# compare all models
check_distribution(uroderma_glmm) # neg. binomial (zero-infl.) 56%

AIC_uroderma <- AICtab(uroderma_glmm_null,
                       uroderma_glmm,
                       uroderma_glmm_nb1,
                       uroderma_glmm_nb2,
                       uroderma_glmm_zi,
                       #uroderma_glmm_nb1_zi,
                       #uroderma_glmm_nb2_zi,
                       uroderma_glmm_hurdle_poisson_zi,
                       #uroderma_glmm_hurdle_bn2_zi,
                       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_uroderma, "AIC_uroderma.csv")

tab_model(uroderma_glmm_hurdle_poisson_zi) ### tab best model
parameters(uroderma_glmm_hurdle_poisson_zi)
r2(uroderma_glmm_hurdle_poisson_zi)
# Check for overdispesion and zero inflation 
summary(uroderma_glmm_hurdle_poisson_zi)
check_overdispersion(uroderma_glmm_hurdle_poisson_zi) # No overdispersion detected.
check_zeroinflation(uroderma_glmm_hurdle_poisson_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(uroderma_glmm))
check_model(uroderma_glmm)

# calculate effect size
glmm_emmeans <-emmeans(uroderma_glmm_hurdle_poisson_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans
# compute predictions for graph using the function predict
uroderma_text <- parse(text = "italic('Uroderma') ~ 'spp.'")
uroderma_data$prediction_uroderma <- predict(uroderma_glmm_hurdle_poisson_zi, uroderma_data, type="response")
glmm_graph_uroderma <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  theme_classic(base_size = 15) +  
  geom_jitter(data = uroderma_data, aes(x = treatment, y = abundance, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab (uroderma_text) +
  xlab ("") 
glmm_graph_uroderma
ggsave(file="uroderma.jpg", 
       plot= glmm_graph_uroderma,
       width=15,height=15,units="cm",dpi=300)

# Ectophylla alba
ectophylla.alba_data <- data_long %>% filter(bat_species %in% c("ectophylla.alba"))
# GLMMs
ectophylla.alba_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = ectophylla.alba_data, family = poisson)
ectophylla.alba_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = ectophylla.alba_data, family = nbinom1(link = "log"))
ectophylla.alba_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = ectophylla.alba_data, family = nbinom2(link = "log"))
ectophylla.alba_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = poisson)
ectophylla.alba_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = nbinom1(link = "log"))
ectophylla.alba_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = nbinom2(link = "log"))
ectophylla.alba_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = ectophylla.alba_data, family = truncated_poisson)
ectophylla.alba_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = ectophylla.alba_data, family=truncated_nbinom2)
ectophylla.alba_glmm_null <- glmmTMB(abundance ~ 1, data = ectophylla.alba_data, family = poisson)
# Model selection
AIC_ectophylla <- AICtab(ectophylla.alba_glmm_null,
                         ectophylla.alba_glmm,
                         ectophylla.alba_glmm_nb1,
                         ectophylla.alba_glmm_nb2,
                         #ectophylla.alba_glmm_zi,
                         #ectophylla.alba_glmm_nb1_zi,
                         #ectophylla.alba_glmm_nb2_zi,
                         ectophylla.alba_glmm_hurdle_poisson_zi,
                         ectophylla.alba_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_ectophylla, "AIC_ectophylla.csv")
tab_model(ectophylla.alba_glmm_hurdle_bn2_zi) ### tab best model
parameters(ectophylla.alba_glmm_hurdle_bn2_zi)
r2(ectophylla.alba_glmm_hurdle_bn2_zi)
# Check for overdispesion and zero inflation 
summary(ectophylla.alba_glmm_hurdle_bn2_zi)
check_overdispersion(ectophylla.alba_glmm_hurdle_bn2_zi) # No overdispersion detected.
check_zeroinflation(ectophylla.alba_glmm_hurdle_bn2_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(ectophylla.alba_glmm_hurdle_bn2_zi))
check_model(ectophylla.alba_glmm_hurdle_bn2_zi)
# calculate effect size
glmm_emmeans <-emmeans(ectophylla.alba_glmm_hurdle_bn2_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)

# compute predictions for graph using the function predict
ectophylla.alba_text <- expression(italic("Ectophylla alba"))
ectophylla.alba_data$prediction_ectophylla.alba <- predict(ectophylla.alba_glmm_hurdle_bn2_zi, ectophylla.alba_data, type="response")
glmm_graph_ectophylla.alba <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = ectophylla.alba_data, aes(x = treatment, y = abundance, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab (ectophylla.alba_text) +
  xlab ("") 
glmm_graph_ectophylla.alba
ggsave(file="ectophylla.alba.jpg", 
       plot= glmm_graph_ectophylla.alba,
       width=15,height=15,units="cm",dpi=300)

# sturnira
sturnira_data <- data_long %>% filter(bat_species %in% c("sturnira_spp"))

# GLMMs
sturnira_glmm <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = sturnira_data, family = poisson)
sturnira_glmm_nb1 <- glmmTMB(abundance ~ treatment + (1|site) +  (1|date), data = sturnira_data, family = nbinom1(link = "log"))
sturnira_glmm_nb2 <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = sturnira_data, family = poisson)
sturnira_glmm_nb1_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = sturnira_data, family = nbinom1(link = "log"))
sturnira_glmm_nb2_zi <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = sturnira_data, family = nbinom2(link = "log"))
sturnira_glmm_hurdle_poisson_zi  <- glmmTMB(abundance ~ treatment + (1|site) + (1|date), zi = ~site, data = sturnira_data, family = truncated_poisson)
sturnira_glmm_hurdle_bn2_zi <- glmmTMB(abundance ~ treatment + (1|site)  + (1|date), zi = ~site, data = sturnira_data, family=truncated_nbinom2)
sturnira_glmm_null <- glmmTMB(abundance ~ 1, data = sturnira_data, family = poisson)
# Model comparison 
AIC_sturnira <- AICtab(sturnira_glmm_null,
                       sturnira_glmm,
                       sturnira_glmm_nb1,
                       sturnira_glmm_nb2,
                       #sturnira_glmm_zi,
                       sturnira_glmm_nb1_zi,
                       #sturnira_glmm_nb2_zi,
                       #sturnira_glmm_hurdle_poisson_zi,
                       #sturnira_glmm_hurdle_bn2_zi,
                       base=TRUE, weights=TRUE, logLik=TRUE)
write.csv(AIC_sturnira, "AIC_sturnira.csv")
tab_model(sturnira_glmm) ### tab best model
parameters(sturnira_glmm)
r2(sturnira_glmm)
# Check for overdispesion and zero inflation 
summary(sturnira_glmm)
check_overdispersion(sturnira_glmm) # No overdispersion detected.
check_zeroinflation(sturnira_glmm) # Model seems ok 
plot(allEffects(sturnira_glmm))
check_model(sturnira_glmm)

# calculate effect size
glmm_emmeans <-emmeans(sturnira_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# compute predictions for graph using the function predict
sturnira_text <- parse(text = "italic('Sturnira') ~ 'spp.'")
sturnira_data$prediction_sturnira <- predict(sturnira_glmm, sturnira_data, type="response")
glmm_graph_sturnira <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = sturnira_data, aes(x = treatment, y = abundance, color = site), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.10) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_color_viridis(option="D", discrete=TRUE, name= "Sites", labels = c("site 1" = "A", "site 2" = "B", "site 3" = "C", "site 4" = "D", "site 5" = "E", "site 6" = "F")) +
  scale_x_discrete(labels = c("control" = "Control", "lures" = "Chemical lures")) + 
  ylab (sturnira_text) +
  xlab ("") 
glmm_graph_sturnira
ggsave(file="sturnira.jpg", 
       plot= glmm_graph_sturnira,
       width=15,height=15,units="cm",dpi=300)

sp_bats <- ggarrange(glmm_graph_carollia,
                     glmm_graph_ectophylla.alba,
                     glmm_graph_uroderma,
                     glmm_graph_sturnira,
                     ncol = 2,
                     nrow = 2,
                     align = "hv",
                     common.legend = TRUE,
                     legend="bottom")
sp_bats
ggsave(file="Figure5.jpg", 
       plot= sp_bats,
       width=16,height=16,units="cm",dpi=300)
