rm(list=ls())

library(tidyverse)
library(ggplot2)
library(effects)
library(ggplot2)
library(glmmTMB)
library(viridis)
library(viridisLite)
library(MetBrewer)
library(dplyr)
library(vegan)
library(ggplot2)
library(MetBrewer)
library(glmmTMB) 
library(ggpubr)
library(parameters)
library(emmeans)

###############################
### VOLATILE EMISSION RATES ###
###############################

### with one filter paper strip
er <- read.csv("emission_rates_1filterpaper.csv") # er = emission rates
head(er)
er$mg_emited <- (er$g_emited)*1000
hep <- er %>%filter(volatile %in% c("2-heptanol"))
car <- er %>%filter(volatile %in% c("b-caryophyllene"))
er$time <- as.numeric(er$time)
glmm_hep <- glmmTMB(mg_emited ~ time + (1|replicate), data = hep)
summary(glmm_hep) # 0.130310 
glmm_car <- glmmTMB(mg_emited ~ time + (1|replicate), data = car)
summary(glmm_car) # 0.017915
hep$predictions <- predict(glmm_hep,type="response", re.form = NA, newdata = hep)
car$predictions <- predict(glmm_car,type="response", re.form = NA, newdata = car)
hep_emmeans <-emmeans(glmm_hep,~ time, type="response")
hep_emmeans <- as.data.frame(hep_emmeans)
er_graph <- ggplot(data = er, aes(x = time, y = mg_emited, color= volatile)) + 
  theme_classic(base_size = 15) + 
  geom_jitter(width = 0.1, size = 4, alpha = 0.7) +
  geom_line(data=hep, aes(x=time ,y=predictions)) +
  geom_line(data=car, aes(x=time ,y=predictions)) + 
  ylab ("Volatile emitted (mg)") +
  xlab ("Time (min)") 
er_graph
ggsave(file="er_graph.jpg", 
       plot= er_graph,
       width=15,height=10,units="cm",dpi=300)

### with three filter paper strips

er <- read.csv("emission_rates_3filterpaper.csv") # er = emission rates
head(er)
hep <- er %>%filter(volatile %in% c("2-heptanol"))
car <- er %>%filter(volatile %in% c("b-caryophyllene"))
er$time <- as.numeric(er$time)
glmm_hep <- glmmTMB(g_emited ~ time + (1|replicate), data = hep)
summary(glmm_hep)
glmm_car <- glmmTMB(g_emited ~ time + (1|replicate), data = car)
summary(glmm_car)
hep$predictions <- predict(glmm_hep,type="response", re.form = NA, newdata = hep)
car$predictions <- predict(glmm_car,type="response", re.form = NA, newdata = car)
hep_emmeans <-emmeans(glmm_hep,~ time, type="response")
hep_emmeans <- as.data.frame(hep_emmeans)
er_graph <- ggplot(data = er, aes(x = time, y = g_emited, color= volatile)) + 
  theme_classic(base_size = 15) + 
  geom_jitter(width = 0.1, size = 4, alpha = 0.7) +
  geom_line(data=hep, aes(x=time ,y=predictions)) +
  geom_line(data=car, aes(x=time ,y=predictions)) + 
  ylab ("Volatile emitted (mg)") +
  xlab ("Time (min)") 
er_graph
ggsave(file="er_graph_3fp.jpg", 
       plot= er_graph,
       width=15,height=10,units="cm",dpi=300)