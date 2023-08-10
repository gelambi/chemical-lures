######################################################
### OBJECTIVE 2: SEED RAIN WITH AND WITHOUT LURES ###
######################################################

### script 3: SEEDS COMMUNITY AND SEED GLMMs ###

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
library(sjPlot)
library(ggvenn)

### read and clean data ###
setwd("~/Desktop/lures/objective2")  # set directory to objective 2 folder
data <- read.csv("seed_data.csv")
data$week <- as.factor(data$week)
data$site_letter <- as.factor(data$site_letter)
data$treatment <- as.factor(data$treatment)
data_cleaned <- data[, c(1, 3, 5, 6:27)]
head(data_cleaned)

### Objective 2.1: characterize the seed rain in the study area

# venn diagram 
species_list <- data_cleaned %>%
  pivot_longer(cols = -c(week, site_letter, treatment), names_to = "species", values_to = "value") %>%
  group_by(treatment) %>%
  filter(value > 0) %>%
  summarise(species_present = list(unique(species))) %>%
  pull(species_present)
names(species_list) <- c("Baseline", "Control", "Chemical lures")

colors <- c("#FBB91FFF", "#9C2964FF", "#300A5BFF")
venn_diagram <- ggvenn(
  species_list,
  fill_color = colors,
  stroke_color = "darkgrey",
  stroke_size = 0.7, set_name_size = 5
)
venn_diagram
ggsave(file="Figure8.jpg",
       plot= venn_diagram,
       width=16,height=16,units="cm",dpi=300)

### evaluate the impact of total seeds: 
totalseeds_glmm <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), data = data, family = poisson)
summary(totalseeds_glmm)
plot(allEffects(totalseeds_glmm))
tab_model(totalseeds_glmm)

### evaluate the impact of seed count*plant specie: 
# pivot longer

data_longer <- data_cleaned %>%
  pivot_longer(cols = -c(week, site_letter, treatment),
               names_to = "plant_species",
               values_to = "seed_count")
head(data_longer)
filtered_df <- data_longer %>% filter(seed_count > 0)

relative_abundance_df <- filtered_df %>%
  group_by(week) %>%
  mutate(relative_abundance = seed_count / sum(seed_count))

allseeds_relative <- ggplot(relative_abundance_df, aes(x = week, y = relative_abundance, fill = plant_species)) +
  theme_classic(base_size = 15) +  
  scale_fill_viridis_d(direction = -1, name = "Seed ID",
                       labels = c("unknown_1" = "Unknown 1",
                                  "unknown_2" = "Unknown 2",
                                  "cecropia_insignis" = parse(text = "italic('Cecropia insignis')"),
                                  "unknown_3" = "Unknown 3",
                                  "unknown_4" = "Unknown 4",
                                  "vismia_sp" = parse(text = "italic('Vismia') ~'sp.'"),
                                  "piper_multiplinervium" = parse(text = "italic('Piper multiplinervium')"),
                                  "cecropia_obtusifolia" = parse(text = "italic('Cecropia obtusifolia')"),
                                  "piper_santifelicis" = parse(text = "italic('Piper santifelicis')"),
                                  "solanum_sp" = parse(text = "italic('Solanum') ~ 'sp.'"),
                                  "anturium_sp" = parse(text = "italic('Anturium') ~'sp.'"),
                                  "ficus_columbrinae" = parse(text = "italic('Ficus columbrinae')"),
                                  "poaceae" = "Poaceae",
                                  "paspalum_conjugatum" = parse(text = "italic('Paspalum conjugatum')"),
                                  "piper_umbricola" = parse(text = "italic('Piper umbricola')"),
                                  "unknown_6" = "Unknown 6",
                                  "solanaceae" = "Solanaceae",
                                  "unknown_5" = "Unknown 5",
                                  "unknown_7" = "Unknown 7",
                                  "unknown_8" = "Unknown 8",
                                  "ficus_insipida" = parse(text = "italic('Ficus insipida')"),
                                  "unknown_9" = "Unknown 9")) +
  geom_bar(stat = "identity") +
  ylab("Relative seed abundance") +
  xlab("Collection week") +
  theme(legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1))
allseeds_relative

head(filtered_df)
supp.labs <- c("Control", "Chemical lures", "Baseline")
names(supp.labs) <- c("control", "treatment", "basal")

allseeds_total <- ggplot(filtered_df, aes(x = site_letter , y = seed_count, fill = plant_species)) +
  theme_classic(base_size = 25) +
  theme(legend.text = element_text(size = 15)) + 
  theme(panel.border = element_blank(), strip.background = element_rect(colour = "NA", fill = "NA")) +
  scale_fill_viridis_d(direction = -1, name = "Seed ID",
                       labels = c("unknown_1" = "Unknown 1",
                                  "unknown_2" = "Unknown 2",
                                  "cecropia_insignis" = parse(text = "italic('Cecropia insignis')"),
                                  "unknown_3" = "Unknown 3",
                                  "unknown_4" = "Unknown 4",
                                  "vismia_sp" = parse(text = "italic('Vismia') ~'sp.'"),
                                  "piper_multiplinervium" = parse(text = "italic('Piper multiplinervium')"),
                                  "cecropia_obtusifolia" = parse(text = "italic('Cecropia obtusifolia')"),
                                  "piper_santifelicis" = parse(text = "italic('Piper santifelicis')"),
                                  "solanum_sp" = parse(text = "italic('Solanum') ~ 'sp.'"),
                                  "anturium_sp" = parse(text = "italic('Anturium') ~'sp.'"),
                                  "ficus_columbrinae" = parse(text = "italic('Ficus columbrinae')"),
                                  "poaceae" = "Poaceae",
                                  "paspalum_conjugatum" = parse(text = "italic('Paspalum conjugatum')"),
                                  "piper_umbricola" = parse(text = "italic('Piper umbricola')"),
                                  "unknown_6" = "Unknown 6",
                                  "achenio_sp" = parse(text = "italic('Achenio') ~'sp.'"),
                                  "solanaceae" = "Solanaceae",
                                  "unknown_5" = "Unknown 5",
                                  "unknown_7" = "Unknown 7",
                                  "unknown_8" = "Unknown 8",
                                  "ficus_insipida" = parse(text = "italic('Ficus insipida')"),
                                  "unknown_9" = "Unknown 9")) +
  geom_bar(stat = "identity") +
  ylab("Total seed number") +
  xlab("") +
  theme(legend.position = "right") + 
  guides(fill = guide_legend(ncol = 1)) + 
  facet_wrap(~ treatment, labeller = labeller(treatment = supp.labs))

allseeds_total
ggsave(file="allseeds_total.jpg", 
       plot= allseeds_total,
       width=18,height=16,units="cm",dpi=300)

#######################################
### ORGANIZE DATA IN PLANT FAMILIES ###
#######################################
data_family_filtered <- data_cleaned %>%
  rowwise() %>%
  mutate(Piperaceae = piper_multiplinervium + piper_umbricola + piper_santifelicis,
         Urticaceae = cecropia_insignis + cecropia_obtusifolia,
         Hypericaceae = vismia_sp,
         Poaceae = paspalum_conjugatum + poaceae,
         Moraceae = ficus_columbrinae + ficus_insipida,
         Araceae = anturium_sp,
         Solanaceae = solanaceae + solanum_sp,
         Unknown = sum(c_across(starts_with("unknown")))) %>%
  filter(!treatment %in% "basal") %>%
  select_if(function(col) any(col != 0))
head(data_family_filtered)
data_family <- data_family_filtered[, c(1:3, 18:23)]
data_family$total <- data_family$Piperaceae + data_family$Urticaceae + data_family$Poaceae + data_family$Moraceae + data_family$Solanaceae + data_family$Unknown

### SUMMARY OF TOTAL SEEDS COLLETED IN CONTROL AND TREATMENT ###
head(data_family)

total_seed_count_per_family <- data_family %>%
  pivot_longer(cols = c(Piperaceae, Urticaceae, Poaceae, Moraceae, Solanaceae, Unknown, total),
               names_to = "Plant_family",
               values_to = "value") %>%
  group_by(treatment, Plant_family) %>%
  summarize(total_seed_count = sum(value)) %>%
  spread(treatment, total_seed_count, sep = "_")
total_seed_count_per_family <- as.data.frame(total_seed_count_per_family)
total_seed_count_per_family
write.csv(total_seed_count_per_family, file = "total_seed_count_per_family.csv")

### TOTAL NUMBER OF SEEDS ###
# GLMMs
total_glmm <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
total_glmm_nb1 <- glmmTMB(total ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
total_glmm_nb2 <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
total_glmm_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
total_glmm_nb1_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
total_glmm_nb2_zi <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
total_glmm_hurdle_poisson_zi  <- glmmTMB(total ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
total_glmm_hurdle_bn2_zi <- glmmTMB(total ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
total_glmm_null <- glmmTMB(total ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_total <- AICtab(total_glmm_null,
                    total_glmm, 
                    total_glmm_nb1,
                    total_glmm_nb2,
                    total_glmm_zi,
                    #total_glmm_nb1_zi,
                    total_glmm_nb2_zi,
                    total_glmm_hurdle_poisson_zi,
                    total_glmm_hurdle_bn2_zi,
                    base=TRUE, weights=TRUE, logLik=TRUE)
AIC_total
write.csv(AIC_total, "AIC_total.csv")
tab_model(total_glmm_zi) ### tab best model
parameters(total_glmm_zi)
r2(total_glmm_zi)

# Check for overdispesion and zero inflation 
summary(total_glmm_zi)
check_overdispersion(total_glmm_zi) # No overdispersion detected.
check_zeroinflation(total_glmm_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(total_glmm_zi))
check_model(total_glmm_zi)

# calculate effect size
glmm_emmeans <-emmeans(total_glmm_zi,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

# control 31.38, lures 72.64
72.64 - 31.38 # 41.26
(41.26*100)/31.38 # 131.485
72.64/31.38

color_palette <- viridis(20, option = "G")
selected_colors <- color_palette[c(5, 10, 15, 18)]
glmm_graph_total <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = total, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Total seed count") +
  xlab ("") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)) + 
  geom_text(x = 1.5, y = 700, label = "***", size = 5)
glmm_graph_total
ggsave(file="Figure9.jpg", 
       plot= glmm_graph_total,
       width=7,height=10,units="cm",dpi=500)

### INTERACTIVE MODEL SEED COUNT * PLANT FAMILY  ###
head(data_family)
pivoted_data_family <- data_family %>%
  pivot_longer(cols = c(Piperaceae, Urticaceae, Poaceae, Moraceae, Solanaceae, Unknown), names_to = "plant_family", values_to = "count")
pivoted_data_family

total_glmm_interactive <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), data = pivoted_data_family, family = poisson)
total_glmm_interactive_nb1 <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) +  (1|week), data = pivoted_data_family, family = nbinom1(link = "log"))
total_glmm_interactive_nb2 <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), data = pivoted_data_family, family = nbinom2(link = "log"))
total_glmm_interactive_zi <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), zi = ~site_letter, data = pivoted_data_family, family = poisson)
total_glmm_interactive_nb1_zi <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), zi = ~site_letter, data = pivoted_data_family, family = nbinom1(link = "log"))
total_glmm_interactive_nb2_zi <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), zi = ~site_letter, data = pivoted_data_family, family = nbinom2(link = "log"))
total_glmm_interactive_hurdle_poisson_zi  <- glmmTMB(count ~ treatment*plant_family + (1|site_letter) + (1|week), zi = ~site_letter, data = pivoted_data_family, family = truncated_poisson)
total_glmm_interactive_hurdle_bn2_zi <- glmmTMB(count ~ treatment*plant_family + (1|site_letter)  + (1|week), zi = ~site_letter, data = pivoted_data_family, family=truncated_nbinom2)
total_glmm_interactive_null <- glmmTMB(count ~ 1, data = pivoted_data_family, family = poisson)
# compare all models
AIC_total_interactive <- AICtab(total_glmm_interactive_null,
                                total_glmm_interactive, 
                                total_glmm_interactive_nb1,
                                total_glmm_interactive_nb2,
                                total_glmm_interactive_zi,
                                total_glmm_interactive_nb1_zi,
                                #total_glmm_interactive_nb2_zi,
                                total_glmm_interactive_hurdle_poisson_zi,
                                total_glmm_interactive_hurdle_bn2_zi,
                                base=TRUE, weights=TRUE, logLik=TRUE)
AIC_total_interactive 
write.csv(AIC_total_interactive , "AIC_total_interactive")
tab_model(total_glmm_interactive_hurdle_bn2_zi) ### tab best model
parameters(total_glmm_interactive_hurdle_bn2_zi)
r2(total_glmm_interactive_hurdle_bn2_zi)
anova_interactive_seeds <- Anova(total_glmm_interactive_hurdle_bn2_zi) ### interaction is not significant
write.csv(anova_interactive_seeds, "anova_interactive_seeds.csv")

# Check for overdispesion and zero inflation 
summary(total_glmm_interactive_hurdle_bn2_zi)
check_overdispersion(total_glmm_interactive_hurdle_bn2_zi) # No overdispersion detected.
check_zeroinflation(total_glmm_interactive_hurdle_bn2_zi) # Model is underfitting zeros (probable zero-inflation).
plot(allEffects(total_glmm_interactive_hurdle_bn2_zi))
check_model(total_glmm_interactive_hurdle_bn2_zi)

### PIPERACEAE ###
# GLMMs
piperaceae_glmm <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
piperaceae_glmm_nb1 <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
piperaceae_glmm_nb2 <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
piperaceae_glmm_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
piperaceae_glmm_nb1_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
piperaceae_glmm_nb2_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
piperaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Piperaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
piperaceae_glmm_hurdle_bn2_zi <- glmmTMB(Piperaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
piperaceae_glmm_null <- glmmTMB(Piperaceae ~ 1, data = data_family, family = poisson)
# Model comparison 
AIC_piperaceae <- AICtab(piperaceae_glmm_null,
                         piperaceae_glmm, 
                         piperaceae_glmm_nb1,
                         #piperaceae_glmm_nb2,
                         piperaceae_glmm_zi,
                         #piperaceae_glmm_nb1_zi,
                         #piperaceae_glmm_nb2_zi,
                         piperaceae_glmm_hurdle_poisson_zi,
                         piperaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_piperaceae
write.csv(AIC_piperaceae, "AIC_piperaceae.csv")
tab_model(piperaceae_glmm) ### tab best model
parameters(piperaceae_glmm)
r2(piperaceae_glmm)

# Check for overdispesion and zero inflation 
summary(piperaceae_glmm)
check_overdispersion(piperaceae_glmm) # No overdispersion detected.
check_zeroinflation(piperaceae_glmm) # Model seems ok 
plot(allEffects(piperaceae_glmm))
check_model(piperaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(piperaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)

glmm_graph_piperaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = Piperaceae, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Piperaceae") +
  xlab ("")

glmm_graph_piperaceae
ggsave(file="piperaceae.jpg", 
       plot= glmm_graph_piperaceae,
       width=15,height=15,units="cm",dpi=300)

### URTUCACEAE ###
# GLMMs
Urticaceae_glmm <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Urticaceae_glmm_nb1 <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Urticaceae_glmm_nb2 <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Urticaceae_glmm_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Urticaceae_glmm_nb1_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Urticaceae_glmm_nb2_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Urticaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Urticaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Urticaceae_glmm_hurdle_bn2_zi <- glmmTMB(Urticaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Urticaceae_glmm_null <- glmmTMB(Urticaceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Urticaceae <- AICtab(Urticaceae_glmm_null,
                         Urticaceae_glmm, 
                         Urticaceae_glmm_nb1,
                         Urticaceae_glmm_nb2,
                         #Urticaceae_glmm_zi,
                         #Urticaceae_glmm_nb1_zi,
                         #Urticaceae_glmm_nb2_zi,
                         #Urticaceae_glmm_hurdle_poisson_zi,
                         #Urticaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Urticaceae
write.csv(AIC_Urticaceae, "AIC_Urticaceae.csv")
tab_model(Urticaceae_glmm_nb2) ### tab best model
parameters(Urticaceae_glmm_nb2)
r2(Urticaceae_glmm_nb2)

# Check for overdispesion and zero inflation 
summary(Urticaceae_glmm_nb2)
check_overdispersion(Urticaceae_glmm_nb2) # No overdispersion detected.
check_zeroinflation(Urticaceae_glmm_nb2) # Model seems ok 
plot(allEffects(Urticaceae_glmm_nb2))
check_model(Urticaceae_glmm_nb2)

# calculate effect size
glmm_emmeans <-emmeans(Urticaceae_glmm_nb2,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans$treatment
glmm_emmeans$response

glmm_graph_Urticaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = response)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = Urticaceae, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = response), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Urticaceae") +
  xlab ("")

glmm_graph_Urticaceae

ggsave(file="Urticaceae.jpg", 
       plot= glmm_graph_Urticaceae,
       width=15,height=15,units="cm",dpi=300)

### POACEAE ###
# GLMMs
Poaceae_glmm <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Poaceae_glmm_nb1 <- glmmTMB(Poaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Poaceae_glmm_nb2 <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Poaceae_glmm_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Poaceae_glmm_nb1_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Poaceae_glmm_nb2_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Poaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Poaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Poaceae_glmm_hurdle_bn2_zi <- glmmTMB(Poaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Poaceae_glmm_null <- glmmTMB(Poaceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Poaceae <- AICtab(Poaceae_glmm_null,
                         Poaceae_glmm, 
                         Poaceae_glmm_nb1,
                         Poaceae_glmm_nb2,
                         Poaceae_glmm_zi,
                         #Poaceae_glmm_nb1_zi,
                         #Poaceae_glmm_nb2_zi,
                         Poaceae_glmm_hurdle_poisson_zi,
                         Poaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Poaceae
write.csv(AIC_Poaceae, "AIC_Poaceae.csv")
tab_model(Poaceae_glmm) ### tab best model
parameters(Poaceae_glmm)
r2(Poaceae_glmm)

# Check for overdispesion and zero inflation 
summary(Poaceae_glmm)
check_overdispersion(Poaceae_glmm) # No overdispersion detected.
check_zeroinflation(Poaceae_glmm) # Model seems ok 
plot(allEffects(Poaceae_glmm))
check_model(Poaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Poaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans$treatment
glmm_emmeans$rate

glmm_graph_Poaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = Poaceae, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Poaceae") +
  xlab ("")

glmm_graph_Poaceae
ggsave(file="Poaceae.jpg", 
       plot= glmm_graph_Poaceae,
       width=15,height=15,units="cm",dpi=300)

### MORACEAE ###
# GLMMs
Moraceae_glmm <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Moraceae_glmm_nb1 <- glmmTMB(Moraceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Moraceae_glmm_nb2 <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Moraceae_glmm_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Moraceae_glmm_nb1_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Moraceae_glmm_nb2_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Moraceae_glmm_hurdle_poisson_zi  <- glmmTMB(Moraceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Moraceae_glmm_hurdle_bn2_zi <- glmmTMB(Moraceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
Moraceae_glmm_null <- glmmTMB(Moraceae ~ 1, data = data_family, family = poisson)
# Model comparison
AIC_Moraceae <- AICtab(Moraceae_glmm_null,
                      Moraceae_glmm, 
                      Moraceae_glmm_nb1,
                      Moraceae_glmm_nb2,
                      #Moraceae_glmm_zi,
                      #Moraceae_glmm_nb1_zi,
                      #Moraceae_glmm_nb2_zi,
                      Moraceae_glmm_hurdle_poisson_zi,
                      #Moraceae_glmm_hurdle_bn2_zi,
                      base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Moraceae
write.csv(AIC_Moraceae, "AIC_Moraceae.csv")
tab_model(Moraceae_glmm) ### tab best model
parameters(Moraceae_glmm)
r2(Moraceae_glmm)

# Check for overdispesion and zero inflation 
summary(Moraceae_glmm)
check_overdispersion(Moraceae_glmm) # No overdispersion detected.
check_zeroinflation(Moraceae_glmm) # Model seems ok 
plot(allEffects(Moraceae_glmm))
check_model(Moraceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Moraceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

glmm_graph_Moraceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = Moraceae, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Moraceae") +
  xlab ("")

glmm_graph_Moraceae
ggsave(file="Moraceae.jpg", 
       plot= glmm_graph_Moraceae,
       width=15,height=15,units="cm",dpi=300)

### SOLANACEAE ###
# GLMMs
Solanaceae_glmm <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = poisson)
Solanaceae_glmm_nb1 <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) +  (1|week), data = data_family, family = nbinom1(link = "log"))
Solanaceae_glmm_nb2 <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), data = data_family, family = nbinom2(link = "log"))
Solanaceae_glmm_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = poisson)
Solanaceae_glmm_nb1_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom1(link = "log"))
Solanaceae_glmm_nb2_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = nbinom2(link = "log"))
Solanaceae_glmm_hurdle_poisson_zi  <- glmmTMB(Solanaceae ~ treatment + (1|site_letter) + (1|week), zi = ~site_letter, data = data_family, family = truncated_poisson)
Solanaceae_glmm_hurdle_bn2_zi <- glmmTMB(Solanaceae ~ treatment + (1|site_letter)  + (1|week), zi = ~site_letter, data = data_family, family=truncated_nbinom2)
# Model comparison 
AIC_Solanaceae <- AICtab(Solanaceae_glmm_null,
                         Solanaceae_glmm, 
                         Solanaceae_glmm_nb1,
                         Solanaceae_glmm_nb2,
                         #Solanaceae_glmm_zi,
                         #Solanaceae_glmm_nb1_zi,
                         #Solanaceae_glmm_nb2_zi,
                         #Solanaceae_glmm_hurdle_poisson_zi,
                         Solanaceae_glmm_hurdle_bn2_zi,
                         base=TRUE, weights=TRUE, logLik=TRUE)
AIC_Solanaceae
write.csv(AIC_Solanaceae, "AIC_Solanaceae.csv")
tab_model(Solanaceae_glmm) ### tab best model
parameters(Solanaceae_glmm)
r2(Solanaceae_glmm)

# Check for overdispesion and zero inflation 
summary(Solanaceae_glmm)
check_overdispersion(Solanaceae_glmm) # No overdispersion detected.
check_zeroinflation(Solanaceae_glmm) # Model seems ok 
plot(allEffects(Solanaceae_glmm))
check_model(Solanaceae_glmm)

# calculate effect size
glmm_emmeans <-emmeans(Solanaceae_glmm,~treatment, type="response")
glmm_emmeans <- as.data.frame(glmm_emmeans)
glmm_emmeans

glmm_graph_Solanaceae <- ggplot(glmm_emmeans, aes(x = treatment, y = rate)) +
  theme_classic(base_size = 15) +  
  geom_jitter(data = data_family, aes(x = treatment, y = Solanaceae, color = site_letter), width = 0.1, size = 2, alpha=0.6) +
  geom_errorbar(data = glmm_emmeans, aes(x = treatment, ymin = asymp.LCL, ymax = asymp.UCL), width = 0.05) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D")) +
  geom_point(data = glmm_emmeans, aes(y = rate), shape = 16, size = 3) + 
  scale_x_discrete(labels = c("control" = "Control", "treatment" = "Chemical lures")) + 
  ylab ("Solanaceae") +
  xlab ("")

glmm_graph_Solanaceae
ggsave(file="Solanaceae.jpg", 
       plot= glmm_graph_Solanaceae,
       width=15,height=15,units="cm",dpi=300)

family_seeds <- ggarrange(glmm_graph_Moraceae,
                          glmm_graph_Poaceae,
                          glmm_graph_piperaceae,
                          glmm_graph_Solanaceae,
                          glmm_graph_Urticaceae,
                          ncol = 2,
                          nrow = 3,
                          align = "hv",
                          common.legend = TRUE,
                          legend="bottom")
family_seeds
ggsave(file="Figure10.jpg", 
       plot= family_seeds,
       width=20,height=20,units="cm",dpi=300)

### NMDS ###
head(data_family)
mmatrix <- data_family[, c(4:9)]
mmatrix <- mmatrix + 1
matrix <- as.matrix(mmatrix) # turn data frame into matrix
nmds_results <- metaMDS(matrix, 
                        distance = "bray",       # Specify distance
                        try = 300)               # Number of iterations
nmds_results$stress
plot(nmds_results, type = "t")

#extract NMDS scores (x and y coordinates)
data.scores <- as.data.frame(scores(nmds_results$points))
data.scores$site <- data_family$site_letter
data.scores$treatment <- data_family$treatment
data.scores

# PERMANOVA
adonis_site <- adonis2(matrix ~ data.scores$site, distance = "bray", perm=999)
adonis_site
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data_family$site_letter
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)
tukey <- TukeyHSD(dispersal)
tukey

color_palette <- viridis(20, option = "G")
selected_colors <- color_palette[c(5, 10, 15, 18)]

nmdsgraph_site <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2, color = site)) +
  geom_point(size = 5) + 
  theme_classic(base_size = 25) +
  scale_color_manual(values = selected_colors, name = "Site", labels = c("1" = "A", "2" = "B", "3" = "C", "4" = "D"), guide = guide_legend(nrow = 2)) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right") + 
  annotate("text", x = -0.1, y = 0.5, size = 6, label = paste("PERMDISP2, P = 0.654\nPERMANOVA, P = 0.001 ***")) +
  stat_ellipse(level = 0.95, aes(color = site)) 
nmdsgraph_site 

# PERMANOVA
adonis_t <- adonis2(matrix ~ data.scores$treatment, distance = "bray", perm=999)
adonis_t
# BETADISP
dist_matrix <- vegdist(matrix, method = "bray")
groups <- data_family$treatment
dispersal <- betadisper(dist_matrix, groups, type = c("centroid"))
anova(dispersal)

nmdsgraph_treatment <- ggplot(data.scores, aes(x = MDS1, y = MDS2, color = treatment)) +
  geom_point(size = 5) + 
  theme_classic(base_size = 25) +
  scale_color_manual(values = selected_colors, name = "Treatment", labels = c("control" = "Control", "treatment" = "Lures")) +
  ylab ("MDS2") +
  xlab ("MDS1") + 
  theme(legend.position = "right") + 
  annotate("text", x = -0.4, y = 0.7, size = 6, label = paste("PERMDISP2, P = 0.170\nPERMANOVA, P = 0.847")) +
  stat_ellipse(level = 0.95)
nmdsgraph_treatment
nmds_seeds <- ggarrange(nmdsgraph_site,
                     nmdsgraph_treatment,
                     align = "h",
                     ncol = 2,
                     nrow = 1)
nmds_seeds 
ggsave(file="nmds_seeds.jpg", 
       plot= nmds_seeds,
       width=35,height=20,units="cm",dpi=300)

### Group all community figures together
nmds_seeds_hor <- ggarrange(nmdsgraph_site,
                            nmdsgraph_treatment,
                            ncol = 1,
                            nrow = 2)

seed_community <- ggarrange(nmds_seeds_hor,
                            allseeds_total,
                            ncol = 2,
                            nrow = 1)
seed_community
ggsave(file="Figure11.jpg", 
       plot= seed_community,
       width=55,height=30,units="cm",dpi=300)

