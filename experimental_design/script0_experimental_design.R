rm(list=ls())

setwd("/Users/marianagelambi/Desktop/lures/experimental_design")
library(ggplot2)
library(viridis)
library(viridisLite)
library(ggpubr)

custom_colors <- c("treatment" = "#440154FF", "control" = "#FDE725FF")
exp <- read.csv("ex_des_lures_obj1.csv")
head(exp)

exp$date <- as.factor(exp$date)
exp$site <- as.factor(exp$site)
table(exp$treatment) # control = 28 ; lures = 27

exp_graph_obj1 <- ggplot(exp, aes(x = date, y = site, color = treatment)) +
  theme_minimal(base_size = 15) +
  geom_point(size = 4) + 
  scale_color_manual(values = custom_colors,
                     name= "Treatment",
                     labels = c("control" = "Control", "treatment" = "Chemical lures")) +
  ylab ("Site") +
  xlab ("Capture night") + 
  ggtitle("Objective 1")
  
exp_graph_obj1

exp <- read.csv("ex_des_lures_obj2.csv")
head(exp)
custom_colors <- c("Baseline" = "#20A386FF", "Treatment" = "#440154FF", "Control" = "#FDE725FF")

exp$date <- as.factor(exp$date)
exp$site <- as.factor(exp$site)

exp$treatment <- factor(exp$treatment, levels = c("Control", "Treatment", "Baseline"))

exp_graph_obj2 <- ggplot(exp, aes(x = date, y = site, color = treatment)) +
  theme_minimal(base_size = 15) +
  geom_point(size = 4) + 
  scale_color_manual(values = custom_colors,
                     name= "Treatment",
                     labels = c("Control" = "Control", "Treatment" = "Chemical lures")) +
  ylab ("Site") +
  xlab ("Collection event") + 
  ggtitle("Objective 2")

exp_graph_obj2

exp_des <- ggarrange(exp_graph_obj1,
                     exp_graph_obj2,
                     align = "hv",
                     ncol = 1,
                     nrow = 2)
exp_des

ggsave(file="FigureS3.jpg", 
       plot= exp_des,
       width=16,height=16,units="cm",dpi=300)
