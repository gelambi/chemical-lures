rm(list=ls())

library(tidyverse)
library(dplyr)

#############################################################
### OBJECTIVE 1: FRUIT BAT CAPTURE WITH AND WITHOUT LURES ###
#############################################################

### organize and clean data 

setwd("/Users/marianagelambi/Desktop/chemical-lures/objective_1") # set directory to objective 1 folder
data <- read.csv("data.csv")
head(data)
summary(data)

data_2 <- data %>%  # create a new dataframe named "allsites" from the "data" dataframe
  select(-c(cperspicillata, csowelli, ccastanea)) %>%  # remove columns named "cperspicillata", "csowelli", and "ccastanea" from the dataframe
  group_by(site, treatment) %>%  # group the data by "site" and "treatment" columns
  summarize(across(bats:insectivorous_bats, sum)) %>%  # for all columns from "bats" to "insectivorous_bats", calculate the sum of each group
  pivot_longer(cols = -c(site, treatment), names_to = "bat_species", values_to = "values")  # reshape the data into long format, with "site" and "treatment" as ID variables, "bat_species" as a new column containing the previously column names, and "values" as a new column containing the corresponding values

# Write output to file
write.csv(data_2, file = "data_nodates.csv", row.names = FALSE)



