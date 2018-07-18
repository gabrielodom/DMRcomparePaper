# Script to Create the Initial Data File
# Zhen Gao
# Gabriel Odom (editor)
# 2018

# source("http://bioconductor.org/biocLite.R")
# biocLite("GEOquery")

library(GEOquery)
library(dplyr)

# set the location to hold gse data
destFolder_char <- "~/2T_Disk/Dropbox (BBSR)/Zhen_Gao/DMR-Saurav/Auto-downloader"
setwd(destFolder_char)

# Get GEO data
gse <- getGEO("GSE41169", GSEMatrix = TRUE, destdir = destFolder_char)
show(gse)

# Subject meta-data
people_meta_data <- pData(gse[[1]])
people_meta_data$GSM <- rownames(people_meta_data)

# Define Data Subset
healthy_male_control_in_20 <-
  people_meta_data %>%
  # select healthy people
  filter(characteristics_ch1.7 == "diseasestatus (1=control, 2=scz patient): 1") %>%
  # select male
  filter(`gender:ch1` == "Male") %>%
  # change data format from "age: 25" to "25"
  mutate(characteristics_ch1.6 = gsub("age: ", "", characteristics_ch1.6)) %>%
  # select larger or equal than 20
  filter(characteristics_ch1.6 >= 20) %>%
  # select less than 30
  filter(characteristics_ch1.6 < 30)

# Select Data and Save
selected_people <-
  healthy_male_control_in_20[, c("GSM",
                                 "gender:ch1",
                                 "characteristics_ch1.6",
                                 "characteristics_ch1.7")]
write.csv(selected_people, "healthy_male_control_in_20.csv")

# Get Expression Data and Save
exp <- exprs(gse[[1]])
data_of_healthy_male_control_in_20 <- exp[, healthy_male_control_in_20$GSM]
write.csv(
  data_of_healthy_male_control_in_20,
  "450K_of_healthy_male_control_in_20.csv"
)

