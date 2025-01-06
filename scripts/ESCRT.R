## Libraries

library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)


## Knn LFQ values for selected genes

# Load selected ESCRT AGI List 
# Based on https://www.sciencedirect.com/science/article/pii/S1360138517301760?via%3Dihub#tbl0005

ESCRT<- read.csv("Data/ESCRT.csv", header = T, stringsAsFactors = F)

# Load Knn LFQ 2019APR data
e <- read.csv("Data/HM_19APR_imp_knn_190925.csv", header = T, stringsAsFactors = F)

ESCRT_LFQ_Knn <- left_join(ESCRT, e, by = "AGI")
