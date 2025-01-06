# Load Libraries

library(tidyverse)
library(data.table)
library(textshape)
library(gplots)


# For Autophagy Genes

c <- read.csv("data/Autophagy_ATGs.csv", header = T, stringsAsFactors = F)
d <- read.csv("data/counts_genes.csv", header = T, stringsAsFactors = F)
#e <- read.csv("data/ED_APR18_0.05_LFC1.5_190502.csv", header = T, stringsAsFactors = F)

#MapBIN_SUBA_LFQ <- left_join(e, c, by = "AGI") %>% left_join(d, by = "AGI")

#MapBIN_SUBA <- dplyr::full_join(c,d, by = "AGI")

#MapBIN_SUBA_LFQ <- dplyr::full_join(MapBIN_SUBA,e, by = "AGI")

Selected_Count_Genes <- left_join(c,d, by = "AGI")
write.csv(Selected_Count_Genes, "ATG_Selected_Count_Genes.csv")


library(stringr)
library(tidyr)
library(dplyr)

df2 <- read.csv("data/ATG_Selected_Count_Genes_ED.csv", header = T, stringsAsFactors = F)

gdf2 <-  gather(df2, "group", "Expression", -AGI)

gdf2$tgroup <-  apply(str_split_fixed(gdf2$group, "_", 3)[, c(1, 2)], 1, paste, collapse ="_")

AVG_Count_Genes <- gdf2 %>% group_by(AGI, tgroup) %>% summarize(expression_mean = mean(Expression)) %>% spread(., tgroup, expression_mean)

write.csv(AVG_Count_Genes, "ATG_AVG_Count_Genes.csv")




# Plotting

library(tidyverse)

df <- read.csv("data/ATG_Selected_Count_Genes_ED.csv", header = T, stringsAsFactors = F)

longdata <- df %>% 
  gather(key, value, -AGI) %>% 
  separate(key, c("Genotype", "Time", "Replicate"))

library(plotrix)
# https://stackoverflow.com/questions/29821841/dplyr-summarise-each-standard-error-function

longdata2 <- longdata %>% 
  group_by(AGI, Genotype, Time) %>% 
  summarise_each(funs(mean,sd,std.error)) %>% 
  ungroup() %>% 
  mutate(Time = factor(Time), Genotype = factor(Genotype))

colnames(longdata2)

# Sorting bars by factor ordering
#https://sebastiansauer.github.io/ordering-bars/

longdata2$Genotype <- factor(longdata2$Genotype,levels = c("WT", "atg2", "atg5", "atg9"))

# Plotting

longdata2 %>% 
  ggplot(aes(x = Time, y = value_mean, fill = Genotype)) + 
  geom_bar(position = position_dodge(), stat = "identity") + 
  geom_errorbar(aes(ymin = value_mean - value_std.error, ymax = value_mean + value_std.error), 
                width = 0.1, 
                position = position_dodge(0.9)) + 
  scale_fill_manual(values=c("#008000","#B8860B","#4169E1","#DC143C"))+
  facet_wrap( ~ AGI,scales = "free_y")


## HeatMap
