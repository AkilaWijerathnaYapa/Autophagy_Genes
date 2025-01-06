library(dplyr)
library(readr)

c <- read.csv("data/bZIP.csv", header = T, stringsAsFactors = F)
d <- read.csv("data/counts_genes.csv", header = T, stringsAsFactors = F)
#e <- read.csv("data/ED_APR18_0.05_LFC1.5_190502.csv", header = T, stringsAsFactors = F)

#MapBIN_SUBA_LFQ <- left_join(e, c, by = "AGI") %>% left_join(d, by = "AGI")

#MapBIN_SUBA <- dplyr::full_join(c,d, by = "AGI")

#MapBIN_SUBA_LFQ <- dplyr::full_join(MapBIN_SUBA,e, by = "AGI")

Selected_Count_bZIP <- left_join(c,d, by = "AGI")
write.csv(Selected_Count_bZIP, "Selected_Count_bZIP.csv")


library(stringr)
library(tidyr)
library(dplyr)

df2 <- read.csv("data/Selected_Count_bZIP.csv", header = T, stringsAsFactors = F)

gdf2 <-  gather(df2, "group", "Expression", -AGI)

gdf2$tgroup <-  apply(str_split_fixed(gdf2$group, "_", 3)[, c(1, 2)], 1, paste, collapse ="_")

AVG_Count_bZIP <- gdf2 %>% group_by(AGI, tgroup) %>% summarize(expression_mean = mean(Expression)) %>% spread(., tgroup, expression_mean)

write.csv(AVG_Count_bZIP, "AVG_Count_bZIP.csv")


# df2<- read.csv("data/AVG_Count_bZIP_Mod.csv", header = T, stringsAsFactors = F)
# 
# 
# # First Column to Raw data
# df3 <- data.frame(df2, row.names = 1)
# 
# 
# df4 <- as.matrix(as.data.frame(df3))
# 
# typeof(df4)
# 
# library(made4)
# 
# heatplot(df4,dend="none", scale="none")
# 
# 
# 


library(tidyverse)

df <- read_csv("data/Selected_Count_bZIP2.csv")

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


# #Subsetting Genotypic Sig data
# 
# longdata3 <- subset(longdata2, AGI %in% c("AT5G53250",
#                                           "AT5G65390",
#                                           "AT3G17160",
#                                           "AT1G58380",
#                                           "AT3G07390",
#                                           "AT5G47700",
#                                           "AT4G00810",
#                                           "AT5G13580",
#                                           "AT2G05765",
#                                           "AT3G14440",
#                                           "AT3G51730",
#                                           "AT5G66930",
#                                           "AT1G10600",
#                                           "AT1G06650",
#                                           "AT4G04620"))
# 
# longdata3$Genotype <- factor(longdata3$Genotype,levels = c("WT", "atg2", "atg5", "atg9"))
# 
# # Plotting
# 
# longdata3 %>% 
#   ggplot(aes(x = Time, y = value_mean, fill = Genotype)) + 
#   geom_bar(position = position_dodge(), stat = "identity") + 
#   geom_errorbar(aes(ymin = value_mean - value_std.error, ymax = value_mean + value_std.error), 
#                 width = 0.1, 
#                 position = position_dodge(0.9)) + 
#   scale_fill_manual(values=c("#008000","#B8860B","#4169E1","#DC143C"))+
#   xlab("Gene Count / TPMs")+
#   ylab("Time (hours)")+
#   facet_wrap( ~ AGI,scales = "free_y")
# 
# 
# 
# #Subsetting Time_Part 1 data
# 
# longdata4 <- subset(longdata2, AGI %in% c("AT3G22840",
#                                           "AT5G54190",
#                                           "AT2G38530",
#                                           "AT5G02270",
#                                           "AT4G14130",
#                                           "AT1G65310",
#                                           "AT2G46790",
#                                           "AT1G53830",
#                                           "AT5G60660",
#                                           "AT5G45690",
#                                           "AT4G02290",
#                                           "AT4G12520",
#                                           "AT5G53250",
#                                           "AT1G22770",
#                                           "AT1G75780",
#                                           "AT5G44440",
#                                           "AT5G44110",
#                                           "AT5G11260"))
# 
# longdata4$Genotype <- factor(longdata4$Genotype,levels = c("WT", "atg2", "atg5", "atg9"))
# 
# # Plotting
# 
# longdata4 %>% 
#   ggplot(aes(x = Time, y = value_mean, fill = Genotype)) + 
#   geom_bar(position = position_dodge(), stat = "identity") + 
#   geom_errorbar(aes(ymin = value_mean - value_std.error, ymax = value_mean + value_std.error), 
#                 width = 0.1, 
#                 position = position_dodge(0.9)) + 
#   scale_fill_manual(values=c("#008000","#B8860B","#4169E1","#DC143C"))+
#   xlab("Gene Count / TPMs")+
#   ylab("Time (hours)")+
#   facet_wrap( ~ AGI,scales = "free_y")
# 
# 
# #Subsetting Time_Part 2 data
# 
# longdata5 <- subset(longdata2, AGI %in% c("AT5G15230",
#                                           "AT1G64110",
#                                           "AT2G25160",
#                                           "AT5G46900",
#                                           "AT1G78440",
#                                           "AT5G37300",
#                                           "AT1G05870",
#                                           "AT5G17050",
#                                           "AT3G22830",
#                                           "AT1G20020",
#                                           "AT5G38430",
#                                           "AT1G76080",
#                                           "AT1G55670",
#                                           "AT5G38410",
#                                           "AT5G48485",
#                                           "AT1G70760",
#                                           "AT1G72600",
#                                           "AT3G21055"))
# 
# longdata5$Genotype <- factor(longdata5$Genotype,levels = c("WT", "atg2", "atg5", "atg9"))
# 
# # Plotting
# 
# longdata5 %>% 
#   ggplot(aes(x = Time, y = value_mean, fill = Genotype)) + 
#   geom_bar(position = position_dodge(), stat = "identity") + 
#   geom_errorbar(aes(ymin = value_mean - value_std.error, ymax = value_mean + value_std.error), 
#                 width = 0.1, 
#                 position = position_dodge(0.9)) + 
#   scale_fill_manual(values=c("#008000","#B8860B","#4169E1","#DC143C"))+
#   xlab("Gene Count / TPMs")+
#   ylab("Time (hours)")+
#   facet_wrap( ~ AGI,scales = "free_y")
