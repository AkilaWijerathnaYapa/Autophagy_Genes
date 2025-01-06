library(tidyverse)

df <- read_csv("data/Selected_Count_Genes_20200215.csv")

df1 <- gather(df,"TP","Values",-AGI) # TP=Sample

library(stringr)

df2 <- cbind(df1,str_split_fixed(df1$TP,"_",3))
colnames(df2)[4:6] <- c("Genotype","time","replicate")
df2$Trt <-  paste(df2$Genotype, df2$time, sep="_")

library(dplyr)
df3 <- select(df2, -c(replicate))

library(Rmisc)
names(df3)
df4 <- summarySE(df3, measurevar="Values", groupvars=c("time","Genotype","Trt","AGI"))
View(df4)
View(df)

library(ggplot2)
ggplot(df4, aes(time, Values, group = Trt, color = Genotype)) +
  geom_line() +
  geom_point() +
  facet_wrap( ~ AGI,scales = "free_y") +
  labs(title = "Gene expression ", x = "Time (hr)", y = "Measurement") +
  theme_linedraw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    strip.text = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  ) +
  geom_ribbon(aes(ymax = Values + sd, ymin = Values - sd),
              alpha = 0.3,
              fill = "grey70",
              colour=NA
  )



ggplot(df4, aes(time, Values, group = Trt, color = Genotype)) +
  geom_bar() +
  geom_point() +
  facet_wrap( ~ AGI,scales = "free_y") 

ggplot(df4, aes(time, Values, group = Trt, color = Genotype)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Values - se, ymax = Values + se)) +
  facet_wrap( ~ AGI,scales = "free_y") 