# Load Libraries
library(tidyverse)
library(data.table)
library(ggplot2)

##############################
##### Data preparation #######
##############################

#read data
df <- fread("data/ATG_AVG_Count_Genes_ED.csv", header = T, stringsAsFactors = F)

df2 <- mutate(df, atg2_0h_vs_WT_0h = atg2_0h/WT_0h,
                  atg2_12h_vs_WT_12h = atg2_12h/WT_12h,
                  atg2_24h_vs_WT_24h = atg2_24h/WT_24h,
                  atg5_0h_vs_WT_0h = atg5_0h/WT_0h,
                  atg5_12h_vs_WT_12h = atg5_12h/WT_12h,
                  atg5_24h_vs_WT_24h = atg5_24h/WT_24h,
                  atg9_0h_vs_WT_0h = atg9_0h/WT_0h,
                  atg9_12h_vs_WT_12h = atg9_12h/WT_12h,
                  atg9_24h_vs_WT_24h = atg9_24h/WT_24h)

write.csv(df2, "ATG_HM_VSWT.csv")

################
### HeatMap ###
###############

df3 <- fread("data/ATG_HM_VSWT_ED.csv")
df4 <- df3 %>% remove_rownames %>% column_to_rownames(var="AGI")

#Working Heatmap for col clustering

cn=colnames(df4)[c(1:9)]
cn
col <- colorRampPalette(c("darkgreen","yellow","red"))(30)
# https://www.worldfullofdata.com/three-ways-create-heatmap-r/
# http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html
# https://www.biostars.org/p/398548/#398696

library(gplots)

row_clust <- hclust(dist(as.matrix(df4[,cn]), method = 'euclidean'), method = 'ward.D2')


par(oma=c(1,1,1,1)); heatmap.2(as.matrix(df4[,cn]), 
                               dendrogram = "row", 
                               Colv = FALSE, 
                               Rowv = as.dendrogram(row_clust),
                               scale = "non", # data scaling ("none","row", "column")
                               col = col,
                               #labRow= NA, # remove row names from plot
                               key = TRUE,
                               keysize=0.75,
                               key.par=list(mar=c(0,0,1,2)),
                               density.info = "none", 
                               key.title = NA, 
                               key.xlab = "Abundance",
                               trace = "none",
                               margins = c(7, 10))
