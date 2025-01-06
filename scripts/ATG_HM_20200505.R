###############################################
## Heat Map for ER ##
###############################################

# For Knn / Imputed Results

# Load Libraries

library(tidyverse)
library(data.table)
library(textshape)
library(gplots)



#####################################
########## Heat Maps ###############
####################################



df1 <- fread("data/ATG_Selected_Count_Genes_ED.csv")

gdf1 <- gather(df1, "group", "Expression",-AGI) %>%
  separate(group, c("sample", "time", "r")) %>%
  unite(tgroup, c("sample", "time"))  %>%
  group_by(AGI, tgroup) %>%
  summarize(expression_mean = mean(Expression)) %>%
  spread(tgroup, expression_mean) %>%
  column_to_rownames(colnames(.)[1])


#Working Heatmap for col clustering

cn=colnames(gdf1)[c(10:12,1:9)]
cn
col <- colorRampPalette(c("darkgreen","yellow","red"))(30)
# https://www.worldfullofdata.com/three-ways-create-heatmap-r/
# http://mannheimiagoesprogramming.blogspot.com/2012/06/drawing-heatmaps-in-r-with-heatmap2.html
# https://www.biostars.org/p/398548/#398696

library(gplots)

row_clust <- hclust(dist(as.matrix(gdf1[,cn]), method = 'euclidean'), method = 'ward.D2')


par(oma=c(1,1,1,1)); heatmap.2(as.matrix(gdf1[,cn]), 
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

# out <- heatmap.2(as.matrix(gdf1[,cn]), 
#                  dendrogram = "row", 
#                  Colv = FALSE, 
#                  Rowv = as.dendrogram(row_clust),
#                  scale = "non", # data scaling ("none","row", "column")
#                  col = col,
#                  labRow=NA, # remove row names from plot
#                  key = TRUE, 
#                  density.info = "none", 
#                  key.title = NA, 
#                  key.xlab = "Abundance",
#                  trace = "none",
#                  margins = c(7, 10))
# 
# plot(out$rowDendrogram)
# rownames(gdf1)[out$rowInd]
# 
# plot(row_clust)

#Cut the dendrogram into groups or specify a height for the cut-off:

#2 groups
#sort(cutree(row_clust, k=2))

#5 groups
#sort(cutree(row_clust, k=5))

# specify a height of 50
# a <- sort(cutree(row_clust, h = 50))
# 
# write.csv(a,"sort50.csv")

# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

# plot(row_clust, cex = 0.5)
# 
# 
# out$rowInd
# 
# write.table(
#   data.frame(gene = rownames(gdf1)[out$rowInd]),
#   'HM_out_20190914.csv',
#   row.names = FALSE,
#   quote = FALSE,
#   sep = ',')


# ATG AGI List data sorting for order, after generating Heat Map

d <- read.csv("data/HMGeneOrder.csv", header = T, stringsAsFactors = F)

#Used this file to filter ER data
e <- read.csv("Data/ATG_Selected_Count_Genes_ED.csv", header = T, stringsAsFactors = F)


MapBIN_SUBA_Knn <- left_join(d, e, by = "AGI")

# e$AGI_ED <- str_sub(e$AGI, end = -3)
# ER_MS1 <- left_join(d,e, by = "AGI_ED")

write.csv(MapBIN_SUBA_Knn, "ATG_Transcript.csv")


# After Editing Columns

########### Two-Way NAOVA ##########



#Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(textshape)
library(gplots)
library(tigerstats)

# SII_MRM / SII(NS)_AVGProt_20190726
df1 <- fread("Data/ATG_Transcript.csv")

gdf1 <- gather(df1, "group", "Expression",-AGI) %>%
  separate(group, c("Genotype","Time", "r")) #%>%
#unite(tgroup, c("sample","Transition", "time"))  #%>%
#unite(tgroup, c("sample","Transition"))  #%>%
# group_by(AGI, tgroup) %>%
# summarize(expression_mean = mean(Expression)) %>%
# spread(tgroup, expression_mean) %>%
# column_to_rownames(colnames(.)[1])

bwplot(Expression ~ Genotype , data=gdf1, main="Protein Expression Lengths by Genotype")

# Two-way ANOVA with interaction effect

#res.aov3 <- aov(Expression ~ Genotype * Time, data = gdf1)

res.aov3 <- aov(Expression ~ Genotype + Time + Genotype:Time, data = gdf1)
summary(res.aov3)


output <- data.frame(AGI = df1$AGI, Genotype  = NA, Time = NA, Interaction = NA)
for (i in 1:nrow(df1)){
  sub_data <- df1[i,]
  gdf1 <- gather(sub_data, "group", "Expression",-AGI) %>%
    separate(group, c("Genotype", "Time", "r"))
  res.aov3 <- summary(aov(Expression ~ Genotype + Time + Genotype:Time, data = gdf1))
  output$Genotype[i] <- res.aov3[[1]]$`Pr(>F)`[1]
  output$Time[i] <- res.aov3[[1]]$`Pr(>F)`[2]
  output$Interaction[i] <- res.aov3[[1]]$`Pr(>F)`[3]
}

write.csv(output, "ATG_ANOVA.csv")


# For Significant Genes

# Box plot with two factor variables
boxplot(Expression ~ Genotype * Time, data=gdf1, frame = FALSE, 
        col = c("#00AFBB", "#E7B800", "#FF0000", "#00FF00", "#0000FF"), ylab="Expression")
# Two-way interaction plot
interaction.plot(x.factor = gdf1$Genotype, trace.factor = gdf1$Time, 
                 response = gdf1$Expression, fun = mean, 
                 type = "b", legend = TRUE, 
                 xlab = "Genotype", ylab="Expression",
                 pch=c(1,19), col = c("#00AFBB", "#E7B800", "#FF0000", "#00FF00", "#0000FF"))

responseList <- gdf1$AGI

modelList    <- lapply(responseList, function(resp) {
  mF <- formula(paste(resp, " ~ AGI"))
  aov(Expression ~ Genotype + Time + Genotype:Time, data = gdf1)
})

lapply(modelList, summary)
