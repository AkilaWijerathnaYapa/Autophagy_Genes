library(tidyverse)
library(data.table)
library(broom)
library(multcompView)

#read data
batch1 <- fread("data/ATG_Selected_Count_Genes_ED.csv", header = T, stringsAsFactors = F)
#long format (gather) (tidyr)
batch2 <- batch1 %>% 
  gather(sample,count,2:ncol(.))

#split sample name into genotype, timepoit and replicate
batch3 <- batch2 %>% 
  mutate(genotype=sapply(strsplit(sample, '_'),'[',1),
         time_point=sapply(strsplit(sample, '_'),'[',2)) 


#group data and calculate arevage and se
batch4 <- batch3 %>% 
  group_by(`AGI`,genotype,time_point) %>% 
  mutate(average_measure=mean(count),n=n(),se=sd(count)/sqrt(n)) 


# # T-test
# #S1 <- filter(batch3, transition %in% c('S1'))
# #all_ttest <- data.frame()
# for (trans in c('S1')){
#   for (time in unique(S1$time_point)){
#     for (AGI1 in unique(S1$AGI)){
#       selection <- filter(S1, AGI==AGI1 & time_point == time)
#       selection_ttest <- pairwise.t.test(selection$count,selection$genotype, p.adjust.method = 'BH')
#       selection_ttest <- as.data.frame(selection_ttest$p.value) %>% 
#         rownames_to_column %>% 
#         gather(comp,pvalue, 2:ncol(.)) %>% 
#         na.omit() %>% 
#         arrange(rowname) %>% 
#         mutate(transiton=trans, time_point=time, AGI=AGI1)
#       all_ttest <- bind_rows(all_ttest,selection_ttest)
#     }
#   }
# }

# T-test### Re-do with adjsting
S1 <- filter(batch3)
all_ttest <- data.frame()
  for (time in unique(S1$time_point)){
    for (AGI1 in unique(S1$AGI)){
      selection <- filter(S1, AGI==AGI1 & time_point == time)
      selection_ttest <- pairwise.t.test(selection$count,selection$genotype, p.adjust.method = 'none', pool.sd=F)
      selection_ttest <- as.data.frame(selection_ttest$p.value) %>% 
        rownames_to_column %>% 
        gather(comp,pvalue, 2:ncol(.)) %>% 
        na.omit() %>% 
        arrange(rowname) %>% 
        mutate(time_point=time, AGI=AGI1)
      all_ttest <- bind_rows(all_ttest,selection_ttest)
    }
  }


write.csv(all_ttest, "ATG_S2_all_ttest.csv")

#Letter 
all_AGIs <- all_ttest %>% mutate(name=paste0(rowname,'-',comp), name=gsub(' ','',name))
all_letters <- data.frame()
for (t in unique(S1$AGI)){
  for(time in unique(all_AGIs$time_point)){
    AGI_frame <- filter(all_AGIs, AGI==t & time_point==time)
    AGI_pvalues <- AGI_frame$pvalue
    names(AGI_pvalues) <- AGI_frame$name
    AGI_letters <- multcompLetters(AGI_pvalues)
    AGI_letters <- as.data.frame(AGI_letters$Letters) %>% rownames_to_column(var='genotype')
    colnames(AGI_letters)[2] <- 'letter'
    AGI_letters <- AGI_letters %>% mutate(AGI=t, time_point=time)
    all_letters <- bind_rows(all_letters,AGI_letters)
  }
}

#plot

#remove single counts from batch4
batch5 <- batch4 %>% distinct(time_point,AGI,gentype, .keep_all = T) %>% left_join(all_letters)

#relevel data and pick colours
batch5$genotype <- factor(batch5$genotype, levels = c("WT","atg2","atg5","atg9"))
colours <- c('#336600','#ffcc00','#3366cc','#990000')

##Plot faceted
library(ggplot2)
ggplot(filter(batch5),aes( x=time_point, y=average_measure, fill=genotype))+
  geom_bar(stat='identity', position='dodge')+scale_fill_manual(values=colours)+
  geom_errorbar(aes(ymin=average_measure-se, ymax=average_measure+se),position=position_dodge(width=0.89),width=0.3,size=0.5)+
  geom_text(aes(y=average_measure+se+.2,label=letter),position=position_dodge(width=.9),size=3)+
  facet_wrap(~`AGI`,scales = "free")+
  labs(title='SII Transition', x='Time', y='average peak area [counts]')#+theme_classic()

#Single Plots

library(ggplot2)
for (prot in unique(batch4$`AGI`)){
  p <- ggplot(filter(batch5,`AGI`==prot),aes( x=time_point, y=average_measure, fill=genotype))+
    geom_bar(stat='identity', position='dodge')+scale_fill_manual(values=colours)+
    geom_errorbar(aes(ymin=average_measure-se, ymax=average_measure+se),position=position_dodge(width=0.89),width=0.3,size=0.8)+
    geom_text(aes(y=average_measure+se+.2,label=letter),position=position_dodge(width=.9),size=3)+
    facet_wrap(~`AGI`)+
    labs(title=paste0('Transition-1',' ',prot), x='Time [d]', y='average peak area [counts]' )+
    theme(strip.text.x =element_text(face='bold', size=20), axis.title = element_text( size=16),
          legend.text=element_text(size=16),legend.title=element_text(size=20), axis.text = element_text(size=16))+theme_classic()
  ggsave(paste0('Plots/single_plot_', prot,'.png'),device = 'png',plot = p, width = 16 ,height = 9,units = 'in', dpi = 600)
}


#with data points
batch3$genotype <- factor(batch3$genotype, levels = c("WT","atg2","atg5","atg7","atg9"))

batch5$genotype <- factor(batch5$genotype, levels = c("WT","atg2","atg5","atg7","atg9"))
colours <- c('#336600','#ffcc00','#3366cc','#ff9900','#990000')

library(ggplot2)
ggplot(filter(batch5,transition=='S1'),aes( x=time_point, y=average_measure, fill=genotype))+
  geom_bar(stat='identity', position='dodge')+scale_fill_manual(values=colours)+
  geom_jitter(data=filter(batch3),aes(x=time_point,y=count),size=1,position = position_dodge(width = 0.9))+
  geom_errorbar(aes(ymin=average_measure-se, ymax=average_measure+se),position=position_dodge(width=0.89),width=0.3,size=0.5)+
  geom_text(aes(y=average_measure+se+.2,label=letter),position=position_dodge(width=.9),size=3)+
  facet_wrap(~`AGI`)+
  labs(title='SI Transition', x='Time', y='average peak area [counts]')


