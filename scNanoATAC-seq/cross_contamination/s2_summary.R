library(ggplot2)
library(pheatmap)
library(dplyr)

setwd("~/SHARE/HYQ_atac/analysis/cross_contamination/")

#raw.df<-read.table('./data/2107_ALC1.cross.txt')
raw.df<-read.table('./data/2107_ALC2.cross.txt')

# load cross alignment data
if (T){
  colnames(raw.df)<-c('cell.id','hg38','mm10')
  outer.bc<-raw.df$cell.id %>% 
    strsplit('_') %>% 
    sapply(function(x) x[1])
  
  raw.df$label[outer.bc %in% c(73:76,89:92)] <- 'mouse'
  raw.df$label[outer.bc %in% 77:80] <- 'human'
  raw.df$label[outer.bc %in% 81:88] <- 'mixed'
  
  df<-raw.df[raw.df$hg38+raw.df$mm10>3e3,]
  
  df$label_2<-df$label
  df$label_2[df$label_2 == 'mixed'] <-
    sapply(which(df$label_2 == 'mixed'), function(x) {
      ratio <- df$hg38[x] / (df$mm10[x] + df$hg38[x])
      if (ratio > .9)
        return('mixed - human')
      else if (ratio < .1)
        return('mixed - mouse')
      else
        return('mixed - doublet')
    })
  
  df$label_2<-factor(df$label_2,levels=c('human','mixed - human','mouse','mixed - mouse','mixed - doublet'))
  df<-df[order(df$label_2),]
}

# statistics
if (T){
  # pass QC
  sum(raw.df$hg38+raw.df$mm10>3e3)/nrow(raw.df) * 100
  
  # baseline
  (1-((df$hg38/(df$mm10 + df$hg38))[df$label=='human'] %>% mean)) * 100
  (1-((df$mm10/(df$mm10 + df$hg38))[df$label=='mouse'] %>% mean)) * 100
  
  # missing reads
  (1-((df$hg38/(df$mm10 + df$hg38))[df$label=='mixed'] %>% subset(.>0.9) %>% mean)) * 100
  (1-((df$mm10/(df$mm10 + df$hg38))[df$label=='mixed'] %>% subset(.>0.9) %>% mean)) * 100
  
  # cell error
  miss<-(df$hg38/(df$mm10 + df$hg38))[df$label=='mixed'] %>% subset(.<0.9 & .>0.1) %>% length
  (miss/sum(df$label=='mixed')) * 100
  
  table(df$label)
  table(df$label_2) %>% t %>% t
}

# make scatter plot of cross contamination
if (T){
  pdf('220112_ALC2_cross_contamination.pdf',height = 6,width = 6)
  
  ggplot(df, aes(x = hg38, y = mm10, color = label)) +
    geom_point(size = 2,alpha=c(1,.3)[(df$label=='mixed')+1]) +
    theme_classic() +
    xlab('Reads Mapped to hg38') +
    ylab('Reads Mapped to mm10') +
    coord_fixed(max(df$hg38) / max(df$mm10)) +
    scale_color_manual(values=c('firebrick2','grey80','dodgerblue2'))
  
  ggplot(df, aes(x = hg38, y = mm10, color = label_2)) +
    geom_point(size = 2,alpha=c(1,.3)[(df$label=='mixed')+1]) +
    theme_classic() +
    xlab('Reads Mapped to hg38') +
    ylab('Reads Mapped to mm10') +
    coord_fixed(max(df$hg38) / max(df$mm10))+
    scale_color_manual(values=c('firebrick2','darkorange','dodgerblue2','deepskyblue2','azure4'))
  
  dev.off()
}
