library(data.table)
library(GenomicRanges)
library(dplyr)
library(InteractionSet)
library(GenomicInteractions)
library(genomation)

library(ggplot2)
library(PerformanceAnalytics)
library(MASS)

setwd('~/SHARE/HYQ_atac/analysis/co_accessibility/')

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# select pairs of co-accessible peaks, 
# make interaction objects, 
# and write them in bedpe files.
if (F) {
  for (ct in c("B","CD4_T","CD8_T","eHAP1","GM12878","HEK293T","HFF1","K562","Monocyte")){
    peaks <-
      readBed(paste0('./data/peak/nanoATAC/',ct,'.bed')) %>% 
      sort
    
    gi.clean<-NULL
    coa.all<-NULL
    
    for (chrom in c(1:22,'X')){
      # get subset of peaks by chromosome
      peaks.chr <- subset(peaks,seqnames==paste0('chr',chrom))
      
      peak.test<-readRDS(paste0('./data/ks_test/nanoATAC.',ct,'.chr',
                                chrom,'.co-accessibility.ks.Rds'))

      # extract co-accessible peak pairs
      criteria <- peak.test$left.p.val.adj[1:(nrow(peak.test) - 1)] < 0.05 &
        peak.test$right.p.val.adj[2:nrow(peak.test)] < 0.05
      sig.left.index <- which(criteria)
      
      # generate a GenomicInteraction object
      gi <- GInteractions(sig.left.index, sig.left.index + 1, peaks.chr)
      
      coa <- data.frame(
        span = pairdist(gi),
        left.end = peak.test[sig.left.index,]$left.length,
        right.end = peak.test[sig.left.index + 1,]$right.length,
        left.control = peak.test[sig.left.index, ]$left.control,
        right.control = peak.test[sig.left.index + 1, ]$right.control
      ) %>%
        mutate(
          left.spring = left.end / left.control,
          right.spring = right.end / right.control,
          res = span / left.end
        )
      
      if (is.null(gi.clean)){
        gi.clean<-gi[coa$res<2]
      }else{
        gi.clean<-c(gi.clean,gi[coa$res<2])
      }
      
      if (is.null(coa.all)){
        coa.all<-coa
      }else{
        coa.all<-rbind(coa.all,coa)
      }
    }
  
    export.bedpe(gi.clean,
      paste0('./data/pairs/',ct,'_co-accsessibility.bedpe'))
    
    # density plot of span ~ linkage
    {
      set.seed(1)
      
      df <- data.frame(
        log10.span = log10(coa.all$span),
        log10.left.end = log10(coa.all$left.end),
        log10.res = log10(coa.all$res)
      ) 
      
      df$density <-
        get_density(df$log10.span,
                    df$log10.left.end,
                    n = 100)
      
      pdf(paste0('./result/',ct,'.span_relationship.pdf'),
          width = 7, height = 7)
      
      # color by density
      p <- ggplot(df,aes(x=log10.span, y=log10.left.end, color = density)) +
        geom_point() +
        stat_density2d(h = c(0.5,0.05), color='grey90') +
        geom_hline(yintercept = median(log10(coa.all$right.control))) +
        geom_abline(slope = 1,intercept = -log10(2))+
        scale_color_viridis()+
        coord_fixed(ratio=diff(range(df$log10.span))/diff(range(df$log10.left.end)))+
        theme_bw()
      
      plot(p)
    
      df <-
        rbind(
          data.frame(dist = peaks %>% start %>% diff, group = 'nanoATAC_peak_distance'),
          data.frame(dist = coa.all$left.end, group = 'median_peak_read_length'),
          data.frame(dist = coa.all$span[coa.all$res<2], group = 'co-accessiblity_span \n (res < 2-fold)')
        )
      
      p <- ggplot(data = df, aes(log10(dist), fill = group)) +
        geom_density(alpha = .5) +
        geom_vline(xintercept = median(log10(coa.all$left.control)))+
        theme_bw()
      
      plot(p)
      
      dev.off()
    }
  }
}

# illustration of co-accessibility detection
# taking SOX4 as an example
if (F) {
  frags <- fread('./data/fragment/GM12878.flanked.bed.gz')
  colnames(frags) <- c('chr', 'start', 'end', 'read', 'mapQ', 'side')
  
  peaks <-
    readBed('./data/peak/nanoATAC/GM12878.bed') %>% 
    sort
  
  chrom <- '6'
  
  # get subset of peaks by chromosome
  peaks.chr <- subset(peaks,seqnames==paste0('chr',chrom))
  
  # get subset of fragments by chromosome and filter supplemntary alignemnts
  frags.chr <- subset(frags, chr==chrom)
  tb <- table(frags.chr$read)
  frags.chr <- frags.chr %>% subset(! read %in% names(tb)[tb>2])
  frags.chr <- frags.chr[order(frags.chr$read,frags.chr$chr,frags.chr$start),]
  rm(tb);rm(frags);gc()
  
  # calculate read length
  rl.chr <- subset(frags.chr, side == '+')$start - 
      subset(frags.chr, side == '-')$start
  names(rl.chr)  <- frags.chr$read %>% unique

  pdf('./result/220626_GM12878_SOX4_co-accessibility_illustration.pdf',
      width = 6,height = 4)

  i <- which(start(peaks.chr) > 21588918 & end(peaks.chr) < 21589991)
  
  result <- list()
  pval <- c()
  for (p in 0:1) {
    peak <- resize(peaks.chr[i+p], width = 1e3, fix = 'center')
    control <- resize(peaks.chr[i+p], width = 1e5, fix = 'center')
    fr <- subset(frags.chr, start >= start(peak) & end <= end(peak))
    fr.control <-
      subset(frags.chr, start >= start(control) & end <= end(control))
    
    rl.chr.peak <- rl.chr[fr$read]
    rl.chr.control <- rl.chr[fr.control$read]
    result[[p+1]] <- list(peak = rl.chr.peak, control = rl.chr.control)
    
    pval[1+p] <- ks.test(rl.chr.peak, rl.chr.control)$p.value
  }
  
  left<-peaks.chr[i] %>% as.character()
  right<-peaks.chr[i+1] %>% as.character()
  
  df <-
    rbind(
      data.frame(length = result[[1]]$peak,
                 type = 'Peak', group=left),
      data.frame(length = result[[1]]$control,
                 type = 'Control\n(100 Kb)', group=left),
      data.frame(length = result[[2]]$peak,
                 type = 'Peak', group=right),
      data.frame(length = result[[2]]$control,
                 type = 'Control\n(100 Kb)', group=right)
    )
  
  annotations <- data.frame(
    xpos = c(Inf, Inf),
    ypos =  c(Inf, Inf),
    p.value = paste0('KolmogorovSmirnov Test\nP-value: ',signif(pval,3)),
    group = c(left, right),
    hjustvar = c(1.1, 1.1) ,
    vjustvar = c(1.4, 1.4)
  )
  
  ggplot() +
    facet_grid(row = vars(group)) +
    geom_density(data = df, aes(log10(length), fill = type), alpha = .5) +
    theme_bw() +
    geom_text(data = annotations,
              aes(
                xpos,
                ypos,
                hjust = hjustvar,
                vjust = vjustvar,
                label = p.value
              )) +
    ylab('density')
  
  dev.off()
}
