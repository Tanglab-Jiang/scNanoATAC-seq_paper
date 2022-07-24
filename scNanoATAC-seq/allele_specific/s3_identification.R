# HYQ long-read atac-seq project

library(dplyr)
library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)

setwd("~/SHARE/HYQ_atac/analysis/asymmetric_haplotype/")

#### functions ####
if (T) {
  # 211208
  # function for loading bed coordinates into GRange objects
  readBedAsGrange <- function(x) {
    df <-
      fread(x,
            header = F,
            sep = '\t') %>% as.data.frame() %>% .[, 1:3] 
    names(df) <-
      c('chr', 'start', 'end')
    df<-df[!grepl('KI|GL|M',df$chr),]
    
    x <- GRanges(seqnames = df$chr,
                 ranges = IRanges(df$start, df$end))
    seqlevelsStyle(x) <- 'UCSC'
    
    x
  }

  # 210926
  # write Grange object to bed file
  writeGrangeAsBed <- function(gr,name=NULL){
    as.data.frame(gr)[,1:3] %>%
      write.table(
        file = name,
        quote = F,
        col.names = F,
        row.names = F,
        sep = '\t'
      )
  }
  
  # unphased reads are filterd in the bam2bed step
  loadFrags<-function(frag.file){
    df <-
      fread(frag.file,
            header = F,
            sep = '\t') %>% as.data.frame()
    names(df) <-
      c('chr', 'start', 'end', 'read', 'cell', 'side', 'contig', 'phase')
    frag <- GRanges(
      seqnames = df$chr,
      ranges = IRanges(df$start, df$end),
      read = df$read,
      cell = df$cell,
      side = df$side,
      contig = df$contig,
      phase = df$phase
    )
    seqlevelsStyle(frag) <- 'UCSC'
    frag
  }
  
  # find overlaps and do hypothesis test
  # only for diploid samples
  # fix 211221: change the cut-off from total frags > 15 to phased frags > 15
  testHaplotypeAsymmetric <- function(peak,frag,mc=NULL) {
    frag.cut.off=15
    
    worker <- function(i) {
      result <- list()
      isec.frag <- frag[unique(hits@from[hits@to == i])]
      
      for (ctg in unique(isec.frag$contig)) {
        isec.ctg.frag <- isec.frag[isec.frag$contig == ctg]
        count.ps1 <-
          isec.ctg.frag$cell[isec.ctg.frag$phase == '1'] %>%
          unique %>%
          length
        count.ps2 <-
          isec.ctg.frag$cell[isec.ctg.frag$phase == '2'] %>%
          unique %>%
          length
        
        if (count.ps1 + count.ps2 >= frag.cut.off){
          binom.p <-
            binom.test(count.ps1, count.ps1 + count.ps2, p = 0.5)$p.value
          result$peak <- c(result$peak, i)
          result$contig <- c(result$contig, ctg)
          result$p.value <- c(result$p.value, binom.p)
          result$ratio <-
            c(result$ratio, count.ps1 / (count.ps1 + count.ps2))
          result$count <- c(result$count, count.ps1 + count.ps2)
        }
      }
      result
    }
    
    hits <- findOverlaps(frag, peak)
    result <- mclapply(unique(hits@to), worker, mc.cores = mc)
    
    result.df <- data.frame(
      peak = result %>% sapply(function(x)
        x$peak) %>% unlist,
      contig = result %>% sapply(function(x)
        x$contig) %>% unlist,
      p.value = result %>% sapply(function(x)
        x$p.value) %>% unlist,
      ratio = result %>% sapply(function(x)
        x$ratio) %>% unlist,
      count = result %>% sapply(function(x)
        x$count) %>% unlist
    )
    result.df$fdr <- result.df$p.value %>% p.adjust(method = 'fdr')
    
    result.df<-result.df[order(result.df$count, decreasing = T), ]
    result.df$chr <-peak[result.df$peak] %>% seqnames %>% as.character
    result.df
  }
}

#### resources ####
if (T) {
  ccres<-readRDS('../../resource/210926_ENCODE_cCREs_GRanges.Rds')
}

#### calling ####
if (F) {
  df <-
    fread('./peak/GM12878.bed',
          header = F,
          sep = '\t') %>% as.data.frame()
  names(df) <-
    c('chr', 'start', 'end')
  peak <- GRanges(seqnames = df$chr,
                  ranges = IRanges(df$start, df$end))
  seqlevelsStyle(peak) <- 'UCSC'
  
  # GIAB Phased
  if (T) {
    frags <-
      loadFrags('./fragment/GM12878.giab_benchmark.whatshap_tag.flank.bed.gz')
    result <- testHaplotypeAsymmetric(peak, frags, mc = 3)
    saveRDS(result, '220625_GM12878_TGS_peaks_biased.Rds')
  }
  
  # rebuttal: hapcut2 de novo phasing from prior GIAB heterozygous SNP
  if (F) {
    frags <-
      loadFrags('./fragment/GM12878.hapcut2_phased.whatshap_tag.flank.bed.gz')
    result <- testHaplotypeAsymmetric(peak, frags, mc = 3)
    saveRDS(result,
            '220625_GM12878_TGS_peaks_biased.hapcut2_phased.Rds')
  }
}

# write peak files
if (F) {
  peak.index <- result$peak[result$fdr < cut.off] %>% unique
  
  peak[peak.index] %>% 
    writeGrangeAsBed('./biased_peak/GM12878.asym.giab.fdr_5e-2.count_15.bed')
  
  peak.index <- result$peak[result$fdr < cut.off & result$ratio < .5] %>% unique
  peak[peak.index] %>% 
    writeGrangeAsBed('./biased_peak/GM12878.asym.giab.fdr_5e-2.count_15.maternal.bed')
  
  peak.index <- result$peak[result$fdr < cut.off & result$ratio > .5] %>% unique
  peak[peak.index] %>% 
    writeGrangeAsBed('./biased_peak/GM12878.asym.giab.fdr_5e-2.count_15.paternal.bed')
}

# visualize
if (F) {
  cut.off <- 5e-2
  
  pdf('220625_GM12878_TGS_peaks_fdr_5e-2.pdf',width = 5,height = 5)
  
  result.pl<-result[order(result$fdr,decreasing = T),]
  
  # color paternal and maternal specific peaks
  col <- c('grey70','deepskyblue2', 'firebrick1')[
    (result.pl$fdr < cut.off) + (result.pl$fdr < cut.off & result.pl$ratio < .5)  + 1
  ] %>% adjustcolor(alpha.f = .5)
  
  plot(
    result.pl$count,
    result.pl$ratio,
    pch = 20,
    col = col,
    main = 'Allele Ratio ~ Fragment Count',
    xlab = 'Count',
    ylab = 'Ratio',
    xlim = c(0, max(result.pl$count))
  )
  
  # highlight chrX
  # put chrX points on the upper layer
  result.pl<-rbind(subset(result.pl,chr!='chrX'),subset(result.pl,chr=='chrX'))
  sig.chrx<-(result.pl$fdr < cut.off & result.pl$chr == 'chrX')
  col <- c('grey70', 'deepskyblue2', 'firebrick1')[
    sig.chrx + (sig.chrx & result.pl$ratio < .5) + 1] %>%
    adjustcolor(alpha.f = .5)
  cex <- c(1,1.5)[(result.pl$fdr < cut.off & result.pl$chr=='chrX')+1]
  
  plot(
    result.pl$count,
    result.pl$ratio,
    cex = cex,
    col = col,
    pch = 20,
    xlab = 'Count',
    ylab = 'Ratio',
    xlim = c(0, max(result.pl$count)),
    main = 'Allele Ratio ~ Fragment Count (chrX only)'
  )
  
  # rank plot
  result.pl<-result[order(result$fdr,decreasing = T),]
  sig<-(result.pl$fdr < cut.off)
  col <- c('grey30','deepskyblue2', 'firebrick1')[
    sig + (sig & result.pl$ratio < .5)  + 1
    ] %>% adjustcolor(alpha.f = .15)
  
  plot(
    result.pl$ratio,
    rank(result.pl$ratio),
    col = col,
    pch = 20,
    cex = c(0.1, 2.5)[(result.pl$fdr < cut.off) + 1],
    ylim = c(-1e3, nrow(result.pl) + 1e3),
    xlab = 'Ratio',
    ylab = 'Rank of Ratio',
    main = 'Rank Plot for Allele Ratio'
  )
  
  plot(
    result.pl$ratio,
    rank(result.pl$ratio),
    col = col,
    pch = 20,
    cex = (-log10(result.pl$fdr)) ,
    ylim = c(-1e3, nrow(result.pl) + 1e3),
    xlab = 'Ratio',
    ylab = 'Rank of Ratio',
    main = 'Rank Plot for Allele Ratio'
  )
  
  dev.off()
}

# maternal specific accessibility is dominant on chrX 
# under the condition of using the GIAB benchmark SNP set which is fully phased
if (F) {
  table(
    (result[(peak[result$peak] %>% seqnames %>% as.character())=='chrX',] %>% 
     subset(fdr < cut.off) %>% .$ratio) < .5
  ) %>% binom.test()
  
  table(
    (result[(peak[result$peak] %>% seqnames %>% as.character())!='chrX',] %>% 
     subset(fdr < cut.off) %>% .$ratio) < .5
  ) %>% binom.test()
}

# chrX asymmetric accessibility is stronger than autosomes in XX cells
if (F) {
  pdf(
    '220625_GM12878_biased_ratio_by_chromosome.pdf',
    width = 4,
    height = 4
  )
  
  tb <- table(result$chr, result$fdr < cut.off)
  df <- data.frame(biased.peak.ratio = tb[, 2] / rowSums(tb),
                   chr = rownames(tb))
  df$chr<-factor(df$chr,paste0('chr',c(1:22,'X')))
  df$col<-'black'
  df$col[df$chr=='chrX']<-'red'
  
  ggplot(df, aes(x = chr, y = biased.peak.ratio)) +
    geom_bar(stat = "identity", fill = df$col)  +
    scale_y_continuous(limits=c(0, 0.2), expand = c(0, 0)) +
    xlab('Chromosome') +
    ylab('Biased Peak Ratio') +
    ggtitle('GM12878 Biased Peak Distribution') +
    theme_bw() +
    theme(
      axis.text.x = element_text(color = 'black', angle=90),
      axis.text.y = element_text(color = 'black'),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", size = 1)
    )
  
  dev.off()
}