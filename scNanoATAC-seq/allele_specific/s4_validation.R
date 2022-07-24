# HYQ long-read atac-seq project
#### validate aymmetric accessibility of long atac-seq ####

library(dplyr)
library(data.table)
library(GenomicRanges)
library(parallel)
library(vcfR)

setwd("~/SHARE/HYQ_atac/analysis/asymmetric_haplotype/")

# note: A/B is not pat/total, but alt/total ("1"/total)
if (F) {
  vcf <- read.vcfR("./ngs_snp/NGS-GE-GM12878-X1.giab_hg001_het.vcf.gz")
  dp <- extract.info(vcf,element = 'DP', as.numeric = T)
  ad <- extract.gt(vcf,element = 'AD')
  ratio <- strsplit(ad[,1],split = ',') %>% 
    sapply(function(x) {
      if (!is.na(x)){
        a<-as.numeric(x[1])
        b<-as.numeric(x[2])
        a/(a+b)
      }else{
        NA
      }
    })
  
  # generate a GRange object for ATAC-seq SNP
  vcf_ref <- read.vcfR('./ngs_snp/ref.hg001.het.vcf.gz')
  gt_ref <-
    extract.gt(vcf_ref)[paste0(getCHROM(vcf), "_", getPOS(vcf)), ]
  
  snp <-
    GRanges(
      seqnames = getCHROM(vcf),
      ranges = IRanges(start = getPOS(vcf), width = 1),
      dp = dp,
      ratio = ratio,
      gt_ref = gt_ref
    )
  seqlevelsStyle(snp) <- 'UCSC'
  snp <- snp[!is.na(snp$ratio)]
  snp <- snp[snp$dp > 15]
  
  # visualize
  hist(snp$ratio, breaks = 100)
  plot(snp$dp, snp$ratio, cex = .1)
}

# load long atac-seq biased peak call set
if (F) {
  long.asym <- readRDS('220625_GM12878_TGS_peaks_biased.hapcut2_phased.Rds')
  
  df <-
    fread('./peak/GM12878.bed',
          header = F,
          sep = '\t') %>% as.data.frame()
  names(df) <-
    c('chr', 'start', 'end')
  peak.gr <- GRanges(seqnames = df$chr,
                     ranges = IRanges(df$start, df$end))
  seqlevelsStyle(peak.gr) <- 'UCSC'
}

# validation of bias by NGS ATAC-seq SNP ratio
if (F) {
  # maternal or paternal specific asymmetric peaks
  mode <- 'mat'
  #mode <- 'pat'
  
  long.cut.off <- 0.05
  ngs.cut.off <- 0.05
  
  if (mode == 'mat') {
    sig <- long.asym$fdr < long.cut.off & long.asym$ratio < .5
    peak.index <-
      long.asym$peak[sig] %>% unique
    minor.gt <- '1|0'
    
    # the count of autosome ASPs
    table(long.asym[sig,]$chr=='chrX')
  }
  
  if (mode == 'pat') {
    sig <- long.asym$fdr < long.cut.off & long.asym$ratio > .5
    peak.index <-
      long.asym$peak[sig] %>% unique
    minor.gt <- '0|1'
    
    # the count of autosome ASPs
    table(long.asym[sig,]$chr=='chrX')
  } 
  
  hit <- findOverlapPairs(peak.gr[peak.index], snp)
  
  binom.p <- c()
  for (i in seq_along(hit@second)) {
    x <- hit@second[i]
    binom.p <- c(binom.p,
                 binom.test(round(x$ratio * x$dp), x$dp, p = 0.5)$p.value)
  }
  binom.p <- p.adjust(binom.p, method = 'fdr')
  
  ##### STATISTICS ####
  # validatable peaks
  hit@first %>% unique() %>% length
  # validation SNPs
  hit@second %>% unique() %>% length
  
  criteria <- (binom.p < ngs.cut.off)
  hit.peak <- hit@first[criteria]
  hit.snp <- hit@second[criteria]
  
  # gross precision (by peak)
  (hit.peak %>% unique %>% length) / (hit@first %>% unique %>% length)
  # gross precision (by SNP, if any SNP in the peak supports)
  (hit.snp %>% unique %>% length) / (hit@second %>% unique %>% length)
  
  # validate consistency of bias direction between long and short atac
  table(hit.snp$gt_ref == minor.gt,
        hit.snp$ratio > .5)
  valid <- !xor(hit.snp$gt_ref == minor.gt,
                hit.snp$ratio > .5)
  # consistent precision in significant subset (by SNP)
  sum(valid) / length(valid)
  
  # overall precision (by peak)
  (hit.peak[valid] %>% unique %>% length) / 
    (hit@first %>% unique %>% length)
  
  # overall precision (by SNP)
  (hit.snp[valid] %>% unique %>% length) / 
    (hit@second %>% unique %>% length)
}

# rebuttal (HapCut2 pahsing)
if (F) {
  long.asym<-readRDS('220625_GM12878_TGS_peaks_biased.hapcut2_phased.Rds')
  
  df <-
    fread('./peak/GM12878.bed',
          header = F,
          sep = '\t') %>% as.data.frame()
  names(df) <-
    c('chr', 'start', 'end')
  peak.gr <- GRanges(seqnames = df$chr,
                     ranges = IRanges(df$start, df$end))
  seqlevelsStyle(peak.gr) <- 'UCSC'
}

# mat/pat is undistinguishable in the de novo phasing mode
# validation of bias by NGS ATAC-seq SNP ratio
if (F) {
  long.cut.off <- 0.05
  ngs.cut.off <- 0.05
  
  peak.index <-
    long.asym$peak[long.asym$fdr < long.cut.off] %>% unique
  
  table(long.asym$chr[long.asym$fdr < long.cut.off]=='chrX')

  hit <- findOverlapPairs(peak.gr[peak.index], snp)
  
  binom.p <- c()
  for (i in seq_along(hit@second)) {
    x <- hit@second[i]
    binom.p <- c(binom.p,
                 binom.test(round(x$ratio * x$dp), x$dp, p = 0.5)$p.value)
  }
  binom.p <- p.adjust(binom.p, method = 'fdr')
  
  ##### STATISTICS ####
  # validatable peaks
  hit@first %>% unique() %>% length
  # validation SNPs
  hit@second %>% unique() %>% length
  
  criteria <- (binom.p < ngs.cut.off)
  hit.peak <- hit@first[criteria]
  hit.snp <- hit@second[criteria]
  
  # gross precision (by peak)
  (hit.peak %>% unique %>% length) / (hit@first %>% unique %>% length)
  # gross precision (by SNP, if any SNP in the peak supports)
  (hit.snp %>% unique %>% length) / (hit@second %>% unique %>% length)
}