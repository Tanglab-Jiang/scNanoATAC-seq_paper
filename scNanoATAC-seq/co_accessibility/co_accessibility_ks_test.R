library(data.table)
library(GenomicRanges)
library(dplyr)
library(genomation)

setwd('~/hyq_atac/tgs/analysis/co_accessibility/')

args <- commandArgs(trailingOnly=T)

# test hypothesis of co-accessibility by chromosome
cell.type <- args[1]
chrom <- args[2]
cores = 5

# nanoATAC fragments
frags <- fread(paste0('../merge_fragment/data/',cell.type,'.flanked.bed.gz'))
colnames(frags) <- c('chr', 'start', 'end', 'read', 'mapQ', 'side')

# nanoATAC peaks
peaks <-
  readBed(paste0('./peaks/',cell.type,'.bed')) %>% 
  sort

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
names(rl.chr) <- frags.chr$read %>% unique

# the worker function for co-accessibility hypothesis test
worker <- function(i) {
  peak <- resize(peaks.chr[i], width = 1e3, fix = 'center')
  control <- resize(peaks.chr[i], width = 1e5, fix = 'center')
  fr <- subset(frags.chr, start >= start(peak) & end <= end(peak))
  fr.control <-
    subset(frags.chr, start >= start(control) & end <= end(control))
  result <- list()
  for (sd in c('-', '+')) {
    index <- subset(fr, side == sd)$read
    index.control <- subset(fr.control, side == sd)$read
    if (length(index) != 0) {
      rl.chr.peak <- rl.chr[index]
      rl.chr.control <- rl.chr[index.control]
      pval <- ks.test(rl.chr.peak, rl.chr.control)$p.value
      result[[sd]] <-
        c(pval, median(rl.chr.peak), median(rl.chr.control))
    } else{
      result[[sd]] <- c(NA, NA, NA)
    }
  }
  result
}

result <- mclapply(seq_along(peaks.chr), worker, mc.cores = cores)

# foramt hypothesis test results
peak.test <- data.frame(
  chr = seqnames(peaks.chr),
  start = start(peaks.chr),
  end = end(peaks.chr),
  left.p.val = sapply(result, function(x) x[['-']][1]),
  left.length = sapply(result, function(x) x[['-']][2]),
  left.control= sapply(result, function(x) x[['-']][3]),
  right.p.val = sapply(result, function(x) x[['+']][1]),
  right.length = sapply(result, function(x) x[['+']][2]),
  right.control = sapply(result, function(x) x[['+']][3])
) %>%
  mutate(left.p.val.adj=p.adjust(left.p.val,method = 'fdr'),
          right.p.val.adj=p.adjust(right.p.val,method = 'fdr'))

saveRDS(peak.test, paste0('./ks_test/nanoATAC.',cell.type,'.chr',chrom,'.co-accessibility.ks.Rds')) 
