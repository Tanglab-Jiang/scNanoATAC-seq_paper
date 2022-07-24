library(dplyr)
library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(genomation)
library(magrittr)
library(ggsci)
library(eulerr)

STD_SV<-c('DEL','INS','DUP','INV','TRA')

setwd('~/SHARE/HYQ_atac/analysis/SV/')

#### Functions ####
if (T) {
  loadSV <-
    function(file = NULL,
             genome = 'hg38',
             chr.used = c(1:22, 'X')) {
      vcf <- readVcf(file, genome)
      vcf <-
        vcf[seqnames(vcf) %in% c(chr.used, paste0('chr', chr.used))] %>% unique
      
      # check chr prefix in seqlevels
      if (all(!grepl('chr', seqlevels(vcf)))) {
        seqlevels(vcf) <- paste0('chr', seqlevels(vcf))
      }
      
      return(vcf)
    }
  
  cuteSV_parser <-
    function(file = NULL,
             genome = 'hg38',
             chr.used = c(1:22, 'X')) {
      vcf <- loadSV(file, genome, chr.used)
      
      # exclude TRA from analysis
      vcf <- vcf[vcf@info$SVTYPE != 'BND']
      
      sv.break <-
        breakpointRanges(vcf, ignoreUnknownSymbolicAlleles = T)
      
      # eliminate trival columns
      if (F) {
        sv.break$REF <- NULL
        sv.break$ALT <- NULL
        sv.break$insSeq <- NULL
        sv.break$QUAL <- NULL
        sv.break$paramRangeID <- NULL
      }
      
      return(sv.break)
    }
  
  survivor_parser <-
    function(file = NULL,
             genome = 'hg38',
             chr.used = c(1:22, 'X')) {
      vcf <- loadSV(file, genome, chr.used)
      
      vcf <- vcf[vcf@info$SVTYPE != 'TRA']
      
      vcf@fixed$QUAL <- as.numeric(vcf@info$SUPP)
      
      # fix wrong object type caused error
      info(vcf)$CIEND <- IntegerList(info(vcf)$CIEND)
      info(vcf)$CIPOS <- IntegerList(info(vcf)$CIPOS)
      info(vcf)$ID <- seq_along(vcf)
      
      sv.break <-
        breakpointRanges(vcf, ignoreUnknownSymbolicAlleles = T,info_columns='ID')
      
      return(sv.break)
    }
  
  gnomad_parser<-function(file = NULL,
                          genome = 'hg38',
                          chr.used = c(1:22, 'X')){
    vcf <- loadSV(file, genome, chr.used)
    
    vcf <- vcf[vcf@info$SVTYPE!='CNV']
    
    vcf@info$SVLEN <- as.numeric(vcf@info$SVLEN)
    
    sv.break <- breakpointRanges(vcf,ignoreUnknownSymbolicAlleles=T)
    
    return(sv.break)
  }
  
  hgsvc2_parser <- function(file = NULL,
                            genome = 'hg38',
                            chr.used = c(1:22, 'X')){
    vcf <- loadSV(file, genome, chr.used)
    sv.break <- breakpointRanges(vcf)
    
    # remove REF and ALT sequence to decrease object size
    if (F) {
      sv.break$REF<-NA
      sv.break$ALT<-NA
    }
    
    return(sv.break)
  }
  
  #  > 50bp
  filter_size <- function(gr){
    gr %>% subset(abs(svLen)>50)
  }
  
  # query: a sv break point object in GRange structure
  filter_region <- function(query, mask = NULL) {
    hits <- findOverlaps(resize(query,1,fix = 'center'), mask)
    query[-unique(hits@from)] %>% .[hasPartner(.)] 
  }
  
  select_region <- function(query, target = NULL) {
    hits <- findOverlaps(resize(query,1,fix = 'center'), target)
    query[unique(hits@from)] %>% .[hasPartner(.)] 
  }
  
  sv_diff <- function(query, subject) {
    for (i in unique(query$svtype)) {
      if (i %in% unique(subject$svtype)) {
        index <- which(query$svtype == i)
        hits <- findBreakpointOverlaps(
          query[index] ,
          subset(subject, svtype == i),
          maxgap = 500,
          sizemargin = 0.5,
          restrictMarginToSizeMultiple = 0.5
        )
        if (length(hits@from)!=0){
          query <- query[-index[hits@from]] %>% .[hasPartner(.)]
        }
      }
    }
    query
  }
  
  # gr: A break point GRanges object
  get_somatic <- function(gr) {
    for (i in ref.beak.list){
      gr <- sv_diff(gr, i)
    }
    gr <- sv_diff(gr, hg001.long.break)
    gr <- filter_region(gr, blackList)
    gr
  }
}

#### Resources ####
if (T) {
  # black list
  blackList<- readGeneric('~/SHARE/resource/hg38/hg38-blacklist.v2.bed.gz')
  
  # germline control
  if (T) {
    ref.dir <- '~/SHARE/resource/hg38/SV/'
    
    gnomad.break <-
      readRDS(file.path(ref.dir, 'gnomAD_v2_hg38_nstd166.hg38.breakPoints.Rds'))
    
    hgsvc2.break <-
      readRDS(file.path(ref.dir, 'HGSVC2.freeze3.sv.alt.hg38.breakPoints.Rds'))
    
    hg001.long.break <- 
      survivor_parser('./long_atac/2110_HG001_subset_1r_2c.sorted.vcf.gz')
    
    ref.beak.list <- list(gnomad.break, hgsvc2.break)
  }
}

#### parse SV reference to break point objects ####
if (F) {
  ref.dir<-'~/SHARE/resource/hg38/SV/'
  
  # gnomad v2
  gnomad.break <- gnomad_parser(file.path(ref.dir,'gnomAD_v2_hg38_nstd166.vcf.gz'))
  saveRDS(gnomad.break, file.path(ref.dir,'gnomAD_v2_hg38_nstd166.hg38.breakPoints.Rds'))
  
  # HGSVC2
  # Note: ID column must be UNIQUE.
  hgsvc2.break <- hgsvc2_parser(file.path(ref.dir,'HGSVC2.freeze3.sv.alt.FIX_ID.vcf.gz'))
  saveRDS(hgsvc2.break, file.path(ref.dir,'HGSVC2.freeze3.sv.alt.hg38.breakPoints.Rds'))
}

#### parse survivor merged SV sets and get somatic mutations ####
if (F) {
  sv.break <-
    survivor_parser('./long_atac/2110_K562_subset_1r_1c.sorted.vcf.gz') %>% filter_size
  sv.break <- get_somatic(sv.break)

  benchmark.break <-
    cuteSV_parser('./benchmark/2201_K562_gDNA_HYQ.cuteSV.10r.ONT.vcf.gz') %>% filter_size
  benchmark.break <- get_somatic(benchmark.break)
}

#### benchmarking ####
if (F) {
  pdf('220612_K562_somatic_SV_benchmark(K562_bulk_ONT).pdf',width = 5,height = 5)
  
  cell.cut.off <- 20
  
  worker <- function(type){
    query <- sv.break %>% subset(svtype==type)
    benchmark <- benchmark.break %>% subset(svtype==type)

    query$benchmark <-
      countBreakpointOverlaps(
        query,
        benchmark,
        maxgap = 500
      ) > 0 
    
    # number of cells supporting SV as quality proxy
    # calculate precision and recall under different quality cut-off
    
    df <- as.data.frame(query) %>%
      mutate(cell_count = QUAL, sv_type=svtype) %>%
      dplyr::select(cell_count, sv_type, benchmark) %>%
      group_by(sv_type, cell_count) %>%
      summarise(calls = n(),
                tp = sum(benchmark)) %>%
      group_by(sv_type) %>%
      arrange(desc(cell_count)) %>%
      mutate(
        cum_tp = cumsum(tp),
        cum_n = cumsum(calls),
        Precision = cum_tp / cum_n,
        Recall = cum_tp / length(benchmark)
      ) %>% 
      subset(cell_count <= cell.cut.off)
    
    df
  }
  
  df <- lapply(c('INS','DEL','DUP','INV'),worker) %>% do.call(rbind,.)
  
  # Precison ~ Quality
  ggplot(df) +
    aes(x = cell_count,
        y = Precision,
        colour = sv_type) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(1, cell.cut.off),
                       limits = c(1, cell.cut.off)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    theme_linedraw() +
    theme(panel.grid.major = element_line(color = 'grey95'),
          panel.grid.minor = element_blank()) +
    coord_fixed(cell.cut.off-1) +
    scale_color_nejm()
  
  # Recall ~ Quality
  ggplot(df) +
    aes(x = cell_count,
        y = Recall,
        colour = sv_type) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(1,cell.cut.off),limits = c(1,cell.cut.off)) +
    scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
    theme_linedraw() +
    theme(panel.grid.major = element_line(color = 'grey95'),
          panel.grid.minor = element_blank()) +
    coord_fixed(cell.cut.off-1) +
    scale_color_nejm()
  
  # Precision ~ Recall curve
  ggplot(df) +
    aes(x = Recall,
        y = Precision,
        colour = sv_type) +
    geom_point(size = 3) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
    scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1)) +
    theme_linedraw() +
    theme(panel.grid.major = element_line(color = 'grey95'),
          panel.grid.minor = element_blank()) +
    coord_fixed() +
    scale_color_nejm()
  
  dev.off()
  
  sv.result <- df %>%
    dplyr::rename(true_positive = cum_tp, count = cum_n) %>%
    dplyr::mutate(tp = NULL, calls = NULL)

  write.csv(sv.result, '220612_K562_somatic_sv_benchmark(ONT).csv', row.names = F)
}

# 220114
#### venn plot for SV benchmark ####
if (F) {
  nc = 5
  venn.df <- sv.result %>% subset(cell_count == nc)
  bp.fold <- c(2, 2, 2, 4)
  names(bp.fold) <- c('DEL', 'DUP', 'INS', 'INV')
  
  pdf('220202_K562_all_sv_ont_benchmark_venn.pdf',width = 5,height = 5)
  for (sv in c('DEL','DUP','INS','INV')){
    recall <- venn.df$Recall[venn.df$sv_type==sv] %>% as.numeric()
    precision <- venn.df$Precision[venn.df$sv_type==sv] %>% as.numeric()
    tp <- venn.df$true_positive[venn.df$sv_type==sv] / bp.fold[sv] %>% as.numeric()

    eu <- euler(c(
      "ONT_bulk&scLongATAC" = tp,
      "scLongATAC" = round(tp/precision)-tp,
      "ONT_bulk" = round(tp/recall) - tp
    ))
    
    p<-plot(
      eu,
      fills = c('steelblue1', 'orange'),
      alpha = 0.7,
      labels = list(cex = .8),
      quantities = list(cex = c(0.8), type = c("counts","percent")),
      main=paste0(sv,'_',nc,'_cells')
    )
    
    plot(p)
  }
  dev.off()
}

#### extract deletions for validation #####
if (F) {
  vcf <- loadSV('./long_atac/2110_K562_subset_1r_1c.sorted.vcf.gz', 'hg38',  c(1:22, 'X'))
  vcf <- vcf[vcf@info$SVTYPE != 'TRA']
  
  index <- sv.break %>% .[.$QUAL>=15 & .$svtype=='DEL'] %>% .$ID %>% unique
  
  writeVcf(vcf[index], 
           '220612_K562_15c.DEL.depop.hg001_2c_controlled.noBalckList.vcf')
}

#### extract insertions for validation ####
if (F) {
  vcf <- loadSV('./long_atac/2110_K562_subset_1r_1c.sorted.vcf.gz', 'hg38',  c(1:22, 'X'))
  vcf <- vcf[vcf@info$SVTYPE != 'TRA']
  
  index <- sv.break %>% .[.$QUAL>=10 & .$svtype=='INS'] %>% .$ID %>% unique
  writeVcf(vcf[index], 
           '220618_K562_10c.INS.depop.hg001_2c_controlled.noBalckList.vcf')
}
