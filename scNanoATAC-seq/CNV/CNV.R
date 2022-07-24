library(Alleloscope) 
library(ArchR)
library(Biobase)
library(forecast)  # ma 
library(genomation)
library(pheatmap)
library(parallel)

setwd("./analysis/CNV/")

# Resources 
if (F) {
  addArchRThreads(1)
  addArchRGenome("hg38")
  
  proj.100kb <- loadArchRProject('../../data/projects/main_100kb_noBinarize/')
  gap <- readBed('../../resource/hg38/hg38_gap.bed')
  
  tile.matrix <- readRDS('./data/220120_cell_line_tile_matrix_100kb.Rds')
}

# Functions
if (T) {
  # 220102
  # get coverage from tile matrix with any bin bigger than original size
  getTileMatrix <- function(obj, chr, archrBin = 1e5) {
    tile.obj <-
      getMatrixFromProject(
        obj,
        useMatrix = 'TileMatrix',
        binarize = F,
        useSeqnames = chr
      )
    
    tile.matrix <- tile.obj@assays@data@listData[["TileMatrix"]]
    tile.gr <- GRanges(
      seqnames = tile.obj@elementMetadata$seqnames,
      ranges = IRanges(start = tile.obj@elementMetadata$start,
                       width = archrBin)
    )
    
    rownames(tile.matrix) <-
      paste0(seqnames(tile.gr), '-', ranges(tile.gr))
    
    return(tile.matrix)
  }
  
  # filter black list and gap regions
  filter_tile_matrix <- function(tile.matrix, filter.region) {
    coo <- rownames(tile.matrix)
    gr <-
      GRanges(
        seqnames = sapply(strsplit(coo, '-'), function(x)
          x[1]),
        ranges = IRanges(
          start = sapply(strsplit(coo, '-'), function(x)
            as.numeric(x[2])),
          end =  sapply(strsplit(coo, '-'), function(x)
            as.numeric(x[3]))
        )
      )
    
    hits <- findOverlaps(gr, filter.region)
    
    tile.matrix <- tile.matrix[-unique(hits@from),]
    
    tile.matrix
  }
  
  # information of single-cell CNA
  get_cell_info <- function() {
    cell.type <<- cbind(proj.100kb$cellNames,proj.100kb$cell.line)
    gr <- proj.100kb@genomeAnnotation$chromSizes
    chr.size <<- cbind(gr@seqinfo@seqnames, gr@seqinfo@seqlengths)
  }
  
  get_pseudo_bulk <- function(){
    get_cell_info()
    tile.matrix.merged <- lapply(unique(cell.type[, 2]), function(x) {
      rowMeans(tile.matrix[, cell.type[, 1][cell.type[, 2] == x]])
    }) %>% do.call(cbind, .)
    colnames(tile.matrix.merged) <- unique(cell.type[, 2])
    tile.matrix.merged <<- tile.matrix.merged[, colnames(tile.matrix.merged) != 'Unidentified']
  }
}

# Divide the genome into 100kb bins
if (F) {
  proj.100kb <-
    addTileMatrix(
      proj,
      tileSize = 1e5,
      binarize = F,
      excludeChr = NULL,
      force = T,
      threads = 6
    )
  
  saveArchRProject(ArchRProj = proj.100kb,
                   outputDirectory = '../../data/projects/main_100kb_noBinarize/')
}

#### long atac-seq tile matrix for CNA plot by Alleloscope ####
if (F) {
  tile.list <-
    lapply(paste0('chr', c(1:22, 'X')), function(chr)
      getTileMatrix(obj = proj.100kb, chr = chr))
  tile.matrix <- do.call(rbind, tile.list)
  
  saveRDS(tile.matrix,'./data/220120_cell_line_tile_matrix_100kb.Rds')
}

# 220105
##### plot smoothed CNA scatter by chromosome #####
if (F) {
  # normalize
  tile.matrix.merged <- t(t(tile.matrix.merged) / colMeans(tile.matrix.merged)) * 2
  
  # control
  tile.matrix.merged<-tile.matrix.merged / rowMedians(tile.matrix.merged) * 2
  
  # plot
  par(mfrow = c(ncol(tile.matrix.merged), 1),mar=rep(2,4))
  for (i in 1:ncol(tile.matrix.merged)) {
    copy.fc <- tile.matrix.merged[,i]
    
    names(copy.fc) <-
      rownames(tile.matrix.merged) %>% 
      strsplit('-') %>% 
      sapply(function(x)x[1])
    copy.fc <- copy.fc[is.finite(copy.fc)]
    
    chr.len <- c()
    copy.fc.chr <- c()
    for (chr in paste0('chr', c(1:22, 'X'))) {
      chr.len <- c(chr.len, sum(names(copy.fc) == chr))
      if (doSmooth) {
        copy.fc.chr <-
          c(copy.fc.chr, copy.fc %>% .[names(.) == chr] %>% ma(10) )
      } else{
        copy.fc.chr <-
          c(copy.fc.chr, copy.fc %>% .[names(.) == chr])
      }
    }
    
    plot(
      x = seq_along(copy.fc.chr),
      y = copy.fc.chr,
      cex = .01,
      main = colnames(tile.matrix.merged)[i],
      xlim = c(0, length(copy.fc.chr)),
      ylim = c(0, 8),
      xaxt = 'n'
    )
    
    abline(v = c(0, cumsum(chr.len)))
    
    text(x = cumsum(chr.len) - 0.5 * chr.len %>% round,
         y = .2,
         c(1:22, 'X'))
  }
}

# 220106
# published
##### single-cell matrix normalized and controlled ####
if (F) {
  get_pseudo_bulk()
  
  cell.line<-unique(proj.100kb$cell.line) %>% .[.!='Unidentified']
  
  for (cl in cell.line){
    tile.matrix.sc <-
      tile.matrix[,proj.100kb$cellNames[proj.100kb$cell.line==cl]]
    
    # normalize
    tile.matrix.sc <-
      t(t(tile.matrix.sc) / colMeans(tile.matrix.sc)) * 2
    
    # control by median of all pseudo bulk
    control <- rowMedians(tile.matrix.merged)
    control <- control / mean(control) * 2
    tile.matrix.sc <- tile.matrix.sc[control>0, ] / control[control>0] * 2
    
    # smooth
    ma_win <- 100
    chr.name <- rownames(tile.matrix.sc) %>%
      strsplit('-') %>%
      sapply(function(x)
        x[1])
    smooth.matrix <-
      mclapply(mc.cores = 6,paste0('chr', c(1:22, 'X')), function(chr) {
        chr.matrix <- tile.matrix.sc[chr.name == chr,]
        lapply(1:ncol(chr.matrix), function(i) {
          ma(chr.matrix[, i], ma_win)
        })  %>% do.call(cbind, .)
      })  %>% do.call(rbind, .)
    rownames(smooth.matrix) <- rownames(tile.matrix.sc)
    colnames(smooth.matrix) <- colnames(tile.matrix.sc)
    tile.matrix.sc <- smooth.matrix
    
    # denoise and clip matrix
    tile.matrix.sc[tile.matrix.sc <= 2.5 & tile.matrix.sc >= 1.5] <- 2
    tile.matrix.sc[tile.matrix.sc > 4] <- 4
    tile.matrix.sc[is.na(tile.matrix.sc)] <- 2
    
    # filter region
    filter.region <- c(proj.100kb@genomeAnnotation$blacklist, gap)
    tile.matrix.sc <-
      filter_tile_matrix(tile.matrix.sc, filter.region)
    
    saveRDS(tile.matrix.sc,paste0('./data/',cl,'_CNA_matrix.Rds'))
  }
}

# plot heatmap
if (F) {
  cell.line <- unique(proj.100kb$cell.line) %>% .[.!='Unidentified']
  for (cl in cell.line) {
    tile.matrix.sc <- readRDS(paste0('./data/', cl, '_CNA_matrix.Rds'))
    if (cl=='eHAP1'){
      tile.matrix.sc<-tile.matrix.sc/2
    }
    
    cell.count <- ncol(tile.matrix.sc)
    
    png(
      filename = paste0(cl, '_long_ATAC_CNA.png'),
      width = 1200,
      height = 1200
    )
    
    chr.name <- rownames(tile.matrix.sc) %>%
      strsplit('-') %>%
      sapply(function(x)
        x[1]) %>%
      factor(levels = paste0('chr', c(1:22, 'X')))
    
    labels_col <- rep(' ', length(chr.name))
    gaps_col <- cumsum(table(chr.name))
    labels_col[gaps_col - round(0.5 * table(chr.name))] <-
      as.character(unique(chr.name))
    
    breaklength <- 50
    pheatmap(
      t(tile.matrix.sc),
      cluster_cols = F,
      cluster_rows = T,
      show_colnames = T,
      show_rownames = F,
      cellheight = 0.25,
      fontsize = 20,
      treeheight_row = 0,
      color = colorRampPalette(c("blue", "white", "red"))(breaklength),
      breaks = seq(0, 4, length.out = breaklength),
      clustering_distance_rows = "correlation",
      clustering_method = "ward.D2",
      gaps_col = gaps_col,
      labels_col = labels_col
    )
    
    dev.off()
  }
}

