# HYQ long-read atac-seq project

library(ArchR)
library(genomation)
library(pheatmap)
library(ggsci)
library(eulerr)

setwd("~/SHARE/HYQ_atac/")

addArchRThreads(1)
addArchRGenome("mm10")
data.path <- './data/arrow_mm10/'
proj.path <- './data/projects/main_mm10/'

# SAVE
# saveArchRProject(ArchRProj = proj,outputDirectory = proj.path)

# LOAD
proj<-loadArchRProject(proj.path)

##### functions #####
if (T) {
  # patches for ArchR
  source('../resource/archr_patches/archr_track.R')
  source('../resource/archr_patches/archr_footprint_no_ribbon.R')
  
  my_theme <- theme_ArchR(
    baseSize = 10,
    baseLineSize = 0.6,
    baseRectSize = 1.3,
    plotMarginCm = 0.5,
    legendPosition = "right",
    legendTextSize = 10
  )
  
  archrReduce <- function(obj) {
    obj <- addIterativeLSI(
      ArchRProj = obj,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      iterations = 1,
      dimsToUse = 1:10,
      force = T
    )
    
    obj <- addUMAP(
      ArchRProj = obj,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      minDist = 0.1,
      force = T
    )
    
    obj <- addClusters(
      input = obj,
      reducedDims = "IterativeLSI",
      method = "Seurat",
      name = "Clusters",
      force = T
    )
    
    obj <- addImputeWeights(obj)
    
    obj
  }
  
  # 210910
  # find overlap function with overlap proportion limitation
  findOverlapRatio <- function(query, subject, rof = NULL) {
    hits <- findOverlaps(query, subject)
    overlaps <- pintersect(query[hits@from], subject[hits@to])
    subject.overlap <- width(overlaps) / width(subject[hits@to])
    query.overlap <- width(overlaps) / width(query[hits@from])
    hits[subject.overlap > rof & query.overlap > rof]
  }
  
  # 210903
  # fix out of boundary in archr footprint
  # gr: Grange object
  fixGrange<-function(gr,flank=NULL){
    genome.gr<-getGenomeAnnotation(proj)@listData[["chromSizes"]]
    start(genome.gr)<-start(genome.gr)+flank
    end(genome.gr)<-end(genome.gr)-flank
    gr %>% .[countOverlaps(.,genome.gr)==1]
  }
  
  # 210925
  # cell metadata annotating function
  # add labels acccording to outer barcodes
  setCellGroup<-function(obj, info=NULL, lib.name=NULL){
    if(is.null(obj$group)){
      obj$group <- obj$Sample
    } 
    cell.id <- rownames(obj)
    index <-grep(lib.name, cell.id)
    obj$group[index] <-
      cell.id[index] %>% strsplit("#") %>%
      sapply(function(x) x[2]) %>% strsplit("_") %>%
      sapply(function(x) info[x[1]] )
    obj
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
  
  # 211009
  # make grange chromosome name standadard
  cleanChr<-function(gr){
    seqlevelsStyle(gr)<-'UCSC'
    subset(gr,seqnames %in% paste0('chr',c(1:22,'X','Y')))
  }
  
  # 211027
  # a wrapper for plotEmbedding
  plotEmbeddingWrap <- function(obj,name=NULL,color=NULL){
    names(color) <- unique(obj@cellColData[[name]]) %>% sort
    plotEmbedding(
      ArchRProj = obj,
      colorBy = "cellColData",
      name = name,
      embedding = 'UMAP',
      size = 0.5,
      rastr = F,
      plotAs = "points",
      pal = color
    ) + 
      my_theme + 
      guides(color = guide_legend(override.aes = list(size = 3))) 
  }
}

##### make the main project object #####
if (F) {
  files<-c('2107_ALC1','2107_ALC2')
  
  proj <- ArchRProject(
    ArrowFiles = paste0(data.path, files,'.arrow'),
    outputDirectory = proj.path,
    copyArrows = F
  )
  
  # annotation
  if (T) {
    # 2107 ALC1
    info <- c(rep('mESC',2),rep('MEF',2),rep('GM12878',2),rep('K562',2),rep('MIX',8),rep('mESC',4))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'ALC1')
    
    # 2107 ALC2
    info <- c(rep('mESC',2),rep('MEF',2),rep('GM12878',2),rep('K562',2),rep('MIX',8),rep('mESC',4))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'ALC2')
    
    # statistics
    cut.off<-1e4
    ((proj$group %in% c('MEF', 'mESC','MIX')) & (proj$nFrags>cut.off)) %>% sum()
    table(proj$group[proj$nFrags>cut.off], proj$Sample[proj$nFrags>cut.off])
    
    # quality
    plot(proj$nFrags %>% log10, proj$TSSEnrichment %>% log2, cex=.2)
    abline(v=4)
  }
  
  # QC and cell type subset
  subset.cells <- proj$cellNames %>% 
    subset(proj$group %in% c('MEF', 'mESC', 'MIX') & proj$nFrags > 1e4)
  
  proj <- subsetArchRProject(
    proj,
    cells = subset.cells,
    outputDirectory = proj.path,
    force = T
  )
  
  # basic processing of dimension reduction
  proj <- archrReduce(proj)
  
  # cluster cells into 5 cell lines
  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    maxClusters = 2,
    force = T
  )
  
  # cell line annotation
  cluster.cell <- c('mESC', 'MEF')
  names(cluster.cell) <- c('C1','C2')
  proj$cell.line <- 'Unidentified'
  for (i in seq_along(cluster.cell)) {
    cell.name<-cluster.cell[i]
    cluster.name<-names(cluster.cell)[i]
    cell.index<-(proj$group==cell.name | proj$group=='MIX')
    cluster.index<-(proj$Clusters==cluster.name)
    proj$cell.line[ cell.index & cluster.index ] <- cluster.cell[i]
  }
  
  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "cell.line",
    force = T
  )
  
  proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "cell.line",
    cutOff = 0.01,
    force = T,
    additionalParams = '--nomodel --llocal 500000'
  )
  
  if (T) {  
    cell.type <- c('MEF','mESC')
    mouse.peak.list <- lapply(cell.type, function(cl) {
      readRDS(file.path(
        './data/projects/main_mm10/PeakCalls/',
        paste0(cl, '-reproduciblePeaks.gr.rds')
      ))
    })
    names(mouse.peak.list) <- cell.type
    
    for (i in names(mouse.peak.list)){
      gr<-mouse.peak.list[[i]]
      writeGrangeAsBed(
        gr=gr,
        name = paste0('./data/projects/main_mm10/PeakCalls/',i,'.bed'))
    }
  }
  
  proj <- addPeakMatrix(proj)
  
  # add motif
  proj <- addBgdPeaks(proj)
  proj <- addMotifAnnotations(ArchRProj = proj)
  
  # Get peak to gene links
  if (F) {
    proj <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "IterativeLSI"
    )
    
    p2g <- getPeak2GeneLinks(
      ArchRProj = proj,
      corCutOff = 0.45,
      resolution = 1,
      returnLoops = FALSE
    )
  }
  
  # marker genes & peaks& motifs
  if (F) {
    markersGS <- getMarkerFeatures(
      ArchRProj = proj %>% subset(proj$cell.line!='Unidentified'),
      useMatrix = "GeneScoreMatrix",
      groupBy = "cell.line",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      maxCells = 1000,
      threads = 6
    )
    
    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj %>% subset(proj$cell.line!='Unidentified'),
      useMatrix = "PeakMatrix",
      groupBy = "cell.line",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      maxCells = 1000,
      threads = 6
    )
    
    markerMotifs <- peakAnnoEnrichment(
      seMarker = markersPeaks,
      ArchRProj = proj %>% subset(proj$cell.line!='Unidentified'),
      peakAnnotation = "Motif",
      cutOff = "FDR <= 0.01 & Log2FC >= 1"
    )
  }
}

##### umap for mouse single cells  #####
if (F) {
  pdf('211030_mouse_umap.pdf', width = 6, height = 6)
  plotEmbeddingWrap(proj, name = 'Sample', color = pal_nejm()(9))
  plotEmbeddingWrap(proj, name = 'group', color = pal_nejm()(9))
  plotEmbeddingWrap(proj, name = 'cell.line', color = pal_nejm()(9))
  dev.off()
}

# 220106
##### cell id susbet for each cell line #####
if (F) {
  for (i in unique(proj$cell.line)) {
    proj$cellNames[proj$cell.line == i] %>% 
      sort %>% 
      write(paste0('./data/subset/', i, '_subset.txt'))
  }
}

# 220108
#### output coverage tracks ####
if (F) {
  getGroupBW(
    proj,
    groupBy = "cell.line",
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e3
  )
}
