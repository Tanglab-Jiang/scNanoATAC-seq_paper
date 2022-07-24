# HYQ long-read atac-seq project

library(ArchR)
library(genomation)
library(pheatmap)
library(mgcv)
library(ggsci)
library(lsa) # cosine similarity
library(eulerr)
library(gtools)
library(ggpubr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(irlba)
library(forecast)
library(rtracklayer)

setwd("~/SHARE/HYQ_atac/")

addArchRThreads(1)
addArchRGenome("hg38")
hg38 <- BSgenome.Hsapiens.UCSC.hg38
data.path <- './data/arrow/'
proj.path <- './data/projects/main/'

# SAVE
# saveArchRProject(ArchRProj = proj,outputDirectory = proj.path)

##### resources #####
if (T) {
  # load the main project
  proj <- loadArchRProject(proj.path)
  
  if (F) {
    df <- fread('./resource/ENCODE_cCREs_full.bed',header = T,sep = '\t') %>% as.data.frame()
    ccres <- GRanges(
      seqnames = df$`#chrom`,
      ranges = IRanges(df$chromStart, df$chromEnd),
      type = df$encodeLabel
    )
    # saveRDS(ccres,'./resource/210926_ENCODE_cCREs_GRanges.Rds')
  }
  
  ccres <- readRDS('./resource/210926_ENCODE_cCREs_GRanges.Rds')
  
  cell.type.order <-
    c("eHAP1", "GM12878", "HEK293T", "HFF1", "K562", "Unidentified")
  
  proj.granja <- loadArchRProject('./data/granja_2021/CellLine/')
  
  ##### load peak called by ArchR by cell line #####
  if (T) {
    # Load TGS atac-seq peak set by cell type 
    cell.type <- c('K562','GM12878','HEK293T','eHAP1','HFF1')
    tgs.peak.list <- lapply(cell.type, function(cl) {
      readRDS(file.path(
        './data/projects/main/PeakCalls/',
        paste0(cl, '-reproduciblePeaks.gr.rds')
      ))
    })
    names(tgs.peak.list) <- cell.type
    
    # Load 10x atac-seq peak set by cell type 
    cell.type <- c('K562','GM12878','HEK293T')
    tenx.peak.list <- lapply(cell.type, function(cl) {
      readRDS(file.path(
        './data/granja_2021/CellLine/PeakCalls/',
        paste0(cl, '-reproduciblePeaks.gr.rds')
      ))
    })
    names(tenx.peak.list) <- cell.type
  }
}

##### functions #####
if (T) {
  # patches for ArchR
  source('../resource/archr_patches/archr_track.R')
  
  # change theme for ArchR plot
  if (T) {
    my_theme <- theme_ArchR
    formals(my_theme)$baseLineSize <- 0.6
    formals(my_theme)$baseRectSize <- 1.3
    formals(my_theme)$plotMarginCm <- 0.5
    formals(my_theme)$legendPosition <- "right"
    formals(my_theme)$legendTextSize <- 10
    
    assignInNamespace('theme_ArchR',my_theme,'ArchR')
  }
  
  archrReduce <- function(obj) {
    obj <- addIterativeLSI(
      ArchRProj = obj,
      useMatrix = "TileMatrix",
      name = "IterativeLSI",
      iterations = 2,
      dimsToUse = 1:30,
      force = T
    )
    
    obj <- addUMAP(
      ArchRProj = obj,
      reducedDims = "IterativeLSI",
      name = "UMAP",
      minDist = 1,
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
  fitGam <- function(x, y) {
    df <- data.frame(x = x, y = y)
    fit <- gam(y ~ s(x),
               data = df)
    plot.gam(fit, shift = fit$coefficients['(Intercept)'])
  }
  
  # 210910
  # find overlap function with overlap proportion limitation
  # deprecated for shrinking consistency among datasets (211202)
  findOverlapRatio <- function(query, subject, rof = 0) {
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
  
  # 210916
  # Caculate enrichment fold change against neighboring basline.
  # Ensure mutual exclusion between flank and center.
  # coo: Grange object
  # window: suppose to be equal to the bin size of tile matrix 
  getEnrichment <-
    function(obj,
             coo = NULL,
             chr = NULL,
             window = 500,
             flank = 2000) {
      
      if (!is.null(chr)){
        coo <- coo[seqnames(coo)==chr]
      }
      
      coo <- resize(coo, 1, fix = "center")
      coo <- unique(coo)
      center.coo <- resize(coo, window, "center") %>% fixGrange(.,flank = window/2)
      
      flank.coo <- c(
        GRanges(seqnames(coo), IRanges(end(coo) + flank - window + 1, end(coo) + flank)),
        GRanges(seqnames(coo), IRanges(start(coo) - flank, start(coo) - flank + window - 1))
      )
      
      exclude.index<-findOverlaps(center.coo,flank.coo)@to
      if (length(exclude.index)!=0){
        flank.coo <- flank.coo %>% .[-exclude.index] %>% fixGrange(.,flank = flank)
      }
      
      # get tile matrix
      tile.obj<-
        getMatrixFromProject(
          obj,
          useMatrix = 'TileMatrix',
          binarize = T,
          useSeqnames = chr)
      
      tile.gr <- GRanges(
        seqnames = tile.obj@elementMetadata$seqnames,
        ranges = IRanges(start = tile.obj@elementMetadata$start,
                         width = window)
      )
      
      # get tile index
      center.index <- findOverlaps(center.coo, tile.gr)@to
      flank.index <- findOverlaps(flank.coo, tile.gr)@to
      
      # count fragments in region of interest
      center.count <- tile.obj@assays@data@listData[["TileMatrix"]][,obj$cellNames][center.index,] %>% colSums
      flank.count <- tile.obj@assays@data@listData[["TileMatrix"]][,obj$cellNames][flank.index,] %>% colSums
      names(center.count) <- obj$cellNames
      names(flank.count) <- obj$cellNames
      
      list(cc=center.count,fc=flank.count,
           lc=length(center.index),lf=length(flank.index))
  }
  
  # 210908
  # get enrichment score of whole genome directly
  # coo: Grange object
  getEnrichmentScore<-function(obj,coo){
    worker <- function(x) {
      getEnrichment(
        obj,
        coo = coo,
        chr = x
      )
    }
    chrs<-coo %>% seqnames %>% unique %>% as.character
    
    result<-lapply(chrs, worker)
    names(result)<-chrs
    
    center<-lapply(chrs, function(x){result[[x]]$cc})
    count.center<-as.data.frame(center) %>% rowSums
    
    flank<-lapply(chrs, function(x){result[[x]]$fc})
    count.flank<-as.data.frame(flank) %>% rowSums
    
    length.center<-sapply(chrs, function(x){result[[x]]$lc}) %>% sum
    length.flank<-sapply(chrs, function(x){result[[x]]$lf}) %>% sum
    
    enrichment.score<-(count.center/length.center)/(count.flank/length.flank)
    enrichment.score[count.flank==0]<-NA
    
    enrichment.score
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
  # run a benchmark test
  ccres_benchmark<-function(peak.set,group=NULL){
    result <- list()
    
    for (i in names(peak.set)) {
      query <- peak.set[[i]]
      seqlevelsStyle(query) <- 'UCSC'
      hits <-
        findOverlaps(query = query, subject = ccres)
      
      precision <- length(unique(hits@from)) / length(query)
      recall <- length(unique(hits@to)) / length(ccres)
      F1 <- 2 / (1 / precision + 1 / recall)
      
      result[[i]]<-do.call(rbind,
                           list(data.frame(name=i,value=precision,metric='Precision'),
                                data.frame(name=i,value=recall,metric='Recall'),
                                data.frame(name=i,value=F1,metric='F1')))
    }
    
    ccres.bench.df<-do.call(rbind, result)
    ccres.bench.df$metric<-factor(ccres.bench.df$metric, c('Precision','Recall','F1'))
    ccres.bench.df$group<-group
    
    ccres.bench.df
  }
  
  # 211201
  # a wrapper for plotEmbedding
  plotEmbeddingWrap <- function(obj,name=NULL,color=NULL){
    names(color) <- obj@cellColData[[name]] %>% unique %>% 
      as.character %>% gtools::mixedsort(.)
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
      my_theme() + 
      guides(color = guide_legend(override.aes = list(size = 3))) 
  }
  
  # 211212
  # get coverage from tile matrix with any bin bigger than original size
  getCoverage <- function(obj,chr,archrBin=500,bin=1e5){
    tile.obj <-
      getMatrixFromProject(obj,
                           useMatrix = 'TileMatrix',
                           binarize = T,
                           useSeqnames = chr)
    
    tile.mean <- tile.obj@assays@data@listData[["TileMatrix"]] %>% rowMeans()
    tile.gr <- GRanges(
      seqnames = tile.obj@elementMetadata$seqnames,
      ranges = IRanges(start = tile.obj@elementMetadata$start,
                       width = archrBin)
    )
    tile.gr$coverage<-tile.mean
    
    gr.windows <-
      tileGenome(seqinfo(hg38),
                 tilewidth = bin) %>% unlist
    gr.data.cov <- coverage(tile.gr, weight="coverage")
    seqlevels(gr.windows, pruning.mode="coarse") <- names(gr.data.cov)
    gr <- binnedAverage(gr.windows, gr.data.cov, "coverage")
    
    return(gr)
  }
}

# make the main project object
if (F) {
  files <- c('2105_AL480','2107_AL960_1','2107_AL960_2',
           '2107_ALC1','2107_ALC2','2108_EHF',
           '2110_GK','2110_HH')
  
  proj <- ArchRProject(
    ArrowFiles = paste0(data.path, files,'.arrow'),
    outputDirectory = proj.path,
    copyArrows = F
  )
  
  # annotation
  if (T) {
    # 2105 AL480
    info <- rep(c('p-HG001', 'K562', 'HEK293T', 'HFF1', 'p-HepG2'), each = 2)
    names(info) <- 81:90
    proj <- setCellGroup(proj, info, 'AL480')
    
    # 2107 AL960_1
    info <- rep(rep(c('K562', 'HEK293T', 'HFF1', 'GM12878', 'eHAP1'), each = 2), 2)
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'AL960_1')
    
    # 2107 AL960_2
    info <- rep(rep(c('K562', 'HEK293T', 'HFF1', 'GM12878', 'eHAP1'), each = 2), 2)
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'AL960_2')
    
    # 2107 ALC1
    info <- c(rep('mESC',2),rep('MEF',2),rep('GM12878',2),rep('K562',2),rep('MIX',8),rep('mESC',4))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'ALC1')
    
    # 2107 ALC2
    info <- c(rep('mESC',2),rep('MEF',2),rep('GM12878',2),rep('K562',2),rep('MIX',8),rep('mESC',4))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, 'ALC2')
    
    # 2108 EHF
    info <- c(rep('eHAP1',2),rep('eHAP1_2n',2),rep('ZSF-NT-EPI',6),rep('HFF1',4),rep('ZSF-NT-FIB',6))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, '2108_EHF')
    
    # 2110 GK
    info <- c(rep('GM12878',10),rep('K562',10))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, '2110_GK')
    
    # 2110 HH
    info <- c(rep('HEK293T',12),rep('HFF1',8))
    names(info) <- 73:92
    proj <- setCellGroup(proj, info, '2110_HH')
    
    # statistics
    cut.off<-1e4
    ((proj$group %in% c('GM12878', 'K562', 'HFF1', 'eHAP1', 'HEK293T','MIX')) & (proj$nFrags>cut.off)) %>% sum()
    table(proj$group[proj$nFrags>cut.off], proj$Sample[proj$nFrags>cut.off])
    
    # quality
    plot(proj$nFrags %>% log10, proj$TSSEnrichment %>% log2, cex=.2)
    abline(v=4)
  }
  
  # QC and cell type subset
  subset.cells <- proj$cellNames %>% 
    subset(proj$group %in% c('GM12878', 'K562', 'HFF1', 'eHAP1', 'HEK293T', 'MIX') & proj$nFrags > 1e4)
  
  proj <- subsetArchRProject(
    proj,
    cells = subset.cells,
    outputDirectory = proj.path,
    force = T)
  
  # basic processing of dimension reduction
  proj <- archrReduce(proj)
  
  # cluster cells into 5 cell lines
  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    maxClusters = 5,
    force = T
  )

  # cell line annotation
  cluster.cell <- c('K562', 'GM12878', 'HFF1', 'HEK293T', 'eHAP1')
  names(cluster.cell) <- paste0('C', 1:5)
  proj$cell.line <- 'Unidentified'
  for (i in seq_along(cluster.cell)) {
    cell.name<-cluster.cell[i]
    cluster.name<-names(cluster.cell)[i]
    cell.index<-(proj$group==cell.name | (proj$group=='MIX' & cell.name %in% c('GM12878','K562')))
    cluster.index<-(proj$Clusters==cluster.name)
    proj$cell.line[ cell.index& cluster.index ] <- cluster.cell[i]
  }
  
  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "cell.line"
  )
  
  if (T) {
    proj <- addReproduciblePeakSet(
      ArchRProj = proj,
      groupBy = "cell.line",
      cutOff = 0.01,
      additionalParams = '--nomodel --llocal 500000',
      threads = 3
    ) 
    
    cell.type <- unique(proj$cell.line)[-6]
    peak.list <- lapply(cell.type, function(cl) {
      readRDS(file.path(
        proj.path,
        'PeakCalls/',
        paste0(cl, '-reproduciblePeaks.gr.rds')
      ))
    })
    names(peak.list) <- cell.type
    
    for (i in names(peak.list)) {
      writeGrangeAsBed(peak.list[[i]],
                       paste0(proj.path, '/PeakCalls/', i,'.bed'))
    }
  }
  
  proj <- addPeakMatrix(proj)
  
  # add motif
  proj <- addBgdPeaks(proj)
  proj <-
    addMotifAnnotations(ArchRProj = proj,
                        motifSet = "cisbp",
                        force = T)
  
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
    
    saveRDS(
      list(
        markersGS = markersGS,
        markersPeaks = markersPeaks,
        markerMotifs = markerMotifs
      ),
      './data/211021_marker.Rds'
    )
  }
}

##### Dimension reduction & QC metrics #####
if (F) {
  pdf('220625_human_cell_line.pdf',
      width = 5,
      height = 5)
  
  plot(
    proj$cell.line %>% factor,
    proj$TSSEnrichment,
    xlab = 'Cell Line',
    ylab = 'TSS Enrichment',
    col = c(pal_lancet()(5),'grey'),
    cex = .25,
    cex.axis = .7,
    pch = 20
  )
  
  plot(
    proj$cell.line %>% factor,
    proj$nFrags %>% log10,
    xlab = 'Cell Line',
    ylab = 'log10(Fragments)',
    col = c(pal_lancet()(5),'grey'),
    cex = .25,
    cex.axis = .7,
    pch = 20
  )
  
  proj$FRIP %>%
  plot(
    proj$cell.line %>% factor,
    .,
    xlab = 'Cell Line',
    ylab = 'FRIP',
    col = c(pal_lancet()(5),'grey'),
    cex = .25,
    cex.axis = .7,
    pch = 20,
    ylim=c(0,max(.))
  )
  
  plotEmbeddingWrap(proj,'cell.line',c(pal_lancet()(5), 'grey'))
  plotEmbeddingWrap(proj,'Sample',c(pal_lancet()(9)))
  
  dev.off()
}

# 210927
##### save bigwig #####
if (F) {
  getGroupBW(
    ArchRProj = proj,
    groupBy = "cell.line",
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e4
  )
}

# 211007
##### cell id susbet for each cell line #####
if (F) {
  for (i in unique(proj$cell.line)) {
    proj$cellNames[proj$cell.line == i] %>% 
      sort %>% 
      write(paste0('./data/subset/', i, '_subset.txt'))
  }
}


#### TSS enrichment of cell lines ####
if (F) {
  pdf('220625_TSS_enrichment_with_ribbon.pdf', width = 5, height = 5)
  
  tss.foot <- getFootprints(
    ArchRProj = proj,
    positions = list(TSS = getTSS(proj)),
    groupBy = "cell.line",
    flank = 2000
  )
  
  tss.plot <- plotFootprints(
    seFoot = tss.foot,
    ArchRProj = proj,
    normMethod = "Subtract",
    smoothWindow = 50,
    plot = F,
    pal = pal_lancet()(9)
  )
  
  grid.newpage()
  grid.draw(tss.plot[1][[1]])
  
  dev.off()
}

##### nanoATAC motif footprint #####
if (F) {
  pdf('220625_CTCF_footprint_with_ribbon.pdf',width = 5,height = 7)
  
  motif.pos <- getPositions(proj,name='Motif')
  
  seFoot <- getFootprints(
    ArchRProj = proj,
    positions = list(ctcf = motif.pos$CTCF_177),
    groupBy = "cell.line"
  )
  
  p<-plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "subtract",
    plot = F,
    smoothWindow = 10,
    pal=pal_lancet()(9)
  )
  
  for (i in 1:length(p)){
    grid::grid.newpage()
    grid::grid.draw(p[[i]],recording = F)
  }
  
  dev.off()
}

##### plot motif foot prints #####
if (F) {
  pdf('220624_cellLine_TF_footprint.pdf',width = 5.5,height = 7)
  
  source('../resource/archr_patches/archr_footprint_no_ribbon.R')
  
  motif.pos <- getPositions(proj,name='Motif')
  
  seFoot <- getFootprints(
    ArchRProj = proj,
    positions = list(
      CTCF = motif.pos$CTCF_177,
      CTCFL = motif.pos$CTCFL_198,
      BCL11A = motif.pos$BCL11A_194,
      STAT2 = motif.pos$STAT2_778,
      JUNB = motif.pos$JUNB_139,
      GATA1 = motif.pos$GATA1_383,
      FOS = motif.pos$FOS_137,
      FOSL1 = motif.pos$FOSL1_142,
      JUND = motif.pos$JUND_124,
      BACH1 = motif.pos$BACH1_130,
      FOSB = motif.pos$FOSB_121,
      FOXP4 = motif.pos$FOXP4_358,
      FOXO3 = motif.pos$FOXO3_354
    ),
    groupBy = "cell.line"
  )
  
  p <- plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj,
    normMethod = "subtract",
    plot = F,
    addDOC = T,
    baseSize = 10,
    smoothWindow = 10,
    pal = pal_lancet()(9)
  )
  
  for (i in 1:length(p)){
    grid::grid.newpage()
    grid::grid.draw(p[[i]],recording = F)
  }
  
  dev.off()
}


# 211028
# Granja, et al., 2021
#### Downsample GM12878, 293T, K562 cells to equivalent numbers with TGS for cCREs benchmark test ####
if (F) {
  proj.granja$cell.line[proj.granja$cell.line=='293T']<-'HEK293T'
  cells<-c()
  for (cl in c('K562','GM12878','HEK293T')){
    cells<-c(cells,
             sample(proj.granja$cellNames[proj.granja$cell.line==cl],
                    sum(proj$cell.line==cl)))
  }
  
  proj.granja.down <- proj.granja[cells]
  
  proj.granja.down <- addGroupCoverages(
    ArchRProj = proj.granja.down,
    groupBy = "cell.line",
    force = T
  )
  
  proj.granja.down <- addReproduciblePeakSet(
    ArchRProj = proj.granja.down, 
    groupBy = "cell.line",
    cutOff = 0.01,
    force = T
  )
  
  #saveRDS(proj.granja.down@peakSet,'./data/211028_granja_downsampled_peakSet.Rds')
  #saveArchRProject(ArchRProj = proj.granja.down, outputDirectory = './data/granja_2021/CellLine_downsampled/')
}


##### cell line peak calling comparasion between 10x and TGS #####
if (F) {
  pdf('220625_cell_line_peak_10x_nanoATAC.pdf',width = 5,height = 5)
  
  for (cl in cell.type){
    query <- tgs.peak.list[[cl]]
    subject <- tenx.peak.list[[cl]]
    hits <- findOverlaps(query = query,subject = subject)
    
    isec.q<-length(unique(hits@from))
    isec.s<-length(unique(hits@to))
    
    eu <- euler(c(
      "ONT&10x" = round((isec.q+isec.s)/2),
      "ONT" = length(query) - isec.q,
      "10x" = length(subject) - isec.s
    ))
    
    p<-plot(
      eu,
      fills = c('steelblue1', 'orange'),
      alpha = 0.7,
      labels = list(cex = .8),
      quantities = list(cex = c(0.8), type = c("counts","percent")),
      main=paste0(cl,' (Granja et al., 2021)')
    )
    
    plot(p)
  }
  
  dev.off()
}

##### take ENCODE cCREs as a benchmark #####
if (F) {
  ccres.bench.df<-rbind(
    ccres_benchmark(tenx.peak.list[c('K562','GM12878','HEK293T')], group = '10X'),
    ccres_benchmark(tgs.peak.list[c('K562','GM12878','HEK293T')], group = 'TGS')
  )
  
  ccres.bench.df$group <- factor(ccres.bench.df$group, c('10X','TGS'))
  
  pdf('220625_llocal_500k_cell_line_cCREs_benchmark.pdf',width = 12,height = 4)
  
  dummy <- ccres.bench.df
  dummy$value[dummy$metric == 'Precision'] <- 1
  dummy$value[dummy$metric == 'Recall'] <- .2
  dummy$value[dummy$metric == 'F1'] <- .2
  
  ggplot(ccres.bench.df) +
    geom_bar(aes(x = name, y = value, fill = group),position = 'dodge', stat = "identity") + 
    geom_blank(data=dummy,aes(x = name, y = value, fill = group)) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values=c('steelblue','orange')) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, color = 'black'),
          axis.text.y = element_text(color = 'black'),
          strip.text = element_text(colour = 'black',size=15))
  
  dev.off()
}

##### cCREs enrichment #####
if (F) {
  peak.list<-tgs.peak.list
  #peak.list<-tenx.peak.list
  
  # intersect peaks with cCREs by category
  result<-list()
  for (i in names(peak.list)){
    for (j in unique(ccres$type)){
      gr <- ccres %>% subset(type==j)
      hits <- findOverlaps(peak.list[[i]], gr)
      count <- length(unique(hits@from))
      
      result <-
        c(result, list(data.frame(
          cell.type = i,
          element.type = j,
          count = count
        )))
    }
    
    # count non cCRE peaks
    hits <- findOverlaps(peak.list[[i]], ccres)
    count <- length(peak.list[[i]])-length(unique(hits@from))
    result<-c(result,list(data.frame(cell.type=i,element.type='None',count=count)))
  }
  ce<-do.call(rbind,result)
  
  # default color palette for cCREs
  ccres.pal <-
    c('#FF0000', '#FFA700', '#FFCD00', '#00B0F0', '#FFAAAA', 'grey70')
  names(ccres.pal) <-
    c('PLS', 'pELS', 'dELS', 'CTCF-only', 'DNase-H3K4me3', 'None')
  ce$element.type <- factor(
    ce$element.type,
    levels = c('PLS', 'pELS', 'dELS', 'CTCF-only', 'DNase-H3K4me3', 'None')
  )
  
  # plot cCREs enrichment
  for (plot_none in c(T,F)){
    pdf(paste0('220625_llocal_500kb_TGS_cCREs_enrichment_(none_',plot_none,').pdf'), width = 6, height = 6)
    
    for (ct in unique(ce$cell.type)){
      df<-ce %>% subset(cell.type == ct)
      df$ratio <- df$count/sum(df$count)
      df$label<-round(df$ratio*100,digits = 1) %>% paste0(.,'%')
      
      if (!plot_none) df<-subset(df,element.type!='None')
      
      p <- ggdonutchart(
        data = df,
        x = 'ratio',
        label = 'label',
        lab.pos = 'out',
        fill = "element.type",
        color = "black",
        lab.font = c(6, 'bold', 'black'),
        palette = ccres.pal,
        ggtheme = theme_pubr(),
        size = .75
      ) +
        annotate(
          'text',
          x = 1,
          y = 0,
          label = ct,
          vjust = 3,
          size = 7
        )
      
      plot(p)
    }
    
    dev.off()
  }
}

#### Peak Annotation by Relationship with Genes ####
if (F) {
  peak.list<-tgs.peak.list
  #peak.list<-tenx.peak.list
  
  pdf('220625_llocal_500kb_TGS_peakAnno.pdf', width = 8, height = 6)
  
  for (ct in names(peak.list)){
    peakAnno <-
      annotatePeak(
        peak.list[[ct]],
        tssRegion = c(-1e3, 1e3),
        TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        annoDb = "org.Hs.eg.db"
      )
    
    df<-peakAnno@annoStat
    df$label<-round(df$Frequency,digits = 1) %>% paste0(.,'%')
    p <- ggdonutchart(
      data = df,
      x = 'Frequency',
      label = 'label',
      lab.pos = 'out',
      fill = "Feature",
      color = "black",
      lab.font = c(6, 'bold', 'black'),
      ggtheme = theme_pubr(),
      size = .75,
      pal=c(pal_tron()(4),pal_startrek()(5))
    ) +
      annotate(
        'text',
        x = 1,
        y = 0,
        label = ct,
        vjust = 3,
        size = 7
      )+theme(legend.position="right")
    plot(p)
  }
  
  dev.off()
}

# 211206
#### Greenleaf K562 10x scATAC-seq group coverage ####
# The maximum of cell number per cell line in long ATAC-seq is 1000
# So we downsample the data from Granja, et al., to 1000
if (F) {
  getGroupBW(
    proj.granja.down,
    groupBy = "cell.line",
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e3
  )
}

# 211208
# save greenleaf 10x cell line peaks
# (have down sampled to the same cell numbers as our long atac-seq)
if (F) {
  for (i in names(tenx.peak.list)){
    gr<-tenx.peak.list[[i]]
    writeGrangeAsBed(
      gr=gr,
      name = paste0('./data/granja_2021/CellLine/PeakCalls/',i,'.bed'))
  }
}


##### benchmark test of single-cell long atac-seq against bulk short ATAC-seq ####
if (F) {
  peak.benchmark.list <- list(
    GM12878 = reduce(c(readBed( './data/ngs_peak/tanglab_bulk/NGS-GE-GM12878-X1_peaks.narrowPeak'),
                      readBed('./data/ngs_peak/tanglab_bulk/NGS-GE-GM12878-X2_peaks.narrowPeak'))),
    K562 = reduce(c(readBed('./data/ngs_peak/tanglab_bulk/AH-K562-1_peaks.narrowPeak'),
                    readBed('./data/ngs_peak/tanglab_bulk/BH-K562-2_peaks.narrowPeak'))),
    HEK293T = reduce(c(readBed('./data/ngs_peak/tanglab_bulk/AH-HEK293T-1_peaks.narrowPeak'),
                       readBed('./data/ngs_peak/tanglab_bulk/BH-HEK293T-2_peaks.narrowPeak'))),
    HFF1 = reduce(c(readBed('./data/ngs_peak/tanglab_bulk/AH-HFF1-1_peaks.narrowPeak'),
                    readBed('./data/ngs_peak/tanglab_bulk/BH-HFF1-2_peaks.narrowPeak')))
  )

  peak.benchmark.list <- lapply(peak.benchmark.list, cleanChr)
  
  # plot euler diagram
  pdf('220628_peak_venn_nanoATAC_NGS.pdf',width = 5,height = 5)
  
  for (i in names(peak.benchmark.list)){
    tgs <- tgs.peak.list[[i]]
    ngs <- tenx.peak.list[[i]]
    bench <- peak.benchmark.list[[i]]
    
    if (any(is.null(tgs),is.null(ngs))){
      n.nt=0
    }else{
      hits.nt <-
        findOverlaps(query = ngs,
                     subject = tgs)
      n.nt <- length(unique(hits.nt@from))
    }

    if (any(is.null(bench),is.null(ngs))){
      n.nb=0
    }else{
    hits.nb <-
      findOverlaps(query = ngs,
                   subject = bench)
    n.nb <- length(unique(hits.nb@from))
    }
    
    if (any(is.null(tgs),is.null(bench))){
      n.tb=0
    }else{
      hits.tb <-
        findOverlaps(query = tgs,
                     subject = bench)
      n.tb <- length(unique(hits.tb@from))
    }

    if (any(is.null(tgs),is.null(bench),is.null(ngs))){
      n.share=0
    }else{
      hits.share <-
        findOverlaps(query = ngs[unique(hits.nt@from)], subject = bench)
      n.share <- length(unique(hits.share@from))
    }

    eu <- euler(
      c(
        "tgs&10x&bulk" = n.share,
        "tgs&10x" = n.nt - n.share,
        "10x&bulk" = n.nb - n.share,
        "tgs&bulk" = n.tb - n.share,
        "tgs" = length(tgs) - n.nt - n.tb + n.share,
        "10x" = length(ngs) - n.nb - n.nt + n.share,
        "bulk" = length(bench) - n.nb - n.tb + n.share
      )
    )
    
    p<-plot(
      eu,
      fills = c('steelblue1', 'orange','purple'),
      alpha = 0.7,
      labels = list(cex = .8),
      quantities = list(cex = c(0.8), type = c("counts","percent")),
      main = i
    )
    
    plot(p)
  }
  
  dev.off()
  
}

#### Peak confidence correlation #####
# score is -log10(p.val.adj)
if (F) {
  peak.benchmark.list <- list(
    GM12878 = readBed('./data/ngs_peak/tanglab_bulk/NGS-GE-GM12878-X2_peaks.narrowPeak'),
    K562 = readBed('./data/ngs_peak/tanglab_bulk/AH-K562-1_peaks.narrowPeak'),
    HEK293T = readBed('./data/ngs_peak/tanglab_bulk/AH-HEK293T-1_peaks.narrowPeak')
  )
  
  # 10x ~ bulk ATAC
  lapply(c('GM12878','K562','HEK293T'), function(i){
    subject <- peak.benchmark.list[[i]]
    seqlevelsStyle(subject) <- 'UCSC'
    
    query <- tenx.peak.list[[i]]
    hits <- findOverlaps(query, subject)
    cor.10x <- cor(query[hits@from]$score, subject[hits@to]$score,method = 'spearman')
    
    query <- tgs.peak.list[[i]]
    hits <- findOverlaps(query, subject)
    cor.nano <- cor(query[hits@from]$score, subject[hits@to]$score,method = 'spearman')
    
    data.frame(cell.type = i,
               cor.10x = cor.10x,
               cor.nano = cor.nano)
  }) %>% do.call(rbind,.) %>% write.csv('220629_peak_confidence_correlation.csv',row.names = F)
} 

#### CTCF and TSS enrichment of 10x scATAC-seq ####
if (F) {
  proj.granja <- addMotifAnnotations(ArchRProj = proj.granja, motifSet = "cisbp")
  # proj.granja <- saveArchRProject(proj.granja, './data/granja_2021/CellLine/')
  motif.pos <- getPositions(proj.granja, name='Motif')

  pdf('220629_10x_CTCF_footprint_with_ribbon.pdf',width = 5,height = 7)
  
  seFoot <- getFootprints(
    ArchRProj = proj.granja,
    positions = list(CTCF = motif.pos$CTCF_177),
    groupBy = "cell.line"
  )
  
  p <- plotFootprints(
    seFoot = seFoot,
    ArchRProj = proj.granja,
    normMethod = "subtract",
    plot = F,
    smoothWindow = 10,
    pal = pal_aaas()(10)
  )
  
  for (i in 1:length(p)){
    grid::grid.newpage()
    grid::grid.draw(p[[i]],recording = F)
  }
  
  dev.off()
  
  pdf('220629_10x_TSS_enrichment_with_ribbon.pdf', width = 5, height = 5)
  
  tss.foot <- getFootprints(
    ArchRProj = proj.granja,
    positions = list(TSS = getTSS(proj.granja)),
    groupBy = "cell.line",
    flank = 2000
  )
  
  tss.plot <- plotFootprints(
    seFoot = tss.foot,
    ArchRProj = proj.granja,
    normMethod = "Subtract",
    smoothWindow = 50,
    plot = F,
    pal = pal_aaas()(10)
  )
  
  grid.newpage()
  grid.draw(tss.plot[1][[1]])
  
  dev.off()
}
