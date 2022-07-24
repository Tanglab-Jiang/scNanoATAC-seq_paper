# HYQ long-read atac-seq project

library(ArchR)
library(genomation)
library(pheatmap)
library(eulerr)
library(ggsci)
library(colorspace)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggpubr)
library(ggsci)

setwd("~/SHARE/HYQ_atac/")

addArchRThreads(1)
addArchRGenome("hg38")
data.path <- './data/arrow/'

# SAVE
# saveArchRProject(ArchRProj = proj, outputDirectory = './data/projects/PBMC/')
# saveArchRProject(ArchRProj = proj.1k, outputDirectory = './data/projects/PBMC.1k/')
# saveArchRProject(ArchRProj = proj.10x.pbmc, outputDirectory = './data/10x_1k_PBMC/project/')

# LOAD
proj <- loadArchRProject('./data/projects/PBMC/')
proj.1k <- loadArchRProject('./data/projects/PBMC.1k/')
proj.10x.pbmc <- loadArchRProject('./data/10x_1k_PBMC/project/')

# resources
if (T) {
  ccres<-readRDS('./resource/210926_ENCODE_cCREs_GRanges.Rds')
}

##### load peak called by ArchR by PBMC cell type #####
if (T) {
  cell.type<-c("CD4+ T","CD8+ T","B","Monocyte")
  
  # Load TGS atac-seq peak set by cell type 
  tgs.peak.list <- lapply(cell.type, function(cl) {
    cl<-cl %>% gsub('\\+','.',.) %>% gsub(' ','.',.)
    readRDS(file.path(
      './data/projects/PBMC.1k/PeakCalls/',
      paste0(cl, '-reproduciblePeaks.gr.rds')
    ))
  })
  names(tgs.peak.list) <- cell.type
  
  # Load 10x atac-seq peak set by cell type 
  tenx.peak.list <- lapply(cell.type, function(cl) {
    cl<-cl %>% gsub('\\+','.',.) %>% gsub(' ','.',.)
    readRDS(file.path(
      './data/10x_1k_PBMC/project/PeakCalls/',
      paste0(cl, '-reproduciblePeaks.gr.rds')
    ))
  })
  names(tenx.peak.list) <- cell.type
}

#### functions ####
if (T) {
  # patches for ArchR
  source('../resource/archr_patches/archr_track.R')
  source('../resource/archr_patches/archr_footprint_no_ribbon.R')

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
  
  # 211027
  # run a benchmark test
  ccres_benchmark<-function(peak.set,group=NULL){
    result <- list()
    
    for (i in names(peak.set)) {
      query <- peak.set[[i]]
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
}

##### make long ATAC-seq PMBC project #####
if (F) {
  proj <- ArchRProject(
    ArrowFiles = paste0(data.path, c('2110_PBMC','2110_CD19','2110_CD48'),'.arrow'),
    outputDirectory = './data/projects/PBMC/',
    copyArrows = T
  )
  
  # annotation
  if (T) {
    # 2110 CD19
    info <- c(rep('CD19',10),
              rep('PBMC',4),
              rep('CD4 Median',6))
    names(info) <- 73:92
    proj <-  setCellGroup(proj, info, '2110_CD19')
    
    # 2110 CD4 & CD8
    info <- c(rep('CD4',10),
              rep('CD8',10))
    names(info) <- 73:92
    proj <-  setCellGroup(proj, info, '2110_CD48')
    
    # 2110 PBMC
    info <- c(rep('PBMC',20))
    names(info) <- 73:92
    proj <-  setCellGroup(proj, info, '2110_PBMC')
  }
  
  proj <- proj %>% subset(proj$nFrags > 1e4)
  proj <- archrReduce(proj)
  
  proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    maxClusters = 4,
    force = T
  )
  
  cell.type<-c('CD8+ T','CD4+ T','Monocyte','B')
  names(cell.type)<-paste0('C',1:4)
  proj$cell.type<-cell.type[proj$Clusters]
  
  proj <- addGroupCoverages(
    ArchRProj = proj,
    groupBy = "cell.type"
  )
  
  proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "cell.type",
    cutOff = 0.01,
    force = T,
    additionalParams = '--nomodel --llocal 500000',
    threads = 4
  )
  
  getGroupBW(
    ArchRProj = proj,
    groupBy = "cell.type",
    normMethod = "nFrags",
    tileSize = 100,
    maxCells = 1e4
  )
  
  # get marker genes & peaks
  if (F) {
    markersGS <- getMarkerFeatures(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      groupBy = "Sample",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon",
      threads = 6
    )
    
    markersPeaks <- getMarkerFeatures(
      ArchRProj = proj, 
      useMatrix = "PeakMatrix", 
      groupBy = "Clusters",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      testMethod = "wilcoxon"
    )
  }
}

#### make a project for 10x PBMC ATAC-seq benchmark ####
if (F) {
  addArchRGenome('hg38')
  
  setwd('~/SHARE/HYQ_atac/data/10x_1k_PBMC/')
  createArrowFiles(
    inputFiles = 'atac_pbmc_5k_nextgem_fragments.tsv.gz',
    minTSS = 4,
    maxFrags = 1e6,
    sampleNames = '10x_pbmc',
    subThreading = F
  )
  
  setwd('~/SHARE/HYQ_atac/')
  proj.10x.pbmc <- ArchRProject(
    ArrowFiles = './data/10x_1k_PBMC/10x_pbmc.arrow',
    outputDirectory = './data/10x_1k_PBMC/project/',
    copyArrows = T
  )
  
  proj.10x.pbmc<-archrReduce(proj.10x.pbmc)
  
  proj.10x.pbmc <- addClusters(
    input = proj.10x.pbmc,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    maxClusters = 4,
    force = T
  )
  
  cell.type<-c('Monocyte','B','CD4+ T','CD8+ T')
  names(cell.type)<-paste0('C',1:4)
  proj.10x.pbmc$cell.type<-cell.type[proj.10x.pbmc$Clusters]
  
  proj.10x.pbmc <- addGroupCoverages(
    ArchRProj = proj.10x.pbmc,
    groupBy = "cell.type",
    force = T
  )
  
  proj.10x.pbmc <- addReproduciblePeakSet(
    ArchRProj = proj.10x.pbmc, 
    groupBy = "cell.type",
    cutOff = 0.01,
    force = T
  )
  
  if (T) {
    rename<-c('CD4_T','CD8_T','B','Monocyte')
    names(rename)<-names(tenx.peak.list)
    for (i in names(tenx.peak.list)) {
      writeGrangeAsBed(
        tenx.peak.list[[i]],
        paste0('./data/10x_1k_PBMC/project/PeakCalls/', rename[i],'.bed'))
    }
  }
  
  proj.10x.pbmc <- addPeakMatrix(
    proj.10x.pbmc,
    force = T
  )
  
  saveArchRProject(ArchRProj = proj.10x.pbmc,
                   outputDirectory = './data/10x_1k_PBMC/project/')
}

##### downsample TGS PBMC to identical cell number with 10x #####
if (F) {
  proj.1k <- subsetArchRProject(
    proj,
    cells=sample(proj$cellNames,1124),
    outputDirectory='./data/projects/PBMC.1k/',
    force = T
  )
  
  proj.1k <- archrReduce(proj.1k)
  
  proj.1k <- addGroupCoverages(
    ArchRProj = proj.1k,
    groupBy = "cell.type",
    force = T
  )
  
  proj.1k <- addReproduciblePeakSet(
    ArchRProj = proj.1k, 
    groupBy = "cell.type",
    cutOff = 0.01,
    force = T,
    additionalParams = '--nomodel --llocal 500000'
  )
  
  # write peaks into bed files by cell type
  if (T) {
    rename<-c('CD4_T','CD8_T','B','Monocyte')
    names(rename)<-names(tgs.peak.list)
    for (i in names(tgs.peak.list)) {
      writeGrangeAsBed(
        tgs.peak.list[[i]],
        paste0('./data/projects/PBMC.1k/PeakCalls/', rename[i],'.bed'))
    }
  }

  proj.1k <- addPeakMatrix(proj.1k)
  
  saveArchRProject(ArchRProj = proj.1k, 
                   outputDirectory = './data/projects/PBMC.1k/')
}

# TGS PBMC data visualization
if (F) {
  # marker
  markerGenes <- 
    c('PTPRC','CD19','CD79A','MS4A1',  # B
      'CD3D','CD4','CD8A',  # T
      'CD14','ITGAM','PECAM1','LYZ',  # monocytpe
      'NCAM1','FCGR3A','FCGR3B'  # NK
    )

  pdf('211201_PBMC_new_wrapper.pdf',width = 6,height = 6)

  # plot UMAP
  plotEmbeddingWrap(proj, name = 'Sample', color = pal_lancet()(9))
  plotEmbeddingWrap(proj, name = 'group', color = pal_lancet()(9))
  plotEmbeddingWrap(proj, name = 'cell.type', color = pal_lancet()(9))
  
  dev.off()
  
  pdf('211201_PBMC_featureplot.pdf',width = 8,height = 8)
  # feature plot
    for ( i in 1:ceiling(length(markerGenes)/9) ) {
      glist <- markerGenes[(9 * (i - 1) + 1):min(length(markerGenes), 9 * i)]
      p <- plotEmbedding(
        ArchRProj = proj,
        colorBy = "GeneScoreMatrix",
        name = glist,
        embedding = 'UMAP',
        rastr = F,
        plotAs = "points",
        size = 0.25
      )
      
      if (class(p) != 'list') p<-list(p)
      
      pl <- lapply(p, function(x) {
        x +
          theme_ArchR(baseSize = 6.5) +
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
      })
      
      do.call(cowplot::plot_grid, c(list(ncol = 3, nrow = 3), pl)) %>% plot()
    }

  dev.off()
}

##### plot marker tracks for TGS PBMC #####
if (F) {
  pdf('211026_TGS_PBMC_marker_track.pdf',width = 8,height = 4)
  
  p <- plotBrowserTrack(
    ArchRProj = proj,
    groupBy = "cell.type",
    geneSymbol = markerGenes
  )
  
  for (g in markerGenes) {
    grid::grid.newpage()
    grid::grid.draw(p[g][[1]])
  }
  
  dev.off()
  
  if (F) {
    plotMarkerHeatmap(
      seMarker = markersGS,
      cutOff = "FDR <= 0.05 & Log2FC >= 1",
      labelMarkers = markerGenes,
      transpose = TRUE
    )
  }
}

##### plot marker tracks for 10x PBMC #####
if (F) {
  pdf('220629_10x_PBMC_marker_track.pdf',width = 8,height = 4)
  
  markerGenes <- 
    c('PTPRC','CD19','CD79A','MS4A1',  # B
      'CD3D','CD4','CD8A',  # T
      'CD14','ITGAM','PECAM1','LYZ'  # monocytpe
    )
  
  plotEmbedding(proj.10x.pbmc,name = 'cell.type',size = 1,rastr = F)
  
  p <- plotBrowserTrack(
    ArchRProj = proj.10x.pbmc,
    groupBy = "cell.type",
    geneSymbol = markerGenes
  )
  
  for (g in markerGenes) {
    grid::grid.newpage()
    grid::grid.draw(p[g][[1]])
  }
  
  dev.off()
  
  pdf('220603_TGS_1k_PBMC_marker_track.pdf',width = 8,height = 4)
  
  markerGenes <- 
    c('PTPRC','CD19','CD79A','MS4A1',  # B
      'CD3D','CD4','CD8A',  # T
      'CD14','ITGAM','PECAM1','LYZ'  # monocytpe
    )
  
  plotEmbedding(proj.1k,name = 'cell.type',size = 1,rastr = F)
  
  p <- plotBrowserTrack(
    ArchRProj = proj.1k,
    groupBy = "cell.type",
    geneSymbol = markerGenes
  )
  
  for (g in markerGenes) {
    grid::grid.newpage()
    grid::grid.draw(p[g][[1]])
  }
  
  dev.off()
}

##### PBMC 10x and TGS venn #####
if (F) {
  pdf('220625_PBMC_peak_benchmark.pdf',width = 5,height = 5)
  for (ct in cell.type){
    query <- tgs.peak.list[[ct]]
    subject <- tenx.peak.list[[ct]]
    hits <-
      findOverlaps(query = query, subject = subject)
    
    isec.q<-length(unique(hits@from))
    isec.s<-length(unique(hits@to))
    
    eu <- euler(c(
      "ONT&10x" = round((isec.q+isec.s)/2),
      "ONT" = length(query) - isec.q,
      "10x" = length(subject) - isec.s
    ))
    
    p <- plot(
      eu,
      fills = c('steelblue1', 'orange'),
      alpha = 0.7,
      labels = list(cex = .8),
      quantities = list(cex = c(0.8), type = c("counts", "percent")),
      main = ct
    )
    
    plot(p)
  }
  
  dev.off()
}

##### take ENCODE cCREs as a benchmark #####
if (F) {
  ccres.bench.df <- rbind(
    ccres_benchmark(tenx.peak.list, group = '10X'),
    ccres_benchmark(tgs.peak.list, group = 'TGS')
  )
  
  ccres.bench.df$group <- factor(ccres.bench.df$group, c('10X','TGS'))
  
  pdf('220625_1k_PBMC_cCREs_benchmark.pdf',width = 12,height = 4)
  
  dummy<-ccres.bench.df
  dummy$value[dummy$metric=='Precision']<-1
  dummy$value[dummy$metric=='Recall']<-.2
  dummy$value[dummy$metric=='F1']<-.2
  
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

# 211208
#### save coverage tracks ####
if (F) {
  getGroupBW(proj.10x.pbmc,'cell.type')
}

# 220106
##### cell id susbet for each cell type #####
if (F) {
  for (i in unique(proj$cell.type)) {
    name<-i %>% gsub(pattern = '\\+',replacement = '_',.) %>% gsub(pattern = '\\ ',replacement = '',.)
    proj$cellNames[proj$cell.type == i] %>% 
      sort %>% 
      write(paste0('./data/subset/', name, '_subset.txt'))
  }
}

#### 10x single-cell UMAP and featureplot ####
if (F) {
  # marker
  markerGenes <- 
    c('PTPRC','CD19','CD79A','MS4A1',  # B
      'CD3D','CD4','CD8A',  # T
      'CD14','ITGAM','PECAM1','LYZ',  # monocytpe
      'NCAM1','FCGR3A','FCGR3B'  # NK
    )
  
  pdf('220629_PBMC_10x.pdf',width = 6,height = 6)
  
  # plot UMAP
  plotEmbeddingWrap(proj.10x.pbmc, name = 'cell.type', color = pal_lancet()(9))
  
  # feature plot
  for ( i in 1:ceiling(length(markerGenes)/9) ) {
    glist <- markerGenes[(9 * (i - 1) + 1):min(length(markerGenes), 9 * i)]
    p <- plotEmbedding(
      ArchRProj = proj.10x.pbmc,
      colorBy = "GeneScoreMatrix",
      name = glist,
      embedding = 'UMAP',
      rastr = F,
      plotAs = "points",
      size = 0.25
    )
    
    if (class(p) != 'list') p<-list(p)
    
    pl <- lapply(p, function(x) {
      x +
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
    })
    
    do.call(cowplot::plot_grid, c(list(ncol = 3, nrow = 3), pl)) %>% plot()
  }
  
  dev.off()
}
