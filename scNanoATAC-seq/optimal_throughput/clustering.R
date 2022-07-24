args = commandArgs(trailingOnly=T)

throughput = args[1]
seed = args[2]
archr.dir = args[3]

library(ArchR)
library(stringr)

addArchRThreads(1)
addArchRGenome("hg38")

setwd(archr.dir)

ct<-c("K562","HEK293T","HFF1","GM12878","eHAP1")

obj <- ArchRProject(
    ArrowFiles = paste(throughput,ct,seed,'arrow',sep='.'),
    outputDirectory = './',
    copyArrows = F
)

obj$cell.line<-sapply(str_split(obj$Sample,'\\.'),function(x)x[2])

obj$cor_off<-runif(length(obj$cellNames))

obj <- addIterativeLSI(
    ArchRProj = obj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    dimsToUse = 1:30,
    corCutOff = 1,
    depthCol = 'cor_off'
)

# cluster cells into 5 cell lines
obj <- addClusters(
    input = obj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    maxClusters = 5
)

tb <- table(obj$Clusters, obj$cell.line)
saveRDS(tb, paste(throughput,seed,'clustering','Rds',sep='.'))
