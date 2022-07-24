args = commandArgs(trailingOnly=T)

ct = args[1]
throughput = args[2]
seed = args[3]
fragment.dir = args[4]

library(ArchR)

addArchRGenome('hg38')
addArchRThreads(threads = 1)

# modified minTSS for TGS atac-seq
createArrowFiles(
  inputFiles = paste0(
    fragment.dir,'/',throughput,'/',
    ct,'.downsampled.seed_',seed,'.fragments.sorted.bed.gz'),
  minTSS = 1,
  maxFrags = 1e+06,
  sampleNames = paste0(throughput,'.',ct,'.',seed),
  subThreading = F
)
