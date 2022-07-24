args = commandArgs(trailingOnly=T)
lib = args[1]
fragment.dir = args[2]
genome = args[3]

library(ArchR)

print(lib)
print(fragment.dir)
print(genome)

# genome initiation
# e.g genome='hg19'
if (F){
  options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
  options(timeout=36000)
  addArchRGenome(genome)
}

addArchRGenome(genome)
addArchRThreads(threads = 1)

# modified minTSS for TGS atac-seq
createArrowFiles(
  inputFiles = file.path(fragment.dir,
                paste0(lib,'.fragments.sorted.bed.gz')),
  minTSS = 1,
  maxFrags = 1e+06,
  sampleNames = lib,
  subThreading = F
)


