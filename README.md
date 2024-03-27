# Scripts for scNanoATAC-seq
`scNanoATAC-seq` is a long-read sequencing (Oxford Nanopore Technologies) based single-cell ATAC-seq method designed to capture more genetic variations per read. Our results demonstrated that scNanoATAC-seq can detect chromatin accessibility and genetic variants simultaneously within an individual cell, and it overcomes the disadvantages of short-read ATAC-seq in haplotyping and the detection of large-scale structural variations (SVs).

# Code Layout
Generally, we used shell scripts to run time-consuming pipelines, which implemented parallel computation on high-performance clusters. After that, we made exploratory data analyses locally in the R environment.

In this project, each folder stands for a module of analysis. If names of scripts start with step numbers (e.g. `s1` and `s2` stands for the first and second step of a pipeline), please **follow the step numbers and run them in sequential**. For those without step numbers, the scripts are parallel and independent.

## scNanoATAC-seq
These scripts in the `scNanoATAC-seq` folder produced the main results of our work, which analyzed scNanoATAC-seq data sequenced from human cell lines, clinical samples (PBMCs) and mouse cell lines.

### **Basic**
This shell pipeline for sequence data preprocessing was in the `basic` folder. The preprocessing converts the raw sequence data to chromatin accessibility signals in the ArchR project.

**Update (Mar 2024):**
**There was an issue with the adapter sequences provided in the Methods section of our paper. We have updated them in the `basic.sh` script.**

*Input*: 

library name (raw data); single-cell ID; threads; working directory; the list of outer barcode.

*Processing*: 

barcode demultiplexing; read trimming; alignment; read quality filtering; sorting; removing PCR duplicates; extraction of read end coordinates; generation of the `ArchR` project; SV detection at the single-cell level (by `cuteSV`); alignment to human-mouse mixed reference genome.

*Output*: 

trimmed and aligned reads in bam files by single cell; an `ArchR` project; single-cell SV; human-mouse mixed alignment.


### **Main**
The `main` folder contains the R scripts used to analyze and illustrate single-cell data produced by scNanoATAC-seq, where most of the analyses were based on `ArchR`. Analyses of human cell lines, human PBMCs and mouse cell lines were in separated scripts.

*Functions*: 

annotation of metadata; quality control; dimensional reduction; unsupervised clustering; peak calling; genome browser style track plot; feature plot; TSS enrichment; TF motif footprint; benchmark test of peak calling on public datasets (including single cell downsampling, peak calling, cCREs benchmark test and bulk ATAC-seq benchmark test); peak annotations.

### **Quality Control**
The pipeline for sequence data quality control is in the `quality_control` folder, which calculates and summarizes metrics at both library and single-cell level.

*Input*: 

library name (log and bam files from `basic`); threads; working directory.

*Processing*: 

generate the summary of the quality of a library; quality control for trimming; quality control for alignment; quality control by genome coverage; merge metrics by single cell; merge metrics by library.

*Output*: 

library quality control table; single-cell quality control table.

### **Cross-contamination**
Scripts specialized for scNanoATAC-seq human-mouse cross-contamination libraries were in the `cross_contamination` folder.

*Input*:

library name; single-cell ID (bam files from `basic`); threads; working directory.

*Processing*:

count reads by reference genome (hg38 or mm10) in each single cell; merge counts by single cell; calculation and illustration in R.

*Output*:

single-cell read counting table; cross-contamination statistics and plots.

### **Optimal Throughput**
The estimation of optimal throughput of scNanoATAC-seq was based on downsampled simulated single-cell data. Scripts are in the `optimal_throughput` folder.

*Functions*:

Simulating single-cell data from scNanoATAC-seq data of human cell lines. We downsampled reads from pseudo bulk of a cell line without replacement at a specific throughput (cells per run). You can set series of random seeds in this step for the reproducibility. Then simulated data of five cell lines at the same condition were combined and put into `ArchR` for unsupervised clustering. The precision of clustering was calculated by batch of simulation and used to evaluating the technical performance under a condition.

### **Allele-specific Accessibility**
Analysis of allele-specific accessibility based on scNanoATAC-seq were performed in the `allele_specific` folder. The cell name list in `subset` folder under this branch recorded single-cell ID belonging to each cell line.

*Functions*:

1. [s1] Merge single-cell bam files by cell line. (assuming GM12878 is analyzed)

2. [s2] Phase prior known heterozygous SNP of GM12878 by `HAPCUT2`. Otherwise, you can download phased SNP of GM12878 from `GIAB` database.

3. [s2] Genotyping (or haplotyping) the reads by `whatshap haplotag` function based on phased SNPs.

4. [s2] Extract chromatin accessibility signals at the both ends of reads like what we did in `basic` part, meanwhile the genotype (or haplotype) was annotated in each signal record.

5. [s3] In R environment, we took the subset of the signals inside peaks called from GM12878. Then, allelic bias of accessibility was tested by peak with `binom.test` (H*0*: π=0.5). P values were adjusted after multiple hypothesis testing by the `fdr` method. Allele-specific peaks (ASPs) with adjusted P values less than 0.05 were kept.

6. [s4] We validated those ASPs detected by scNanoATAC-seq using NGS based bulk ATAC-seq data of GM12878. Heterozygous SNP of GM12878 were intersected with ASPs and their allele frequency bias (derived from the `Short-read Bulk ATAC-seq` part) was tested by `binom.test` (H*0*: π=0.5) with P value adjusting. scNanoATAC-seq ASPs with any significantly biased heterozygous SNP detected by bulk ATAC-seq were regarded as true positives.

7. [s5] (Optional) If you phased SNPs based on long reads of scNanoATAC-seq in `s2`, you may evaluate the switch error of SNP phasing. You need to provide pedigree information and SNP sets of the trio.

8. [s6] Haplotyped chromatin accessibility of GM12878 by scNanoATAC-seq can be divided into paternal and maternal coverage tracks in `bigwig` format. This made a better visualization of allele-specific accessibility in genome browser.

### **Structural Variation (SV)**

In this section, single-cell SV files called by `basic` pipeline will be further analyzed. Scripts in the `SV` folder will merge VCF files from the single-cell level to the cell line level, and perform breakpoint analysis of SV based on `StructuralVariantAnnotation` R implementation.

*Functions*:

1. Call SURVIVOR to merge single-cell VCF files recording SVs by cell line. 

2. Fix the merged VCF file by adding unique ID to each record and sorting according to ID. In this way, we made the VCF file acceptable for `StructuralVariantAnnotation` package.

3. Load VCF files recording SVs into R environment, including benchmark SVs of K562 by bulk ONT sequencing (by `Long-read Bulk WGS` pipeline), SVs of K562 and GM12878 called by scNanoATAC-seq, germline SV datasets (HGSVC2 and gnomAD v2). Using the R pipeline, we performed SV benchmark test, SV filtering (to get somatic SVs in K562), extraction of the validation set of K562 somatic SV.
  
### **Copy Number Variation (CNV)**

The R script in the `CNV` folder is used to make single-cell copy number estimation matrix from scNanoATAC-seq. We binned the human genome by 100-kb tile and removed the tiles intersecting with the black list regions of hg38. With controlling and smoothing of single-cell CNV matrix, we finally made a heatmap to illustrate CNVs in five human cell lines at single-cell level from scNanoATAC-seq data.


### **Co-accessibility**

Scripts in the `co_accessibility` folder conducted the co-accessibility detection in a way unique to scNanoATAC-seq data. The main idea was that if two neighboring peaks were co-accessible in each single cell, the length of reads between them will change relative to the background distribution of read length. 

In the step one, we called R script to detect the alteration of read length distribution by `ks.test` function. And in the step two, we summarized the parallel computation results with multiple hypothesis test correction, co-accessible peak span filtering. Finally, we got co-accessibiltiy pairs of neighboring peaks by cell line and wrote them into bedpe files for IGV visualization.

## Short-read Bulk ATAC-seq

Validation of the technical reliability of scNanoATAC-seq was supported by short-read bulk ATAC-seq data. In this section, pipelines were used to implement basic analyses of bulk ATAC-seq data, including trimming, alignment, sorting, read quality filtering, PCR duplicate removing and peak calling. To support the validation of ASPs called by scNanoATAC-seq, we used `freebayes` to calculate allele frequency of known heterozygous SNPs of GM12878 in the context of ATAC-seq. All the scripts are in the `NGS_bulk_ATAC-seq` folder.

## Long-read Bulk WGS

To evaluate SVs detected by scNanoATAC-seq, we performed long-read sequencing on K562 sample without PCR amplification of genomic DNA, and got a long-read bulk WGS library. The pipeline in the `ONT_bulk_gDNA` folder for processing the data is relatively simple. We aligned raw reads directly to hg38 reference genome, and filtered them according to mapQ (greater than 30). The SV caller for long-read bulk WGS was `cuteSV`, which is consistent with that previously used in the `basic` pipeline of scNanoATAC-seq.

# Requirements
Install the following software before running the scripts.
```
    python (3.8.5)
    R (4.1.0)
    ArchR (1.0.1)
    nanoplexer (0.1)
    cutadapt (3.2)
    minimap2 (2.17-r941)
    bowtie2 (2.4.2)
    samtools (1.14)
    bedtools (2.30.0)
    MACS2 (2.2.7.1)
    MultiQC (1.9)
    freebayes (1.3.5)
    bcftools (1.14)
    cuteSV (1.0.10)
    SURVIVOR (1.0.7)
    HapCUT2 (1.3.1)
    whatshap (1.2.1)
    Picard
    bedGraphToBigWig

    (and all the packages declared at the head of R scripts...)
```

Please note that the generation of arrow files depends on **ArchR (v1.0.1)**, because newer versions filter the fragment by size. Newer versions of ArchR may work well on further downstream analyses after generation of arrow files.

# To Get Help
Please go to the [issues page][issue] for help or contact with me directly by zhenhuan.jiang@stu.pku.edu.cn.

[issue]: https://github.com/Tanglab-Jiang/scNanoATAC-seq_paper/issues/

# Cite Me
Hu, Y., Jiang, Z., Chen, K. et al. scNanoATAC-seq: a long-read single-cell ATAC sequencing method to detect chromatin accessibility and genetic variants simultaneously within an individual cell. Cell Res (2022). https://doi.org/10.1038/s41422-022-00730-x

# License
MIT © Zhenhuan Jiang
