#! /bin/bash
genome_chr_size=~/ref/hg38/genome_noChr_hg38.txt

#strategy=giab_benchmark
strategy=hapcut2

root_dir=~/hyq_atac/tgs/analysis/haplotype/
mkdir $root_dir/fragment/GM12878/${strategy}_coverage/
cd $root_dir/fragment/GM12878/${strategy}_coverage/

#name=pat
#hap=1

name=mat
hap=2

cat ../${strategy}/GM12878.*.flank.bed \
| awk -vOFS='\t' '{ if ($8==hap) print}' hap=$hap \
| bedtools sort -i - \
| bedtools genomecov \
    -i - \
    -bga \
    -g $genome_chr_size \
| sort -k 1,1 -k 2,2n \
> GM12878.$name.bedGraph

bedtools makewindows \
    -g $genome_chr_size \
    -w 100 \
| sort -k 1,1 -k 2,2n \
> hg38.100bp.windows.bed

bedtools map \
    -a hg38.100bp.windows.bed \
    -b GM12878.$name.bedGraph \
    -c 4 \
    -o sum \
> GM12878.$name.100bp.bedGraph

bedGraphToBigWig \
    GM12878.$name.100bp.bedGraph \
    $genome_chr_size \
    GM12878.$name.100bp.bw
