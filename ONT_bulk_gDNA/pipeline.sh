#! /bin/bash

sample=$1
threads=$2
root_dir=$3

script=$root_dir/pipeline.sh

raw_dir=$root_dir/raw.data/$sample/
alignment_dir=$root_dir/alignment/$sample/
coverage_dir=$root_dir/coverage/$sample/

hg38(){
  genome_ref=~/ref/hg38/minimap2/GRCh38_ONT.mmi
  ref_genome_fa=~/ref/hg38/GRCh38.fa
}

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

align(){
  init_dir $alignment_dir
  
  if [ ! -f $sample.mapQ30.bam ] || \
     [ ! -s $sample.mapQ30.bam ] || \
     [[ `wc -c $sample.mapQ30.bam | grep -o ^[0-9]*` -lt 30 ]]
   then
    echo "align $sample"
    minimap2 \
      --MD \
      -ax map-ont \
      -t $threads \
      $genome_ref \
      $raw_dir/$sample.*f*q.gz \
    2> ./log/$sample.log \
    | samtools view -bS -q 30 - \
    > $sample.mapQ30.bam 
  fi
}

# run on fat nodes
sort_index(){
  cd $alignment_dir
  
  if [ ! -f $sample.mapQ30.sorted.bam ] || \
     [ ! -s $sample.mapQ30.sorted.bam ] || \
     [ ! -f $sample.mapQ30.sorted.bam.bai ] ; then
    echo "sort $sample"
    samtools sort \
      -@ $threads \
      -m 2G \
      -o $sample.mapQ30.sorted.bam \
      $sample.mapQ30.bam
    samtools index \
      -@ $threads \
      $sample.mapQ30.sorted.bam
  fi
}

call_cutesv(){
  local r=$1
  
  cutesv_dir=$root_dir/cutesv/$sample/
  
  if [ ! -f $cutesv_dir/$sample.cuteSV.${r}r.ONT.vcf ] || \
    [ ! -s $cutesv_dir/$sample.cuteSV.${r}r.ONT.vcf ] ; then
    echo "call $sample sv"
    init_dir $cutesv_dir
    cuteSV \
      $alignment_dir/$sample.mapQ30.sorted.bam \
      $ref_genome_fa                          \
      $sample.cuteSV.${r}r.ONT.vcf            \
      .                                       \
      -t $threads                             \
      --max_cluster_bias_INS	100             \
      --diff_ratio_merging_INS	0.3           \
      --max_cluster_bias_DEL  	100           \
      --diff_ratio_merging_DEL	0.3           \
      --min_support $r                        \
      --sample ${sample}

    cat $sample.cuteSV.${r}r.ONT.vcf \
    | sed 's/chr//g' \
    > $sample.cuteSV.${r}r.ONT.vcf.clean
  
    mv $sample.cuteSV.${r}r.ONT.vcf.clean \
      $sample.cuteSV.${r}r.ONT.vcf
  fi
}

hg38

align
sort_index
call_cutesv 10

