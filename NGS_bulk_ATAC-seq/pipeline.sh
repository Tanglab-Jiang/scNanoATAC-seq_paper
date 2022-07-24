#! /bin/bash

lib=$1
sample=$2
threads=$3
root_dir=$4

script=$root_dir/pipeline.sh

raw_dir=$root_dir/raw.data/$lib/$sample/
trim_dir=$root_dir/trim/$lib/
align_dir=$root_dir/alignment/$lib/
peak_dir=$root_dir/peak/$lib/

EXT_SIZE=150

hg38(){
  genome_name=hg38
  g_size=hs
  ref=~/ref/hg38/bowtie2/hg38
  chr_size=~/ref/hg38/genome_hg38.txt
  genome_fa=~/ref/hg38/GRCh38.fa
}

mm10(){
  genome_name=mm10
  g_size=mm
  ref=~/ref/mm10/bowtie2/mm10
  chr_size=~/ref/mm10/mm10.chrom.sizes
}

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

trim(){
  init_dir $trim_dir
  if [ ! -f $align_dir/$sample.bam ]; then
  cutadapt \
    -q 20 \
    -g AGATGTGTATAAGAGACAG \
    -G AGATGTGTATAAGAGACAG \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    --times 1 \
    --minimum-length 30 \
    -j $threads \
    -o $sample.trimed.R1.fq.gz \
    -p $sample.trimed.R2.fq.gz \
    $raw_dir/$sample.R1.fq.gz \
    $raw_dir/$sample.R2.fq.gz \
    > ./log/$sample.txt
  fi
}

align(){
  init_dir $align_dir
  
  if [ ! -f $sample.bam ]; then 
  bowtie2 \
    -p $threads \
    --mm --no-unal \
    -x $ref \
    -1 $trim_dir/$sample.trimed.R1.fq.gz \
    -2 $trim_dir/$sample.trimed.R2.fq.gz \
  2> ./log/$sample.log \
  | samtools view -bS -q 30 - \
  > $sample.mapQ30.bam
  fi
}

sort_index(){
  cd $align_dir
  
  if [ ! -f $sample.mapQ30.sorted.bam ]; then
    samtools sort \
        -@ $threads \
        $sample.mapQ30.bam \
        -o $sample.mapQ30.sorted.bam
  fi

  if [ ! -f $sample.mapQ30.sorted.bam.bai ]; then
    samtools index \
        -@ $threads \
        $sample.mapQ30.sorted.bam
  fi
}

dedup(){
  cd $align_dir
  
  if [ ! -f $sample.mapQ30.sorted.dedup.bam ];then
    picard MarkDuplicates \
        I=$sample.mapQ30.sorted.bam \
        O=$sample.mapQ30.sorted.dedup.bam  \
        M=./log/$sample.dedup.metrics \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true
  fi
}

call_narrow_peak(){
  init_dir $peak_dir

  macs2 callpeak \
    -t $align_dir/$sample.mapQ30.sorted.dedup.bam \
    -n $sample \
    -f BAM \
    -g $g_size \
    -q 0.01 \
    --nomodel \
    --shift -$(($EXT_SIZE/2)) \
    --extsize $EXT_SIZE \
    --bdg \
    --outdir . \
  2> ./log/$sample.log

  rm -f ${sample}_control_lambda.bdg
}

bg2bw(){
  cd $peak_dir
  
  if [ ! -f ${sample}.bw ]
  then
    bg_temp=$(mktemp)
    chr_temp=$(mktemp)

    cat ${sample}_treat_pileup.bdg \
    | sed 's/^chr//g' \
    | grep -E "^[0-9]|^X|^Y" \
    > $bg_temp

    cat $chr_size \
    | sed 's/^chr//g' \
    > $chr_temp

    bedGraphToBigWig \
      $bg_temp \
      $chr_temp \
      ${sample}.bw

    rm $bg_temp $chr_temp
  fi
}

hg38
#mm10

trim
align
sort_index
dedup
call_narrow_peak
bg2bw