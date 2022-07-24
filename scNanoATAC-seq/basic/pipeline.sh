#! /bin/bash

library=$1
sample=$2
threads=$3
root_dir=$4

outer_list=$5
inner_list=`seq 1 48`

script=$root_dir/pipeline.sh
dual_bc_script=~/ref/barcode/nanopore_dual/barcode.sh
archr_script=$root_dir/archr.R

raw_dir=$root_dir/raw.data/$library/
barcode_dir=$root_dir/barcode/$library/
trim_dir=$root_dir/trim/$library/
fragment_dir=$root_dir/fragment/$library/
coverage_dir=$root_dir/coverage/
archr_dir=$root_dir/archr/$library/

FILT_SIZE=1000

run(){
  hg38
  #mm10

  load_dual
  dual_barcode
  
  trim
  align
  sort_index
  bam_rmdup
  bam2bed_fragment
  flank_fragment
  add_barcode

  archr_fragment_file
  create_arrow_file

  call_cutesv 1
}

run_cross(){
  hg38_mm10
  
  load_dual
  dual_barcode

  trim
  align_cross
  sort_index
  bam_rmdup
}

# human reference
hg38(){
  genome_ref=~/ref/hg38/minimap2/GRCh38_ONT.mmi
  genome_chr_size=~/ref/hg38/genome_noChr_hg38.txt
  ref_genome_fa=~/ref/hg38/GRCh38.fa
  g_size=hs
  genome_name=hg38

  alignment_dir=$root_dir/alignment/$library/
}

# mouse reference
mm10(){
  genome_ref=~/ref/mm10/minimap2/mm10_ONT.mmi
  genome_chr_size=~/ref/mm10/mm10.chrom.sizes
  ref_genome_fa=~/mm10/mm10.fa
  g_size=mm
  genome_name=mm10

  alignment_dir=$root_dir/alignment_mm10/$library/
}

# cross reference
hg38_mm10(){
  genome_ref=~/ref/cross/minimap2/hg38_mm10_mix_ONT.I10G.mmi
  ref_genome_fa=~/ref/cross/hg38_mm10_mix.fa

  alignment_dir=$root_dir/alignment_cross/$library/
}

# 21.8.1
# dual barcode generator
load_dual(){
  outer_barcode=$(mktemp)
  bash $dual_bc_script "$outer_list" > $outer_barcode

  inner_barcode=$(mktemp)
  bash $dual_bc_script "$inner_list" > $inner_barcode
}

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

dual_barcode(){
  init_dir $barcode_dir

  # outer barcode demultiplexing
  nanoplexer \
    -b $outer_barcode \
    -t $threads \
    -p . \
    $raw_dir/$library.pass.f*q.gz

  # inner barcode demultiplexing
  for bc1 in $outer_list
  do
    mkdir $bc1
    nanoplexer \
      -b $inner_barcode \
      -t $threads \
      -p ./$bc1/ \
      $bc1.fastq
    rm -f $bc1.fastq
  done

  # clean up temporal files
  for bc1 in $outer_list
  do
    for bc2 in $inner_list
    do
      mv \
        ./$bc1/$bc2.fastq \
        ${bc1}_${bc2}.fastq
    done
    rm -rf ./$bc1/
  done
}

# trim adaptor
# library structure from 21.04
# symmetirc structure
trim(){
  init_dir $trim_dir

  if [ ! -s $library.$sample.trimed.fq.gz ]
  then
  cutadapt \
    -g TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG \
    -a CTGTCTCTTATACACATCTCCGAGCCCACGAGA \
    --times 1 \
    --minimum-length $FILT_SIZE \
    -e 0.2 \
    -O 12 \
    -j $threads \
    $barcode_dir/$sample.fastq \
    2> ./log/$library.$sample.txt \
  | gzip \
  > $library.$sample.trimed.fq.gz
  fi
}

align(){
  init_dir $alignment_dir
  
  if [ ! -s $library.$sample.mapQ30.bam ] || \
     [[ `wc -c $library.$sample.mapQ30.bam | grep -o ^[0-9]*` -lt 30 ]] 
  then
    minimap2 \
      --MD \
      -ax map-ont \
      -t $threads \
      $genome_ref \
      $trim_dir/$library.$sample.trimed.fq.gz \
      2> ./log/$library.$sample.log \
    | samtools view -bS -q 30 - \
    > $library.$sample.mapQ30.bam 
  fi
}

align_cross(){
  init_dir $alignment_dir
  
  if [ ! -s $library.$sample.mapQ30.bam ] || \
     [[ `wc -c $library.$sample.mapQ30.bam | grep -o ^[0-9]*` -lt 30 ]] 
  then
    minimap2 \
      --MD \
      -ax map-ont \
      -t $threads \
      $genome_ref \
      $trim_dir/$library.$sample.trimed.fq.gz \
      2> ./log/$library.$sample.log \
    | samtools view -bS -q 30 - \
    > $library.$sample.mapQ30.bam 
  fi
}

sort_index(){
  cd $alignment_dir
  
  if [ ! -s $library.$sample.mapQ30.sorted.bam ] && \
     [ -s $library.$sample.mapQ30.bam ]
  then
    samtools sort \
      -@ $threads \
      $library.$sample.mapQ30.bam \
      -o $library.$sample.mapQ30.sorted.bam
  fi

  if [ ! -s $library.$sample.mapQ30.sorted.bam.bai ]
  then
    samtools index \
      -@ $threads \
      $library.$sample.mapQ30.sorted.bam
  fi
}

bam_rmdup(){
  cd $alignment_dir

  if [ ! -s $library.$sample.mapQ30.rmdup.sorted.bam ] && \
     [ -s $library.$sample.mapQ30.sorted.bam ]
  then
    samtools rmdup -s \
      $library.$sample.mapQ30.sorted.bam \
      $library.$sample.mapQ30.rmdup.sorted.bam
    
    samtools index \
      -@ $threads \
      $library.$sample.mapQ30.rmdup.sorted.bam
  fi
}

bam2bed_fragment(){
  init_dir $fragment_dir
  
  if [ ! -s $library.$sample.bed ]
  then
    bedtools bamtobed \
      -i $alignment_dir/$library.$sample.mapQ30.rmdup.sorted.bam \
    | awk -vOFS='\t' \
      '{ if($3-$2>=fs) print }' \
      fs=$FILT_SIZE \
    | sort -k1,1 -k2,2n \
    > $library.$sample.bed
  fi
}

# flank fragment with read end bias
flank_fragment(){
  cd $fragment_dir

  if [ ! -s $library.$sample.flank.bed ]
  then
    cat \
    <(cat $library.$sample.bed \
      | grep -E "^[0-9]|^X|^Y|^chr[0-9]|^chrX|^chrY" \
      | bedtools flank \
          -i - \
          -r 1 \
          -l 0 \
          -g $genome_chr_size \
      | awk -vOFS='\t' '{$6="+"; print}') \
    <(cat $library.$sample.bed \
      | grep -E "^[0-9]|^X|^Y|^chr[0-9]|^chrX|^chrY" \
      | bedtools flank \
          -i - \
          -r 0 \
          -l 1 \
          -g $genome_chr_size \
      | awk -vOFS='\t' '{$6="-"; print}') \
    > $library.$sample.flank.bed
  fi
}

# run by library
fragment_bw(){
  init_dir $coverage_dir

  cat $fragment_dir/$library.*.flank.bed \
  | bedtools sort -i - \
  | bedtools genomecov \
      -i - \
      -bga \
      -g $genome_chr_size \
  | sort -k 1,1 -k 2,2n \
  > $library.flank.bedGraph

  bedGraphToBigWig \
    $library.flank.bedGraph \
    $genome_chr_size \
    $library.flank.bw
}

# run by single cell
call_cutesv(){
  local r=$1
  
  cutesv_dir=$root_dir/cutesv/${r}r/$library/
  init_dir $cutesv_dir
  
  if [ ! -s $library.$sample.cuteSV.ONT.vcf ] ; then
    echo "call sv $library.$sample"
    init_dir $cutesv_dir/$sample/
    cuteSV \
      $alignment_dir/$library.$sample.mapQ30.sorted.bam \
      $ref_genome_fa                          \
      $library.$sample.cuteSV.ONT.vcf         \
      .                                       \
      -t $threads                             \
      --max_cluster_bias_INS	100             \
      --diff_ratio_merging_INS	0.3           \
      --max_cluster_bias_DEL  	100           \
      --diff_ratio_merging_DEL	0.3           \
      --min_support $r                        \
      --sample ${library}.${sample} 
    
    cd $cutesv_dir
    mv ./$sample/$library.$sample.cuteSV.ONT.vcf ./
    rm -rf ./$sample/

    cat $library.$sample.cuteSV.ONT.vcf \
    | sed 's/chr//g' \
    > $library.$sample.cuteSV.ONT.vcf.clean

    mv \
      $library.$sample.cuteSV.ONT.vcf.clean \
      $library.$sample.cuteSV.ONT.vcf
  fi
}

# run by library
archr_fragment_file(){
  cd $(dirname $fragment_dir)

  cell_list=`ls ./$library/ | grep -v log | \
	awk -F "." '{print $2}' | sort | uniq`

  for cell in $cell_list
  do
    cat ./$library/$library.$cell.flank.bed \
    | awk -vOFS='\t' \
      '{
        if ($1!~/chr/) $1="chr"$1
        print $1,$2,$3,cell,1
      }' \
      cell=$cell
  done \
  | sort -k 1,1 -k 2,2n -k 4,4 \
  > $library.fragments.sorted.bed

  bgzip $library.fragments.sorted.bed
  tabix $library.fragments.sorted.bed.gz
}

# run by library
# activate r40 first
create_arrow_file(){
  source activate r40  

  init_dir $archr_dir
  
  Rscript $archr_script \
    $library $(dirname $fragment_dir) \
    $genome_name 
  
  mv $library.arrow ../
}

# run by single cell
add_barcode(){
  cd $alignment_dir

  if [[ -s $library.$sample.mapQ30.rmdup.sorted.bam ]] && \
     ([[ ! -s $library.$sample.barcode.bam ]] || \
      [[ `wc -c $library.$sample.barcode.bam | grep -o ^[0-9]*` -lt 30 ]])
  then
    cat <( samtools view $library.$sample.mapQ30.rmdup.sorted.bam -H ) \
        <( samtools view $library.$sample.mapQ30.rmdup.sorted.bam \
           | awk \
              -vOFS='\t' \
              '{ print $0,"CB:Z:"cell }' \
              cell="${library}#${sample}" ) \
    | samtools view -bS - \
    > $library.$sample.barcode.bam
  fi

  if [[ ! -s $library.$sample.barcode.bam.bai ]]
  then
    samtools index \
      -@ $threads \
      $library.$sample.barcode.bam
  fi
}

run
#run_cross