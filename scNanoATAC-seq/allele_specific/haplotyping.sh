#! /bin/bash

subset=$1
threads=$2
root_dir=$3
chr=$4

script=$root_dir/haplotyping.sh

subset_dir=$root_dir/subset/
alignment_dir=$root_dir/alignment/
merge_dir=$root_dir/merged_bam/$subset/

#whatshap_dir=$root_dir/snp_whatshap/$subset/
hapcut2_dir=$root_dir/snp_hapcut2/$subset/

tag_dir=$root_dir/tag/$subset/
fragment_dir=$root_dir/fragment/$subset/
extract_dir=$root_dir/extract/$subset/

FILT_SIZE=1000

hg38(){
  genome_chr_size=~/ref/hg38/genome_noChr_hg38.txt
  ref_genome_fa=~/ref/hg38/GRCh38.fa
  HG001_SNP=~/ref/hg38/snp/GIAB/hg001_na12878/HG001_GRCh38_GIAB_nochr_noPATMAT.vcf
}

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

# suppose cell barcodes already added into bam files with CB flag
merge_bam(){
    init_dir $merge_dir

    bam_list=`cat $subset_dir/${subset}_subset.txt \
    | awk \
        '{
          split($0,x,"#")
          lib=x[1]
          cell=x[2]
          file=file" "dir"/"lib"/"lib"."cell".barcode.bam"
        }END{
          print file
        }' \
      dir="$alignment_dir"`

    samtools merge \
      -@ $threads \
      $subset.sorted.bam \
      $bam_list

    samtools index \
      -@ $threads \
      $subset.sorted.bam
}

hapcut2_phase(){
  init_dir $hapcut2_dir

  cat \
    <(cat $HG001_SNP | grep "^#") \
    <(cat $HG001_SNP | grep "^$chr[^0-9]" | sed 's/|/\//g') \
  > $chr.snp.vcf

  extractHAIRS \
    --bam $merge_dir/$subset.sorted.bam \
    --VCF $chr.snp.vcf \
    --ONT 1 \
    --out $subset.$chr.frag \
    --region $chr \
    --ref $ref_genome_fa
    
  HAPCUT2 \
    --error_analysis_mode 1 \
    --fragments $subset.$chr.frag \
    --VCF $chr.snp.vcf \
    --output $subset.$chr

  bgzip $subset.$chr.phased.VCF
  tabix $subset.$chr.phased.VCF.gz

  rm -f $chr.snp.vcf $subset.$chr.frag $subset.$chr
}

hapcut2_tag(){
  init_dir $tag_dir/hapcut2/

  whatshap haplotag \
    -o $subset.$chr.bam \
    --reference $ref_genome_fa \
    --ignore-read-groups \
    --regions $chr \
    $hapcut2_dir/$subset.$chr.phased.VCF.gz \
    $merge_dir/$subset.sorted.bam

  samtools index \
    -@ $threads \
    $subset.$chr.bam
}

benchmark_tag(){
  init_dir $tag_dir/giab_benchmark/

  whatshap haplotag \
    -o $subset.$chr.bam \
    --reference $ref_genome_fa \
    --ignore-read-groups \
    --regions $chr \
    ${HG001_SNP}.gz \
    $merge_dir/$subset.sorted.bam

  samtools index \
    -@ $threads \
    $subset.$chr.bam
}

# convert phased bam files to fragment bed files
bam2bed_fragment(){
  init_dir $fragment_dir/$strategy/

  if [ ! -s $subset.$chr.bed ]
  then
  bam=$tag_dir/$strategy/$subset.$chr.bam
  cat \
    <(samtools view -H $bam) \
    <(samtools view $bam) \
  | sam2bed \
  | awk -vOFS='\t' \
    '{
      if ($3-$2>fs){
        match($0,/CB:Z:(\S+)/,CB)
        match($0,/PS:i:([0-9]+)/,PS)
        match($0,/HP:i:([0-9])/,HP)
        cb=CB[1]
        ps=PS[1]
        hp=HP[1]
        if (ps != "" && hp != "" && cb != ""){
          print $1,$2,$3,$4,cb,ps,hp
        }
      }
    }' \
    fs=$FILT_SIZE \
  > $subset.$chr.bed
  fi
}

# flank with side bias
flank_fragment(){
  cd $fragment_dir/$strategy/

  if [ ! -s $subset.$chr.flank.bed ]
  then
    cat \
      <(cat $subset.$chr.bed \
        | bedtools flank \
            -i - \
            -r 1 \
            -l 0 \
            -g $genome_chr_size \
        | awk -vOFS='\t' '{$6="+\t"$6; print}') \
      <(cat $subset.$chr.bed \
        | bedtools flank \
            -i - \
            -r 0 \
            -l 1 \
            -g $genome_chr_size \
        | awk -vOFS='\t' '{$6="-\t"$6; print}') \
    > $subset.$chr.flank.bed
  fi
}


hg38

merge_bam

strategy=giab_benchmark
benchmark_tag
bam2bed_fragment
flank_fragment

# --- 220606 rebuttal protocol ---
strategy=hapcut2
hapcut2_phase
hapcut2_tag
bam2bed_fragment
flank_fragment
