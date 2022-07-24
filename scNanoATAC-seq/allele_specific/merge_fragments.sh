#! /bin/bash

subset=$1

root_dir=~/hyq_atac/tgs/
analysis_dir=$root_dir/analysis/co_accessibility/
subset_dir=$root_dir/subset/
fragment_dir=$root_dir/fragment/
merged_fragment_dir=$analysis_dir/merged_fragment/

script=$analysis_dir/merge_fragments.sh

genome_chr_size=~/ref/hg38/genome_noChr_hg38.txt

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

merge_frag(){
    init_dir $merged_fragment_dir

    frag_list=`cat $subset_dir/${subset}_subset.txt \
    | awk \
        '{
          split($0,x,"#")
          lib=x[1]
          cell=x[2]
          file=file" "dir"/"lib"/"lib"."cell".flank.bed"
        }END{
          print file
        }' \
      dir="$fragment_dir"`

    cat $frag_list \
    | sort -k 1,1 -k 2,2n \
    | gzip \
    > $subset.flanked.bed.gz
}

merge_frag