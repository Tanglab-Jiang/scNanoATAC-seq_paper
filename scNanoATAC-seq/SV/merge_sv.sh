#! /bin/bash

sample=$1
subset=$2
raw_dir=$3
r=$4
c=$5
root_dir=$6

merge_dir=$root_dir/merge/$sample/
script=$root_dir/merge_sv.sh

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

init_dir $merge_dir

sample_file=$sample.$r.$c.txt

get_subset(){
  if [[ $subset == "." ]]
  then
    ls $raw_dir/${r}r/$sample/*.vcf > $sample_file
  else
    if [ -f $sample_file ]
    then
      rm -f $sample_file
    fi
    touch $sample_file
    
    for cell in `cat $subset | tr -s '\n'`
    do
      lib=`echo $cell | cut -d# -f1`
      id=`echo $cell | cut -d# -f2`
      echo \
        "$raw_dir/${r}r/$lib/$lib.$id.cuteSV.ONT.vcf" \
        >> $sample_file
    done
  fi
}

run_survivor(){
  max_dist=500
  min_len=50
  distinguish_type=1
  distinguish_strand=1
  estimate_distance=1

  SURVIVOR merge \
    $sample_file \
    $max_dist \
    $c \
    $distinguish_type \
    $distinguish_strand \
    $estimate_distance \
    $min_len \
    ${sample}_${r}r_${c}c.vcf
  
  rm -f $sample_file
}

get_subset
run_survivor
