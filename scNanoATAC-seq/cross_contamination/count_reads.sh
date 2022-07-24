#! /bin/bash

library=$1
sample=$2
threads=$3
root_dir=$4

alignment_dir=$root_dir/alignment_cross/$library/
analysis_dir=$root_dir/analysis/cross_contamination/
count_dir=$analysis_dir/count/$library/

init_dir(){
  mkdir -p $1
  cd $1
}

# by cell
count_cross(){
    init_dir $count_dir

    samtools view \
        $alignment_dir/$library.$sample.mapQ30.rmdup.sorted.bam \
    | awk -vOFS='\t' \
      'BEGIN{hg=0;mm=0}
      {
         if($3~/hg38/){ hg+=1 }else if($3~/mm10/){ mm+=1 }
      }
      END{print hg,mm}' \
    > $sample.txt
}

# run one
merge_count(){
    cd $(dirname $count_dir)

    for i in `ls $count_dir`
    do
        sample=$(echo $i | sed 's/.txt//g')
        cat $count_dir/$i \
        | awk -vOFS='\t' \
            '{ print sample,$0}' \
	  sample=$sample
    done > $library.cross.txt
}

count_cross
merge_count
