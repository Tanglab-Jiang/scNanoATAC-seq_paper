#! /bin/sh

root_dir=~/hyq_atac/analysis/sv/survivor/
subset_dir=~/hyq_atac/tgs/subset/
raw_dir=~/hyq_atac/tgs/cutesv/

declare -A sample_list
sample_list=(\
[2110_K562_subset]="$subset_dir/2110_K562_subset.txt" \
[2110_HG001_subset]="$subset_dir/2110_GM12878_subset.txt" \
)

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

_info(){
  job=survivor_${sample}_${r}r${c}c
  name=jzh
  script=$root_dir/merge_sv.sh
}

run(){
    for sample in ${!sample_list[*]}
    do
    subset=${sample_list[$sample]}
    for r in `seq 1 6`
    do
    for c in `seq 1 15`
    do
      _info
      sbatch \
          -J ${name}_${job} \
          -A tangfuchou_g1 \
          -p cn-long \
          --qos tangfuchoucnl \
          -c 2 \
          -o $root_dir/log/stdout_$job.log \
          -e $root_dir/log/stderr_$job.log \
          $script $sample $subset $raw_dir $r $c $root_dir
    done
    done
    done
}

run
