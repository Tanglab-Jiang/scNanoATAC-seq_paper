#! /bin/bash

root_dir=~/hyq_atac/tgs/analysis/haplotype/
script=$root_dir/switch_error.sh

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

subset_list="GM12878"

info(){
    job=hyq_${subset}_${chr}
    name=jzh
}

run_long_chr(){
    threads=2

    for subset in $subset_list
    do
    for chr in `seq 1 22` X
    do
    info
    sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p cn-long \
        --qos tangfuchoucnl \
        -c $threads \
        -o $root_dir/log/stdout_$job.log \
        -e $root_dir/log/stderr_$job.log \
        $script $chr $root_dir
    done
    done
}

run_long_chr

