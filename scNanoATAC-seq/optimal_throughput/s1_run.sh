#!/bin/bash

root_dir=~/hyq_atac/tgs/analysis/optimal_throughput/
script=$root_dir/pipeline.sh

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

cell_type="eHAP1  GM12878  HEK293T  HFF1  K562"
throughput="1000 2000 3000 4000 5000 1500 2500 3500 4500"
seed=`seq 1 15`

info(){
    job=hyq_atac_${ct}_${tp}_${sd}
    name=jzh
}

run_by_cell_type(){
    threads=2
    
    for ct in $cell_type
    do
    for tp in $throughput
    do
    for sd in $seed
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
        $script $ct $tp $sd $root_dir
    done
    done
    done
}

run_by_throughput(){
    threads=2
    
    for tp in $throughput
    do
    for sd in $seed
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
        $script . $tp $sd $root_dir
    done
    done
}

run_by_cell_type
run_by_throughput
