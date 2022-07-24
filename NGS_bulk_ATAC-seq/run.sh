#! /bin/bash

root_dir=~/hyq_atac/ngs/
script=$root_dir/pipeline.sh

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

lib_list=""

info(){
    job=hyq_atac_${lib}_${sample}
    name=jzh
}

run_one(){
    threads=6
    
    info
    for lib in $lib_list
    do
    sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p cn-long \
        --qos tangfuchoucnl \
        -c $threads \
        -o $root_dir/log/stdout_$job.log \
        -e $root_dir/log/stderr_$job.log \
        $script $lib . $threads $root_dir
    done
}

run_sample(){
    threads=2

    for lib in $lib_list
    do
    sample_list=`ls $root_dir/raw.data/$lib/`
    for sample in $sample_list
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
            $script $lib $sample $threads $root_dir
    done
    done
}


run_one
run_sample