#! /bin/bash

root_dir=~/hyq_atac/tgs/
script=$root_dir/analysis/cross_contamination/count_reads.sh

if [ ! -f $root_dir/log/ ]
then
    mkdir -p $root_dir/log/
fi

lib_list="2107_ALC1  2107_ALC2"
group=`seq 73 92`

info(){
    job=hyq_atac_$lib
    name=jzh
}

run_long_one(){
    threads=2
    
    for lib in $lib_list
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
            $script $lib . $threads $root_dir
    done
}

run_long_ht(){
    threads=1

    for lib in $lib_list
    do
        for g in $group
        do
        for id in `seq 1 48`
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
                $script $lib ${g}_${id} $threads $root_dir
        done
        done
    done
}

run_long_ht
run_long_one
