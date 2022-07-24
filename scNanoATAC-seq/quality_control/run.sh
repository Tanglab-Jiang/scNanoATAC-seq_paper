#!/bin/sh

root_dir=~/hyq_atac/tgs/
script=$root_dir/quality.sh

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

lib_list="SMA-HG01  SMA-K562  20210222_ATAC_ONT 2104_293T  2104_HePG2  2104_HFF1 \
  2105_H293T-96  2105_HG01-96  2105_K562-96  2105_HepG2-96  2105_HFF1-96  2105_AL480 \
  2106_X-HG001_12  2106_X-HG001_34  2107_AL960_1  2107_AL960_2  2107_ALC1  2107_ALC2 \
  2107_scTAG_X-HG001  2108_BT  2108_EHF  2110_NH_novogene"

info(){
  job=hyq_atac_$lib
  name=jzh
}

run_lib(){
    threads=4
    
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
        $script $lib $threads $root_dir
    done
}

run_one(){
    threads=4
    
    info
    sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p cn-long \
        --qos tangfuchoucnl \
        -c $threads \
        -o $root_dir/log/stdout_$job.log \
        -e $root_dir/log/stderr_$job.log \
        $script . $threads $root_dir
}

run_lib
run_one
