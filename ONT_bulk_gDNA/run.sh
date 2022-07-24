#! /bin/bash

root_dir=~/hyq_atac/tgs_gDNA/
script=$root_dir/pipeline.sh

if [ ! -f $root_dir/log/ ]
then
    mkdir -p $root_dir/log/
fi

_info(){
    job=hyq_gDNA_$sample
    name=jzh
}

sample_list="2201_K562_gDNA_HYQ"

run_long_single(){
    threads=15
    
    _info
    sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p cn-long \
        --qos tangfuchoucnl \
        -c $threads \
        -o $root_dir/log/out_$job.log \
        -e $root_dir/log/err_$job.log \
        $script $sample $threads $root_dir
}

run_fat_single(){
    threads=2

    _info
    sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p fat4way \
        --qos tangfuchouf4w \
        -c $threads \
        -o $root_dir/log/out_$job.log \
        -e $root_dir/log/err_$job.log \
        $script $sample $threads $root_dir
}

for sample in $sample_list
do
   run_long_single 
   #run_fat_single
done
