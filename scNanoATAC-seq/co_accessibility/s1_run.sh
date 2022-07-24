#! /bin/bash

root_dir=~/hyq_atac/tgs/analysis/co_accessibility
script=$root_dir/co_accessibility_ks_test.R

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

info(){
    job=hyq_atac_${ct}_${chr}
    name=jzh
}

run_long_one(){
    threads=20
    
    cell_type="B CD4_T CD8_T eHAP1 K562 GM12878 HEK293T HFF1 Monocyte"

    for ct in $cell_type
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
        <(echo -e "#! /bin/bash \n Rscript $script $ct $chr")
    done
    done
}

run_long_one
