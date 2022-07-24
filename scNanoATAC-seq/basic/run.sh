#!/bin/bash

root=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/jiangzh/
root_dir=$root/hyq_atac/tgs/
script=$root_dir/pipeline.sh

if [ ! -f $root_dir/log ]
then
    mkdir -p $root_dir/log
fi

# 96-cell library
#lib_list="SMA-HG01  SMA-K562  20210222_ATAC_ONT 2104_293T  2104_HePG2  2104_HFF1 \
#	2105_H293T-96  2105_HG01-96  2105_K562-96  2105_HepG2-96  2105_HFF1-96" 

#lib_list="2106_X-HG001_12"
#barcode="81 82"

#lib_list="2106_X-HG001_34"
#barcode="83 84"

#lib_list="2105_AL480"
#barcode=`seq 81 90`

#lib_list="2107_AL960_1  2107_AL960_2  2107_ALC1  2107_ALC2  2108_EHF  2110_GK  2110_HH  2110_CD19  2110_CD48  2110_PBMC"
#barcode=`seq 73 92`

# 220111
# quality control redo
lib_list="2107_ALC2"
barcode=`seq 73 92`

info(){
    job=hyq_atac_$lib
    name=jzh
}

run_long_one(){
    threads=10
    
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
        $script $lib . $threads $root_dir "$barcode"
    done
}

run_long_96(){
    threads=4

    for lib in $lib_list
    do
        for bc in `seq 1 96`
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
            $script $lib sc$bc $threads $root_dir
        done
    done
}

run_long_ht(){
    threads=4

    for lib in $lib_list
    do
        for bc1 in $barcode
        do
        for bc2 in `seq 1 48`
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
            $script $lib ${bc1}_${bc2} $threads $root_dir
        done
        done
    done
}

run_fat_one(){
    threads=5

    for lib in $lib_list
    do
        info
        sbatch \
        -J ${name}_${job} \
        -A tangfuchou_g1 \
        -p fat8way \
        --qos tangfuchouf8w \
        -c $threads \
        -o $root_dir/log/out_$job.log \
        -e $root_dir/log/err_$job.log \
        $script $lib . $threads $root_dir
    done
}

run_one_directly(){
  for lib in $lib_list
  do
    bash $script $lib . 4 $root_dir &
  done
}

run_96_directly(){
    for lib in $lib_list
    do
      for cell in `seq 1 96`
      do
        bash $script $lib sc$cell 1 $root_dir & 
      done
    done
}

run_ht_directly(){
    threads=1

    for lib in $lib_list
    do
        for bc1 in $barcode
        do
        for bc2 in `seq 1 48`
        do
            bash $script $lib ${bc1}_${bc2} $threads $root_dir 
        done
        done
    done
}

#run_long_one
run_long_ht
#run_fat_one
#run_one_directly
#run_ht_directly
