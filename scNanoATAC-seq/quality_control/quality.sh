#!/bin/bash

library=$1
threads=$2
root_dir=$3

raw_dir=$root_dir/raw.data/$library/
trim_dir=$root_dir/trim/$library/
align_dir=$root_dir/alignment/$library/
quality_dir=$root_dir/quality/$library/

#align_dir=$root_dir/alignment_mm10/$library/
#quality_dir=$root_dir/quality_mm10/$library/

script=$root_dir/quality.sh

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

raw_qc(){
    init_dir $quality_dir/fastq/
    Rscript \
        $(which MinIONQC.R) \
        -i $raw_dir/*summary.txt.gz \
        -p $threads \
        -o . 
    multiqc -f .
}

trim_qc(){
    init_dir $quality_dir/trim/
    multiqc -f $trim_dir
}

align_qc(){
    init_dir $quality_dir/align/

    echo -e cell"\t"count_q30"\t"count_dedup"\t"count_MT \
        > align.txt

    cell_list=`\
    cat $quality_dir/trim/multiqc_data/multiqc_cutadapt.txt \
    | awk '{ if(NR!=1) print $1 }'`

    for cell in $cell_list
    do
        q30_bam=$align_dir/$library.$cell.mapQ30.sorted.bam
        dedup_bam=$align_dir/$library.$cell.mapQ30.rmdup.sorted.bam

        if [[ -s $q30_bam ]]
        then
            n_q30=`samtools view -F 2048,2064 -c -@ $threads $q30_bam`
        else
            n_q30=0
        fi

        if [[ -s $dedup_bam ]]
        then
            n_dedup=`samtools view -F 2048,2064 -c -@ $threads $dedup_bam`
            n_MT=`samtools view -c -@ $threads $dedup_bam MT`
        else
            n_dedup=0
            n_MT=0
        fi

        echo -e $cell"\t"$n_q30"\t"$n_dedup"\t"$n_MT \
            >> align.txt
    done
}

coverage_qc(){
    init_dir $quality_dir/coverage/genomecov/

    echo -e cell"\t"coverage \
        > ../$library.coverage.txt

    cell_list=`\
    cat $quality_dir/trim/multiqc_data/multiqc_cutadapt.txt \
    | awk '{ if(NR!=1) print $1 }'`

    for cell in $cell_list
    do
        dedup_bam=$align_dir/$library.$cell.mapQ30.rmdup.sorted.bam

        if [[ -s $dedup_bam ]]
        then
            if [[ ! -s $library.$cell.coverage.bedGraph ]]
            then
                bedtools genomecov \
                    -ibam $dedup_bam \
                    -bga \
                | sed 's/^chr//g' \
                | grep -vE '^KI|^GL|^MT' \
                > $library.$cell.coverage.bedGraph
            fi

            cat $library.$cell.coverage.bedGraph \
            | awk -vOFS='\t' \
                '{
                    len=$3-$2; s+=len
                    if($4>0) c+=len
                } END {
                    print cell, c/s
                }' \
            library=$library cell=$cell \
            >> ../$library.coverage.txt
        else
	  echo -e "$cell\t0" >> ../$library.coverage.txt
	fi
    done
}

# run by library
collect_sc(){
    cd $(dirname $quality_dir)

    echo -e \
    "cell\traw_reads\treads_after_trimming\tbases_after_trimming\tmean_insert_size\ttrimming_rate\tQ30_reads\tQ30_mapping_rate\tunique_reads\tunique_rate\tMT_reads\tcoverage" \
    > $library.single-cell.quality.txt

    trim_report=$quality_dir/trim/multiqc_data/multiqc_cutadapt.txt
    align_report=$quality_dir/align/align.txt
    coverage_report=$quality_dir/coverage/${library}.coverage.txt

    paste $trim_report $align_report $coverage_report \
    | awk -vOFS='\t' \
    '{
        if(NR!=1){
            if ($5==0){
                print $1,$2,$5,$7,0,$8/100,$10,0,$11,0,$12,$14
            }else if ($10==0){
                print $1,$2,$5,$7,$7/$5,$8/100,$10,$10/$5,$11,0,$12,$14
            }else{
                print $1,$2,$5,$7,$7/$5,$8/100,$10,$10/$5,$11,$11/$10,$12,$14
            }
        }
    }' \
    >> $library.single-cell.quality.txt
}

clean_bg(){
  rm -f $quality_dir/coverage/genomecov/*.coverage.bedGraph
}

# run one
collect_lib(){
    cd $(dirname $quality_dir)
    
    echo -e \
    "library\ttotal_gigabases\ttotal_reads\tN50_length\tmean_length\tmedian_length\tmedian_Q\tdemultiplex_rate" \
    > library.quality.txt

    for lib in `ls $(dirname $raw_dir)`
    do
        fq_report=$(dirname $quality_dir)/$lib/fastq/multiqc_data/multiqc_minionqc.txt
        trim_report=$(dirname $quality_dir)/$lib/trim/multiqc_data/multiqc_cutadapt.txt
        
        trimmed_reads=`awk '{if(NR!=1){sum+=$2}}END{print sum}' $trim_report`
        awk -vOFS='\t' \
            '{ if (NR!=1) print $1,$2,$3,$4,$5,$6,$9,tr/$3 }' \
            tr=$trimmed_reads \
            $fq_report
    done \
    >> library.quality.txt
}

raw_qc
trim_qc
align_qc
coverage_qc
collect_sc

collect_lib

clean_bg

