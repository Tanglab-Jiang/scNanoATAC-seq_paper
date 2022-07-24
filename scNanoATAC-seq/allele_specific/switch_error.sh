#! /bin/bash

chr=$1
root_dir=$2

snp_dir=$root_dir/snp_hapcut2/GM12878/
accuracy_dir=$root_dir/switch_error/
ref_1000g=~/ref/hg38/snp/snp_1000g_phased/

init_dir(){
    mkdir -p $1
    cd $1
}

init_dir $accuracy_dir

# extract hg001 trio SNP
ref_trio(){
    cd $ref_1000g
    zcat *chr${chr}.*phased.vcf.gz \
    | sed 's/chr//g' \
    | bcftools view \
        -O b \
        -s NA12891,NA12892 \
        - \
    > extract/chr${chr}.hg001_parents.bcf

    bcftools index extract/chr${chr}.hg001_parents.bcf
}

# create a ped fammily information file
create_ped(){
    echo -e \
    "1  HG001  NA12891  NA12892  2  0
    1  NA12891  0  0  1  0
    1  NA12892  0  0  2  0" \
    > GM12878.ped
}

# convert longshot vcf to bcf
snp2bcf(){
    tabix $snp_dir/GM12878.$chr.phased.VCF.gz

    bcftools view \
        -O b \
        $snp_dir/GM12878.$chr.phased.VCF.gz \
    > GM12878.$chr.bcf

    bcftools index GM12878.$chr.bcf
}

merge_query_trio(){
    bcftools merge \
        -O b \
        $ref_1000g/extract/chr${chr}.hg001_parents.bcf \
        GM12878.$chr.bcf \
    > GM12878.$chr.trio.bcf

    bcftools index GM12878.$chr.trio.bcf
}

get_switch_error(){
    bcftools +trio-switch-rate \
        GM12878.$chr.trio.bcf \
        -- -p GM12878.ped \
    > GM12878.chr$chr.report.txt
}

clean(){
    rm -f GM12878.$chr.trio.bcf GM12878.$chr.trio.bcf \
    GM12878.$chr.trio.bcf.csi GM12878.$chr.trio.bcf.csi
}

summary(){
    awk -vOFS='\t' '{if(NR==4)print "chr",$0}' GM12878.chr1.report.txt > summary.txt
    for chr in `seq 1 22`
    do
        awk -vOFS='\t' '{if(NR==5)print chr,$0}' chr=$chr GM12878.chr$chr.report.txt
    done \
    >> summary.txt
}

create_ped
snp2bcf
merge_query_trio
get_switch_error
