#! /bin/bash
# sort survivor VCF

f=$1

f_name=`echo $f | sed 's/.vcf//g'`
cat \
    <(cat $f | grep "##") \
    <(cat $f | grep -v "##" \
    | cut -f 1-8 \
    | awk -vOFS='\t' '{if(NR!=1) $3=NR-1; print}' \
    | sort -k1,1 -k2,2n) \
> $f_name.sorted.vcf

bgzip $f_name.sorted.vcf
tabix $f_name.sorted.vcf.gz
