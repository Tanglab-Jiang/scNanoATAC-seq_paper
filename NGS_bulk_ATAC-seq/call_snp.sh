#! /bin/bash

ref_genome_fa=~/ref/hg38/GRCh38.fa
hg001=~/ref/hg38/giab_benchmark/hg001_na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz

get_ref(){
  zcat $hg001 \
  | sed 's/chr//g' \
  | grep -E "#|0\|1|1\|0" \
  | bgzip -c \
  > ref.hg001.het.vcf.gz

  tabix ref.hg001.het.vcf.gz
}

# fileds in VCF to be analyzed: AB(alternate allele / total count), DP(depth)
call_snp(){
  mkdir freebayes
  for chr in `seq 1 22` X
  do
  {
    freebayes \
        -f $ref_genome_fa \
        -@ ref.hg001.het.vcf.gz \
        -l \
        -r $chr \
        --min-alternate-fraction 0 \
        --min-alternate-count 0 \
        --min-alternate-total 0 \
        ~/hyq_atac/ngs/alignment/210819_GM12878/NGS-GE-GM12878-X1.mapQ30.sorted.dedup.bam \
    | bgzip -c \
    > ./freebayes/NGS-GE-GM12878-X1.giab_hg001_het.chr$chr.vcf.gz

    tabix ./freebayes/NGS-GE-GM12878-X1.giab_hg001_het.chr$chr.vcf.gz
  } & 
  done
}

concat_snp(){
  bcftools concat \
    ./freebayes/NGS-GE-GM12878-X1.giab_hg001_het.chr*.vcf.gz \
    -O z \
    -o NGS-GE-GM12878-X1.giab_hg001_het.vcf.gz
  
  tabix NGS-GE-GM12878-X1.giab_hg001_het.vcf.gz
}

get_ref
call_snp
concat_snp
