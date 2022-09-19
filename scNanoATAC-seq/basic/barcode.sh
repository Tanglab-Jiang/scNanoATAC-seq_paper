#! /bin/bash

barcode_ref=~/ref/barcode/nanopore_dual/96_barcode.txt

index=$1
for i in $index
do
    bc=`cat $barcode_ref \
    | sed -n ${i}p \
    | cut -f 2`
    echo -e ">$i\n${bc}"
done
