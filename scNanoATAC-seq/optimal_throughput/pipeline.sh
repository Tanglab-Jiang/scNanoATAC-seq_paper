#! /bin/bash

cell_type=$1
throughput=$2
seed=$3
root_dir=$4

script=$root_dir/simulate_cell.sh
arrow_script=$root_dir/create_arrow.R
clustering_script=$root_dir/clustering.R

raw_fragment=$(dirname $root_dir)/merge_fragment/data/$cell_type.flanked.bed.gz
downsample_dir=$root_dir/downsample/$throughput/
archr_dir=$root_dir/archr/${throughput}.${seed}/

init_dir(){
  mkdir -p $1/log/
  cp $script $1/log/
  cd $1
}

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

# run by cell type
downsample(){
    init_dir $downsample_dir

    zcat $raw_fragment \
    | shuf --random-source=<(get_seeded_random $seed) \
    | awk -vOFS='\t' '{
        if(NR<=20*10**9/tp) {
            id+=1
            if (id > 500) id=1
            if ($1 !~ /chr/) $1="chr"$1
            cell=ct"_"id
            print $1,$2,$3,cell,1
        }
    }' ct=$cell_type tp=$throughput \
    | sort -k 1,1 -k 2,2n \
    > $cell_type.downsampled.seed_$seed.fragments.sorted.bed

    bgzip $cell_type.downsampled.seed_$seed.fragments.sorted.bed
    tabix $cell_type.downsampled.seed_$seed.fragments.sorted.bed.gz
}

# run by cell type
# activate r40 first
create_arrow_file(){
  init_dir $archr_dir
  
  if [ ! -s $throughput.$cell_type.$seed.arrow ]
  then
    source activate r40  
    Rscript $arrow_script \
      $cell_type $throughput $seed $(dirname $downsample_dir)
  fi
}

# run by throughput
clustering(){
  cd $archr_dir

  if [ ! -s $throughput.$seed.clustering.Rds ]
  then
    source activate r40  
    Rscript $clustering_script \
      $throughput $seed $archr_dir
  fi
}

downsample
create_arrow_file

clustering
