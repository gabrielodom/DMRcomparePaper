#!/bin/bash
 
cd "Comb-P-Folders"
 
SEED_FOLDERS=*
 
for seed in $SEED_FOLDERS
do
    echo "Processing $seed seed folder"
    cd $seed
     
    #str="seed-001"
    prefix1="seed-"
    seed_num=$(echo "$seed" | sed -e "s/^$prefix1//")
    echo $seed_num
     
    SEED_DIST_FOLDERS=*
 
    for dist in $SEED_DIST_FOLDERS
    do
        echo "Processing $dist seed dist folder"
        cd $dist
        
        prefix2="dist-"
        dist_num=$(echo "$dist" | sed -e "s/^$prefix2//")
        echo $dist_num
         
        FOLDERS=*
 
        for f in $FOLDERS
        do
            echo "Processing $f file folder"
            cd $f
            
            run_start=`date +%s%N`
            comb-p pipeline -c 4 --seed 1e-1 --dist $dist_num -p Train --region-filter-p $seed_num --annotate hg19 *.bed  
            run_end=`date +%s%N`
            run_runtime=$((run_end-run_start))
            run_runtime_milli=$((run_runtime/1000000))
            echo $run_runtime_milli >> run_time.txt 2>&1 # in millisecond
            #echo $seed $seed_dist $f
             
            cd ..
        done
        
        cd ..
    done
    
    cd ..    
done
