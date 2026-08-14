#!/bin/bash

# 处理第一个文件
awk 'BEGIN {FS=OFS="\t"} 
NR==1 {
    # 找到C_ID1列的位置
    for(i=1; i<=NF; i++) {
        if($i == "C_ID1") {
            start_col = i;
            break;
        }
    }
    print $0;
    next;
}
{
    # calculate sum from C_ID1
    sum = 0;
    for(i=start_col; i<=NF; i++) {
        sum += $i;
    }
    
    # get porportion
    for(i=start_col; i<=NF; i++) {
        if(sum != 0) {
            $i = $i / sum;
        } else {
            $i = 0;
        }
    }
    
    print $0;
}' /public/data/Topography_analysis/mutation_data/all.indel.ID83.partial.credit.txt \
> /public/data/Topography_analysis/mutation_data/all.indel.ID83.normalized.txt
