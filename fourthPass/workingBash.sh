#!/usr/bin/bash 

for k in `ls /ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/*_L1_1.fq.gz`; do echo $k; echo $k | sed 's/L1_1/L1_2/g' ; done; 

 for k in `ls /ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/TR_* | grep "_L1_1.fq.gz" `; do echo "bash buildingSalmonscripts_mm10.sh $k"; done >> createdSalmon.sh



