#!/usr/bin/bash 




for k in `ls -1 ../../RNAseqData/H202SC19100044/Rawdata/`; do
	i=`ls  "../../RNAseqData/H202SC19100044/Rawdata/"$k `
	array=(`echo $i | sed 's/\s/\n/g'`)
	echo ${array[0]}
	echo
done 







