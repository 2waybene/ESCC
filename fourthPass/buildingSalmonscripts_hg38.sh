#!/bin/bash/

# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis_on_allfiles.sh <name of fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file

fq=$1

# grab base of filename for naming outputs

file=`basename $fq`
if [[ $file =~ ^(.*)(_1\.fq.gz)$ ]]; 
then 
  base=${BASH_REMATCH[1]} ; 
else 
  echo "Not proper format"; 
fi


# specify the number of cores to use


##	trouble shoot, using Sara gtf and STAR_index

#	set up output filenames and locations
#mm10_salmon_index=/ddn/gs1/home/li11/refDB/mm10/salmonIndex/mm10_salmon_index
hg38_salmon_index=/ddn/gs1/home/li11/refDB/hg38/salmonIndex/hg38_salmon_index
salmon_out=/ddn/gs1/home/li11/project2019/RNAseqProj/results/Salmon/${base}".salmon"
fq1=/ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/${base}"_1.fq.gz"
fq2=/ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/${base}"_2.fq.gz"
script=/ddn/gs1/home/li11/project2019/RNAseqProj/scripts/executeScripts/$base"_analysis.sh"

# Produce  STAR script

echo "salmon quant -i $hg38_salmon_index  -l MU -1 <(zcat -c $fq1 ) -2 <(zcat -c $fq2 ) -o $salmon_out "  > $script

echo "" >>  $script





