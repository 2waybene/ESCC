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

cores=6

# directory with genome reference FASTA and index files + name of the gene annotation file
#genome=/ddn/gs1/home/li11/refDB/mm10/STAR_index_mm10_RefSeq06Nov2019/
#gtf=/ddn/gs1/home/li11/refDB/mm10/mm10_ucsc_NCBIrefSeq_09Nov2019.gtf

##	trouble shoot, using Sara gtf and STAR_index

genome=/ddn/gs1/shared/dirib/reference_genomes/mm10/refseq20140406/STAR_index/STAR_index_refseq20140406/
gtf=/ddn/gs1/shared/fargod/reference_genomes/mm10/refseq20140406/mm10_RefSeq20140406.gtf
#	set up output filenames and locations

fastqc_out=/ddn/gs1/home/li11/project2019/RNAseqProj/results/fastqc/
align_out=/ddn/gs1/home/li11/project2019/RNAseqProj/results/STAR/${base}_
counts_input_bam=/ddn/gs1/home/li11/project2019/RNAseqProj/results/STAR/${base}_Aligned.sortedByCoord.out.bam
counts_output_bam=/ddn/gs1/home/li11/project2019/RNAseqProj/results/STAR/${base}_Aligned.sortedByName.out.bam
counts=/ddn/gs1/home/li11/project2019/RNAseqProj/results/counts/${base}_featurecounts.txt
fq1=/ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/${base}"_1.fq.gz"
fq2=/ddn/gs1/home/li11/project2019/RNAseqProj/raw_data/${base}"_2.fq.gz"
script=/ddn/gs1/home/li11/project2019/RNAseqProj/scripts/executeScripts/$base"_analysis.sh"

# Produce  STAR script

echo "STAR --runThreadN $cores --genomeDir $genome --readFilesIn <(zcat -c $fq1 ) <(zcat -c $fq2 ) --outFileNamePrefix $align_out --outFilterMultimapNmax 10 --outSAMstrandField intronMotif --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS"  > $script

echo "" >>  $script

# Create BAM index

#echo "samtools index $counts_input_bam"  >>  $script
echo "samtools sort  $counts_input_bam  -o  $counts_output_bam -O BAM "  >>  $script
echo ""  >>  $script

echo "samtools index $counts_output_bam"  >>  $script
echo ""  >>  $script

# Count mapped reads
#echo "featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_input_bam"  >> $script
echo "featureCounts -T $cores -s 2 -a $gtf -o $counts $counts_output_bam"  >> $script
echo "" >>  $script




