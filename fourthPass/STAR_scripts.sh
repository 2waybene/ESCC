#!/usr/bin/bash

##	Don't have the zipped files, and need to re-write
#	STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_1/K70_1_CRRA190004025-1a_HMHVVDSXX_L1_1.fq) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_1/K70_1_CRRA190004025-1a_HMHVVDSXX_L1_2.fq) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70_1_w_GTF_ &

##	use the following instead

STAR  --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn  ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_1/K70_1_CRRA190004025-1a_HMHVVDSXX_L1_1.fq ~/project2019/
ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_1/K70_1_CRRA190004025-1a_HMHVVDSXX_L1_2.fq  --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf   --outFileNamePrefix starAln/K70_1_w_GTF_  &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_2/K70_2_CRRA190004026-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_2/K70_2_CRRA190004026-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70_2_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_3/K70_3_CRRA190004027-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70_3/K70_3_CRRA190004027-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70_3_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_1/K70MD_1_CRRA190004033-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_1/K70MD_1_CRRA190004033-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70MD_1_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_2/K70MD_2_CRRA190004034-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_2/K70MD_2_CRRA190004034-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70MD_2_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_3/K70MD_3_CRRA19H000008-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70MD_3/K70MD_3_CRRA19H000008-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70MD_3_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_1/K70P10_1_CRRA190004029-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_1/K70P10_1_CRRA190004029-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70P10_1_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_2/K70P10_2_CRRA190004030-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_2/K70P10_2_CRRA190004030-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70P10_2_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_3/K70P10_3_CRRA190004031-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/K70P10_3/K70P10_3_CRRA190004031-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/K70P10_3_w_GTF_ &


STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_1/NKO70_1_CRRA190004021-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_1/NKO70_1_CRRA190004021-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/NKO70_1_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_3/NKO70_3_CRRA190004023-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_3/NKO70_3_CRRA190004023-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/NKO70_3_w_GTF_ &

STAR --genomeDir /ddn/gs1/home/li11/refDB/hg38/STAR_index_hg38noAlt_RefSeq10oct2019/ --readFilesIn <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_4/NKO70_4_CRRA190004024-1a_HMHVVDSXX_L1_1.fq.gz) <(zcat -c ~/project2019/ThuyAi_track_stuff/bamDir/RNAseqData/H202SC19100044/Rawdata/NKO70_4/NKO70_4_CRRA190004024-1a_HMHVVDSXX_L1_2.fq.gz) --runThreadN 6 --sjdbGTFfile /ddn/gs1/home/li11/refDB/hg38/hg38_RefSeq10oct2019.gtf --outFileNamePrefix ~/project2019/ThuyAi_track_stuff/bamDir/analysisDir/starAln/NKO70_4_w_GTF_ &








