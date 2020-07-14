#!/bin/bash

#This script regenerates test data from replicate RID Q-Y5F6

module load singularity/3.5.3
module load pigz/2.4

mkdir -p NEW_test_data

ln -sfn `readlink -e ./test_data/auth/credential.json` ~/.deriva/credential.json

mkdir -p ./NEW_test_data/bagit
singularity run 'docker://bicf/gudmaprbkfilexfer:1.3' deriva-download-cli dev.gudmap.org --catalog 2 ./workflow/conf/replicate_export_config.json . rid=Q-Y5F6
cp Replicate_Q-Y5F6.zip ./NEW_test_data/bagit/Replicate_Q-Y5F6.zip

mkdir -p ./NEW_test_data/fastq
unzip ./test_data/bagit/Replicate_Q-Y5F6.zip
singularity run 'docker://bicf/gudmaprbkfilexfer:1.3' bash ./workflow/scripts/bdbagFetch.sh Replicate_Q-Y5F6 Replicate_Q-Y5F6
cp Replicate_Q-Y5F6.R1.fastq.gz ./NEW_test_data/fastq/Replicate_Q-Y5F6.R1.fastq.gz
cp Replicate_Q-Y5F6.R2.fastq.gz ./NEW_test_data/fastq/Replicate_Q-Y5F6.R2.fastq.gz

mkdir -p ./NEW_test_data/fastq/small
singularity exec 'docker://bicf/seqtk:2.0.0' seqtk sample -s100 ./NEW_test_data/fastq/Replicate_Q-Y5F6.R1.fastq.gz 1000000 1> Q-Y5F6_1M.R1.fastq
singularity exec 'docker://bicf/seqtk:2.0.0' seqtk sample -s100 ./NEW_test_data/fastq/Replicate_Q-Y5F6.R2.fastq.gz 1000000 1> Q-Y5F6_1M.R2.fastq
pigz Q-Y5F6_1M.R1.fastq
pigz Q-Y5F6_1M.R2.fastq
cp Q-Y5F6_1M.R1.fastq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.R1.fastq.gz
cp Q-Y5F6_1M.R2.fastq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.R2.fastq.gz

mkdir -p ./NEW_test_data/meta
singularity run 'docker://bicf/trimgalore:1.1' trim_galore --gzip -q 25 --illumina --length 35 --basename Q-Y5F6_1M.se -j 20 ./NEW_test_data/fastq/small/Q-Y5F6_1M.R1.fastq.gz
singularity run 'docker://bicf/trimgalore:1.1' trim_galore --gzip -q 25 --illumina --length 35 --paired --basename Q-Y5F6_1M.pe -j 20 ./NEW_test_data/fastq/small/Q-Y5F6_1M.R1.fastq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.R2.fastq.gz
cp Q-Y5F6_1M.se_trimmed.fq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.se_trimmed.fq.gz
cp Q-Y5F6_1M.pe_R1_val_1.fq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.pe_R1_val_1.fq.gz
cp Q-Y5F6_1M.pe_R2_val_2.fq.gz ./NEW_test_data/fastq/small/Q-Y5F6_1M.pe_R2_val_2.fq.gz
cp Q-Y5F6_1M.R1.fastq.gz_trimming_report.txt ./NEW_test_data/meta/Q-Y5F6_1M.R1.fastq.gz_trimming_report.txt
cp Q-Y5F6_1M.R2.fastq.gz_trimming_report.txt ./NEW_test_data/meta/Q-Y5F6_1M.R2.fastq.gz_trimming_report.txt

touch metaTest.csv
echo 'Replicate_RID,Experiment_RID,Study_RID,Paired_End,File_Type,Has_Strand_Specific_Information,Used_Spike_Ins,Species' > metaTest.csv
echo 'Replicate_RID,Experiment_RID,Study_RID,uk,FastQ,no,no,Homo sapiens' >> metaTest.csv
cp metaTest.csv ./NEW_test_data/meta/metaTest.csv

mkdir -p ./NEW_test_data/bam
mkdir -p ./NEW_test_data/bam/small
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' hisat2 -p 20 --add-chrname --un-gz Q-Y5F6_1M.se.unal.gz -S Q-Y5F6_1M.se.sam -x /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/hisat2/genome --rna-strandness F -U ./NEW_test_data/fastq/small/Q-Y5F6_1M.se_trimmed.fq.gz --summary-file Q-Y5F6_1M.se.alignSummary.txt --new-summary
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools view -1 -@ 20 -F 4 -F 8 -F 256 -o Q-Y5F6_1M.se.bam Q-Y5F6_1M.se.sam
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools sort -@ 20 -O BAM -o Q-Y5F6_1M.se.sorted.bam Q-Y5F6_1M.se.bam
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools index -@ 20 -b Q-Y5F6_1M.se.sorted.bam Q-Y5F6_1M.se.sorted.bam.bai
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' hisat2 -p 20 --add-chrname --un-gz Q-Y5F6_1M.pe.unal.gz -S Q-Y5F6_1M.pe.sam -x /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/hisat2/genome --rna-strandness FR --no-mixed --no-discordant -1 ./NEW_test_data/fastq/small/Q-Y5F6_1M.pe_R1_val_1.fq.gz -2 ./test_data/fastq/small/Q-Y5F6_1M.pe_R2_val_2.fq.gz --summary-file Q-Y5F6_1M.pe.alignSummary.txt --new-summary
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools view -1 -@ 20 -F 4 -F 8 -F 256 -o Q-Y5F6_1M.pe.bam Q-Y5F6_1M.pe.sam
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools sort -@ 20 -O BAM -o Q-Y5F6_1M.pe.sorted.bam Q-Y5F6_1M.pe.bam
singularity run 'docker://bicf/gudmaprbkaligner:2.0.0' samtools index -@ 20 -b Q-Y5F6_1M.pe.sorted.bam Q-Y5F6_1M.pe.sorted.bam.bai
cp Q-Y5F6_1M.se.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.bam
cp Q-Y5F6_1M.pe.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.pe.bam
cp Q-Y5F6_1M.se.sorted.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.bam
cp Q-Y5F6_1M.se.sorted.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.bam.bai
cp Q-Y5F6_1M.pe.sorted.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.pe.sorted.bam
cp Q-Y5F6_1M.pe.sorted.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.pe.sorted.bam.bai
cp Q-Y5F6_1M.se.alignSummary.txt ./NEW_test_data/meta/Q-Y5F6_1M.se.alignSummary.txt
cp Q-Y5F6_1M.pe.alignSummary.txt ./NEW_test_data/meta/Q-Y5F6_1M.pe.alignSummary.txt

singularity run 'docker://bicf/gudmaprbkdedup:2.0.0' java -jar /picard/build/libs/picard.jar MarkDuplicates I=./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.bam O=Q-Y5F6_1M.se.deduped.bam M=Q-Y5F6_1M.se.deduped.Metrics.txt REMOVE_DUPLICATES=true
singularity run 'docker://bicf/gudmaprbkdedup:2.0.0' samtools sort -@ 20 -O BAM -o Q-Y5F6_1M.se.sorted.deduped.bam Q-Y5F6_1M.se.deduped.bam
singularity run 'docker://bicf/gudmaprbkdedup:2.0.0' samtools index -@ 20 -b Q-Y5F6_1M.se.sorted.deduped.bam Q-Y5F6_1M.se.sorted.deduped.bam.bai
cp Q-Y5F6_1M.se.sorted.deduped.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.bam
cp Q-Y5F6_1M.se.sorted.deduped.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.bam.bai
cp Q-Y5F6_1M.se.deduped.Metrics.txt /NEW_test_data/meta/Q-Y5F6_1M.se.deduped.Metrics.txt
cp Q-Y5F6_1M.se.deduped.Metrics.txt ./NEW_test_data/meta/Q-Y5F6_1M.se.deduped.Metrics.txt

for i in {"chr8","chr4","chrY"}; do 
      echo "samtools view -b ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.bam ${i} > Q-Y5F6_1M.se.sorted.deduped.${i}.bam; samtools index -@ 20 -b Q-Y5F6_1M.se.sorted.deduped.${i}.bam Q-Y5F6_1M.se.sorted.deduped.${i}.bam.bai;";
      done | singularity run 'docker://bicf/gudmaprbkdedup:2.0.0' parallel -j 20 -k
cp Q-Y5F6_1M.se.sorted.deduped.chr4.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chr4.bam
cp Q-Y5F6_1M.se.sorted.deduped.chr4.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chr4.bam.bai
cp Q-Y5F6_1M.se.sorted.deduped.chr8.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chr8.bam
cp Q-Y5F6_1M.se.sorted.deduped.chr8.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chr8.bam.bai
cp Q-Y5F6_1M.se.sorted.deduped.chrY.bam ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chrY.bam
cp Q-Y5F6_1M.se.sorted.deduped.chrY.bam.bai ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.chrY.bam.bai

mkdir -p ./NEW_test_data/counts
mkdir -p ./NEW_test_data/counts/small
singularity run 'docker://bicf/subread2:2.0.0' featureCounts -T 20 -a /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/genome.gtf -G /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/genome.fna -g 'gene_name' -o Q-Y5F6_1M.se.featureCounts -s 1 -R SAM --primary --ignoreDup ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.bam 
singularity run 'docker://bicf/subread2:2.0.0' Rscript ./workflow/scripts/calculateTPM.R --count Q-Y5F6_1M.se.featureCounts
cp Q-Y5F6_1M.se.featureCounts ./NEW_test_data/counts/small/Q-Y5F6_1M.se.featureCounts
cp Q-Y5F6_1M.se.countTable.csv ./NEW_test_data/counts/small/Q-Y5F6_1M.se.countTable.csv

mkdir -p ./NEW_test_data/bw
mkdir -p ./NEW_test_data/bw/small
singularity run 'docker://bicf/deeptools3.3:2.0.0' bamCoverage -p 20 -b ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.bam -o Q-Y5F6_1M.se.bw
cp Q-Y5F6_1M.se.bw ./NEW_test_data/bw/small/Q-Y5F6_1M.se.bw

mkdir -p ./NEW_test_data/fastqc
mkdir -p ./NEW_test_data/fastqc/small
singularity run 'docker://bicf/fastqc:2.0.0' ./NEW_test_data/fastq/small/Q-Y5F6_1M.R1.fastq.gz -o .
cp Q-Y5F6_1M.R1_fastqc.html ./NEW_test_data/fastqc/small/Q-Y5F6_1M.R1_fastqc.html
cp Q-Y5F6_1M.R1_fastqc.zip ./NEW_test_data/fastqc/small/Q-Y5F6_1M.R1_fastqc.zip

echo -e  "geneID\tchrom\ttx_start\ttx_end\tTIN" > Q-Y5F6_1M.se.sorted.deduped.tin.xls
for i in {"chr8","chr4","chrY"}; do
echo "tin.py -i ./NEW_test_data/bam/small/Q-Y5F6_1M.se.sorted.deduped.${i}.bam -r /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/bed/genome.bed; cat Q-Y5F6_1M.se.sorted.deduped.${i}.tin.xls | tr -s \"\\w\" \"\\t\" | grep -P \"\\t${i}\\t\";"; done | singularity run 'docker://bicf/rseqc3.0:2.0.0' parallel -j 20 -k >> Q-Y5F6_1M.se.sorted.deduped.tin.xls
cp Q-Y5F6_1M.se.sorted.deduped.tin.xls ./NEW_test_data/meta/Q-Y5F6_1M.se.sorted.deduped.tin.xls

