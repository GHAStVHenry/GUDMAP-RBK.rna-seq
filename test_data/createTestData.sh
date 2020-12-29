#!/bin/bash

#This script regenerates test data from replicate RID 16-DNBG

module load singularity/3.5.3
module load pigz/2.4

ln -sfn /project/BICF/BICF_Core/shared/gudmap/test_data/* ../test_data/

mkdir -p NEW_test_data

ln -sfn ./test_data/auth/credential.json ~/.deriva/credential.json

mkdir -p ./NEW_test_data/bag
singularity run 'docker://gudmaprbk/deriva1.3:1.0.0' deriva-download-cli staging.gudmap.org --catalog 2 ../workflow/conf/Replicate_For_Input_Bag.json . rid=16-DNBG
cp 16-DNBG_inputBag.zip ./NEW_test_data/bag/16-DNBG_inputBag.zip

mkdir -p ./NEW_test_data/fastq
unzip ./NEW_test_data/bag/16-DNBG_inputBag.zip
singularity run 'docker://gudmaprbk/deriva1.3:1.0.0' bash ../workflow/scripts/bdbagFetch.sh 16-DNBG_inputBag 16-DNBG
cp 16-DNBG.R1.fastq.gz ./NEW_test_data/fastq/16-DNBG.R1.fastq.gz
cp 16-DNBG.R2.fastq.gz ./NEW_test_data/fastq/16-DNBG.R2.fastq.gz

mkdir -p ./NEW_test_data/fastq/small
singularity exec 'docker://gudmaprbk/seqtk1.3:1.0.0' seqtk sample -s100 ./NEW_test_data/fastq/16-DNBG.R1.fastq.gz 1000000 1> 16-DNBG_1M.R1.fastq
singularity exec 'docker://gudmaprbk/seqtk1.3:1.0.0' seqtk sample -s100 ./NEW_test_data/fastq/16-DNBG.R2.fastq.gz 1000000 1> 16-DNBG_1M.R2.fastq
pigz 16-DNBG_1M.R1.fastq
pigz 16-DNBG_1M.R2.fastq
cp 16-DNBG_1M.R1.fastq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.R1.fastq.gz
cp 16-DNBG_1M.R2.fastq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.R2.fastq.gz

mkdir -p ./NEW_test_data/meta
singularity run 'docker://gudmaprbk/trimgalore0.6.5:1.0.0' trim_galore --gzip -q 25 --illumina --length 35 --basename 16-DNBG_1M.se -j 20 ./NEW_test_data/fastq/small/16-DNBG_1M.R1.fastq.gz
singularity run 'docker://gudmaprbk/trimgalore0.6.5:1.0.0' trim_galore --gzip -q 25 --illumina --length 35 --paired --basename 16-DNBG_1M.pe -j 20 ./NEW_test_data/fastq/small/16-DNBG_1M.R1.fastq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.R2.fastq.gz
cp 16-DNBG_1M.se_trimmed.fq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.se_trimmed.fq.gz
cp 16-DNBG_1M.pe_val_1.fq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.pe_val_1.fq.gz
cp 16-DNBG_1M.pe_val_2.fq.gz ./NEW_test_data/fastq/small/16-DNBG_1M.pe_val_2.fq.gz
cp 16-DNBG_1M.R1.fastq.gz_trimming_report.txt ./NEW_test_data/meta/16-DNBG_1M.R1.fastq.gz_trimming_report.txt
cp 16-DNBG_1M.R2.fastq.gz_trimming_report.txt ./NEW_test_data/meta/16-DNBG_1M.R2.fastq.gz_trimming_report.txt

touch metaTest.csv
echo 'Replicate_RID,Experiment_RID,Study_RID,Paired_End,File_Type,Has_Strand_Specific_Information,Used_Spike_Ins,Species' > metaTest.csv
echo 'Replicate_RID,Experiment_RID,Study_RID,uk,FastQ,no,no,Homo sapiens' >> metaTest.csv
cp metaTest.csv ./NEW_test_data/meta/metaTest.csv

mkdir -p ./NEW_test_data/bam
mkdir -p ./NEW_test_data/bam/small
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' hisat2 -p 20 --add-chrname --un-gz 16-DNBG_1M.se.unal.gz -S 16-DNBG_1M.se.sam -x /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/hisat2/genome --rna-strandness F -U ./NEW_test_data/fastq/small/16-DNBG_1M.se_trimmed.fq.gz --summary-file 16-DNBG_1M.se.alignSummary.txt --new-summary
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools view -1 -@ 20 -F 4 -F 8 -F 256 -o 16-DNBG_1M.se.bam 16-DNBG_1M.se.sam
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools sort -@ 20 -O BAM -o 16-DNBG_1M.se.sorted.bam 16-DNBG_1M.se.bam
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools index -@ 20 -b 16-DNBG_1M.se.sorted.bam 16-DNBG_1M.se.sorted.bam.bai
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' hisat2 -p 20 --add-chrname --un-gz 16-DNBG_1M.pe.unal.gz -S 16-DNBG_1M.pe.sam -x /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/hisat2/genome --rna-strandness FR --no-mixed --no-discordant -1 ./NEW_test_data/fastq/small/16-DNBG_1M.pe_val_1.fq.gz -2 ./NEW_test_data/fastq/small/16-DNBG_1M.pe_val_2.fq.gz --summary-file 16-DNBG_1M.pe.alignSummary.txt --new-summary
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools view -1 -@ 20 -F 4 -F 8 -F 256 -o 16-DNBG_1M.pe.bam 16-DNBG_1M.pe.sam
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools sort -@ 20 -O BAM -o 16-DNBG_1M.pe.sorted.bam 16-DNBG_1M.pe.bam
singularity run 'docker://gudmaprbk/hisat2.2.1:1.0.0' samtools index -@ 20 -b 16-DNBG_1M.pe.sorted.bam 16-DNBG_1M.pe.sorted.bam.bai
cp 16-DNBG_1M.se.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.bam
cp 16-DNBG_1M.pe.bam ./NEW_test_data/bam/small/16-DNBG_1M.pe.bam
cp 16-DNBG_1M.se.sorted.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.bam
cp 16-DNBG_1M.se.sorted.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.bam.bai
cp 16-DNBG_1M.pe.sorted.bam ./NEW_test_data/bam/small/16-DNBG_1M.pe.sorted.bam
cp 16-DNBG_1M.pe.sorted.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.pe.sorted.bam.bai
cp 16-DNBG_1M.se.alignSummary.txt ./NEW_test_data/meta/16-DNBG_1M.se.alignSummary.txt
cp 16-DNBG_1M.pe.alignSummary.txt ./NEW_test_data/meta/16-DNBG_1M.pe.alignSummary.txt

singularity run 'docker://gudmaprbk/picard2.23.9:1.0.0' java -jar /picard/build/libs/picard.jar MarkDuplicates I=./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.bam O=16-DNBG_1M.se.deduped.bam M=16-DNBG_1M.se.deduped.Metrics.txt REMOVE_DUPLICATES=true
singularity run 'docker://gudmaprbk/picard2.23.9:1.0.0' samtools sort -@ 20 -O BAM -o 16-DNBG_1M.se.sorted.deduped.bam 16-DNBG_1M.se.deduped.bam
singularity run 'docker://gudmaprbk/picard2.23.9:1.0.0' samtools index -@ 20 -b 16-DNBG_1M.se.sorted.deduped.bam 16-DNBG_1M.se.sorted.deduped.bam.bai
cp 16-DNBG_1M.se.deduped.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.deduped.bam
cp 16-DNBG_1M.se.sorted.deduped.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.bam
cp 16-DNBG_1M.se.sorted.deduped.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.bam.bai
cp 16-DNBG_1M.se.deduped.Metrics.txt ./NEW_test_data/meta/16-DNBG_1M.se.deduped.Metrics.txt

for i in {"chr8","chr4","chrY"}; do 
      echo "samtools view -b ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.bam ${i} > 16-DNBG_1M.se.sorted.deduped.${i}.bam; samtools index -@ 20 -b 16-DNBG_1M.se.sorted.deduped.${i}.bam 16-DNBG_1M.se.sorted.deduped.${i}.bam.bai;";
      done | singularity run 'docker://gudmaprbk/picard2.23.9:1.0.0' parallel -j 20 -k
cp 16-DNBG_1M.se.sorted.deduped.chr4.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chr4.bam
cp 16-DNBG_1M.se.sorted.deduped.chr4.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chr4.bam.bai
cp 16-DNBG_1M.se.sorted.deduped.chr8.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chr8.bam
cp 16-DNBG_1M.se.sorted.deduped.chr8.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chr8.bam.bai
cp 16-DNBG_1M.se.sorted.deduped.chrY.bam ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chrY.bam
cp 16-DNBG_1M.se.sorted.deduped.chrY.bam.bai ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.chrY.bam.bai

mkdir -p ./NEW_test_data/counts
mkdir -p ./NEW_test_data/counts/small
ln -s /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/geneID.tsv
ln -s /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/Entrez.tsv
singularity run 'docker://gudmaprbk/subread2.0.1:1.0.0' featureCounts -T 20 -a /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/genome.gtf -G /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/genome.fna -g 'gene_name' --extraAttributes 'gene_id' -o 16-DNBG_1M.se_countData -s 1 -R SAM --primary --ignoreDup ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.bam 
singularity run 'docker://gudmaprbk/subread2.0.1:1.0.0' Rscript ../workflow/scripts/calculateTPM.R --count 16-DNBG_1M.se_countData
singularity run 'docker://gudmaprbk/subread2.0.1:1.0.0' Rscript ../workflow/scripts/convertGeneSymbols.R --repRID 16-DNBG_1M.se
cp 16-DNBG_1M.se_countData ./NEW_test_data/counts/small/16-DNBG_1M.se_countData
cp 16-DNBG_1M.se.countTable.csv ./NEW_test_data/counts/small/16-DNBG_1M.se.countTable.csv
cp 16-DNBG_1M.se_tpmTable.csv ./NEW_test_data/counts/small/16-DNBG_1M.se_tpmTable.csv

mkdir -p ./NEW_test_data/bw
mkdir -p ./NEW_test_data/bw/small
singularity run 'docker://gudmaprbk/deeptools3.5.0:1.0.0' bamCoverage -p 20 -b ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.bam -o 16-DNBG_1M.se.bw
cp 16-DNBG_1M.se.bw ./NEW_test_data/bw/small/16-DNBG_1M.se.bw

mkdir -p ./NEW_test_data/fastqc
mkdir -p ./NEW_test_data/fastqc/small
singularity run 'docker://gudmaprbk/fastqc0.11.9:1.0.0' fastqc ./NEW_test_data/fastq/small/16-DNBG_1M.R1.fastq.gz -o .
cp 16-DNBG_1M.R1_fastqc.html ./NEW_test_data/fastqc/small/16-DNBG_1M.R1_fastqc.html
cp 16-DNBG_1M.R1_fastqc.zip ./NEW_test_data/fastqc/small/16-DNBG_1M.R1_fastqc.zip

echo -e  "geneID\tchrom\ttx_start\ttx_end\tTIN" > 16-DNBG_1M.se.sorted.deduped.tin.xls
for i in {"chr8","chr4","chrY"}; do
echo "tin.py -i ./NEW_test_data/bam/small/16-DNBG_1M.se.sorted.deduped.${i}.bam -r /project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12.v31/bed/genome.bed; cat 16-DNBG_1M.se.sorted.deduped.${i}.tin.xls | tr -s \"\\w\" \"\\t\" | grep -P \"\\t${i}\\t\";"; done | singularity run 'docker://gudmaprbk/rseqc4.0.0:1.0.0' parallel -j 20 -k >> 16-DNBG_1M.se.sorted.deduped.tin.xls
cp 16-DNBG_1M.se.sorted.deduped.tin.xls ./NEW_test_data/meta/16-DNBG_1M.se.sorted.deduped.tin.xls

chgrp -R BICF_Core ./NEW_test_data
chmod -R 750 ./NEW_test_data
