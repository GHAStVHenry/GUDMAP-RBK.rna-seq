#!/bin/bash

bdbag --resolve-fetch all --fetch-filter filename\$*fastq.gz $1
for i in $(find */ -name "*.R*.fastq.gz"); do
    path=${2}$(echo ${i##*/} | grep -o "\.R.\.fastq\.gz");
    mv ${i} ./${path}
done;