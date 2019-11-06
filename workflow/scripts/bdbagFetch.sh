#!/bin/bash

bdbag --resolve-fetch all --fetch-filter filename\$*fastq.gz $1 &&
for i in $(find */ -name "*.R*.fastq.gz"); do
  mv ${i} .;
done;
