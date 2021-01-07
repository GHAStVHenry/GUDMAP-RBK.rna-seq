#!/bin/bash

bdbag --materialize ${1} --debug
validateError=true
bdbag --validate full ${1} && validateError=false
if validateError
then
    n=0
    until [ "${n}" -ge "3" ]
    do
        bdbag --resolve-fetch missing --validate full ${1} --debug && validateError=false && break
        n=$((n+1)) 
        sleep 15
    done
fi
if validateError
then
    exit 1
fi
for i in $(find */ -name "*R*.fastq.gz")
do
    path=${2}.$(echo ${i##*/} | grep -o "R[1,2].fastq.gz")
    cp ${i} ./${path}
done