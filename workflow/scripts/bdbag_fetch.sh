#!/bin/bash

bdbag --materialize ${1} --debug
validate=""
bdbag --validate full ${1} 2> validate.txt
validate=$(tail -n1 validate.txt | grep -o 'is valid')
if [ "${validate}" != "is valid" ]
then
    n=0
    until [ "${n}" -ge "3" ]
    do
        bdbag --resolve-fetch missing --validate full ${1} --debug --config-file bdbag.json && validate=$(tail -n validate.txt | grep -o 'is valid') && break
        n=$((n+1)) 
        sleep 15
    done
fi
if [ "${validate}" != "is valid" ]
then
    exit 1
fi
for i in $(find */ -name "*[_.]R[1-2].fastq.gz")
do
    path=${2}.$(echo ${i##*/} | grep -o "R[1,2].fastq.gz")
    cp ${i} ./${path}
done
