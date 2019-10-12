#!/bin

while read loc checksum fileLocation
do
    file=$(echo ${fileLocation##*/})
    fileName=$(echo ${file%.R*.fastq.gz})
    fileExt=$(echo ${file##${fileName}.})
    while IFS="," read RID Study_RID Experiment_RID Replicate_RID Caption File_Type File_Name URI File_size MD5 GEO_Archival_URL dbGaP_Accession_ID Processed Notes Principal_Investigator Consortium Release_Date RCT RMT Legacy_File_RID GUDMAP_NGF_OID GUDMAP_NGS_OID
    do
        if [ ${file} == ${File_Name} ]
        then
            find . -type f -name ${file} -execdir mv {} ${Replicate_RID}.${fileExt} ';'
        fi
    done < $1/data/File.csv
done < $1/fetch.txt