#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
current_pipeline_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_pipeline_versionMajor=$(echo ${current_pipeline_version} | cut -f1 -d".")
current_pipeline_versionMajor=$(echo ${current_pipeline_versionMajor}".")
echo "Major pipeline version for search: "${current_pipeline_versionMajor}

echo "collecting workflow RIDs from servers"
dev_workflow_RID=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")
staging_workflow_RID=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")
prod_workflow_RID=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")

echo "collecting unique replicates with successful execution runs"
dev_total=0
dev_success=0
dev_error=0
for rid in ${dev_workflow_RID}
do
    temp_total=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/\!Execution_Status=In-progress/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    dev_total=$(expr ${dev_total} + ${temp_total})
    temp_sucess=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Success/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    dev_success=$(expr ${dev_success} + ${temp_success})
    temp_error=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Error/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    dev_error=$(expr ${dev_error} + ${temp_error})
done
staging_total=0
staging_success=0
staging_error=0
for rid in ${staging_workflow_RID}
do
    temp_total=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/\!Execution_Status=In-progress/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    staging_total=$(expr ${staging_total} + ${temp_total})
    temp_sucess=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Success/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    staging_success=$(expr ${staging_success} + ${temp_success})
    temp_error=$(curl -s https://dstagingev.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Error/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    staging_error=$(expr ${staging_error} + ${temp_error})
done
prod_total=0
prod_success=0
prod_error=0
for rid in ${prod_workflow_RID}
do
    temp_total=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/\!Execution_Status=In-progress/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    prod_total=$(expr ${prod_total} + ${temp_total})
    temp_sucess=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Success/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    prod_success=$(expr ${prod_success} + ${temp_success})
    temp_error=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/Q:=RNASeq:mRNA_QC/E:=\(Execution_Run\)=\(RNASeq:Execution_Run:RID\)/Workflow=${rid}/Execution_Status=Error/\$Q | grep -o \"Replicate\".*,\"Paired_End | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Paired_End)" | sort | uniq | wc -l)
    prod_error=$(expr ${prod_error} + ${temp_error})
done

echo "collecting badges"
mkdir -p ./badges/counts
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${dev_success}/${dev_total}-green?style=flat > ./badges/counts/dev_success.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${dev_error}/${dev_total}-red?style=flat > ./badges/counts/dev_error.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${staging_success}/${staging_total}-green?style=flat > ./badges/counts/staging_success.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${staging_error}/${staging_total}-red?style=flat > ./badges/counts/staging_error.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${prod_success}/${prod_total}-green?style=flat > ./badges/counts/prod_success.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${prod_error}/${prod_total}-red?style=flat > ./badges/counts/prod_error.svg

curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${dev_success}-lightgrey?style=flat > ./badges/counts/dev_counts.svg
curl --request GET https://img.shields.io/badge/Staging%20Replicate%20Count-${staging_success}-lightgrey?style=flat > ./badges/counts/staging_counts.svg
curl --request GET https://img.shields.io/badge/Production%20Replicate%20Count-${prod_success}-lightgrey?style=flat > ./badges/counts/prod_counts.svg
