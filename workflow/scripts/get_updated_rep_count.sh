#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
current_pipeline_version=$(git show ${latest_release_tag}:nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_pipeline_versionMajor=$(echo ${current_pipeline_version} | cut -f1 -d".")
current_pipeline_versionMajor=$(echo ${current_pipeline_versionMajor}".")
echo "Major pipeline version for search: "${current_pipeline_versionMajor}

echo "collecting workflow RIDs from servers"
dev_workflow_RID=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")
staging_workflow_RID=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")
prod_workflow_RID=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version::ciregexp::%5E${current_pipeline_versionMajor} | grep -o '\"RID\":\".*\",\"RCT' | cut -f4 -d"\"")

echo "collecting unique replicates with successful execution runs"
dev_count=0
for rid in ${dev_workflow_RID}
do
    temp_count=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    dev_success=$(expr ${dev_success} + ${temp_count})
    temp_count=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Error/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    dev_error=$(expr ${dev_error} + ${temp_count})
done
staging_count=0
for rid in ${staging_workflow_RID}
do
    temp_count=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    staging_success=$(expr ${staging_success} + ${temp_count})
    temp_count=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Error/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    staging_error=$(expr ${staging_error} + ${temp_count})
done
prod_count=0
for rid in ${prod_workflow_RID}
do
    temp_count=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    prod_success=$(expr ${prod_success} + ${temp_count})
    temp_count=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Error/Workflow=${rid} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
    prod_error=$(expr ${prod_error} + ${temp_count})
done

echo "collecting badges"
mkdir -p ./badges/counts
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Success-${dev_success}-lightgrey?style=flat > ./badges/counts/dev_success.svg
curl --request GET https://img.shields.io/badge/Staging%20Replicate%20Success-${staging_success}-lightgrey?style=flat > ./badges/counts/staging_success.svg
curl --request GET https://img.shields.io/badge/Production%20Replicate%20Success-${prod_success}-lightgrey?style=flat > ./badges/counts/prod_success.svg
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Error-${dev_error}-lightgrey?style=flat > ./badges/counts/dev_error.svg
curl --request GET https://img.shields.io/badge/Staging%20Replicate%20Error-${staging_error}-lightgrey?style=flat > ./badges/counts/staging_error.svg
curl --request GET https://img.shields.io/badge/Production%20Replicate%20Error-${prod_error}-lightgrey?style=flat > ./badges/counts/prod_error.svg
