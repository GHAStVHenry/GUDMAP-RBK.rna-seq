#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
current_pipeline_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")

echo "collecting workflow RIDs from servers"
dev_workflow_RID=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version=${current_pipeline_version} | grep -o '\"RID\":\".*\",\"RCT')
dev_workflow_RID=${dev_workflow_RID:7:-6}
staging_workflow_RID=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version=${current_pipeline_version} | grep -o '\"RID\":\".*\",\"RCT')
staging_workflow_RID=${staging_workflow_RID:7:-6}
prod_workflow_RID=$(curl -s https://www.gudmap.org/ermrest/catalog/2/entity/RNASeq:Workflow/Version=${current_pipeline_version} | grep -o '\"RID\":\".*\",\"RCT')
prod_workflow_RID=${prod_workflow_RID:7:-6}

echo "collecting unique replicates with successful execution runs"
dev_count=$(curl -s https://dev.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${dev_workflow_RID} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
staging_count=$(curl -s https://staging.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${staging_workflow_RID} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)
prod_count=$(curl -s https://prod.gudmap.org/ermrest/catalog/2/entity/RNASeq:Execution_Run/Execution_Status=Success/Workflow=${prod_workflow_RID} | grep -o \"Replicate\".*,\"Workflow | grep -oP "(?<=\"Replicate\":\").*(?=\",\"Workflow)" | sort | uniq | wc -l)

echo "collecting badges"
mkdir -p ./badges/counts
curl --request GET https://img.shields.io/badge/Development%20Replicate%20Count-${dev_count}-lightgrey?style=flat > ./badges/counts/dev_counts.svg
curl --request GET https://img.shields.io/badge/Staging%20Replicate%20Count-${staging_count}-lightgrey?style=flat > ./badges/counts/staging_counts.svg
curl --request GET https://img.shields.io/badge/Production%20Replicate%20Count-${prod_count}-lightgrey?style=flat > ./badges/counts/prod_counts.svg
