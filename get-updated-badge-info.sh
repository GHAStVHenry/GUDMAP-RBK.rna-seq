#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
current_pipeline_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_nextflow_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
master_pipeline_version=$(git show master:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
master_nextflow_version=$(git show master:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
develop_pipeline_version=$(git show develop:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
develop_nextflow_version=$(git show develop:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")

echo "collecting badges"
mkdir badges
curl --request GET https://img.shields.io/badge/Latest%20Release-${latest_release_tag}-green?style=plastic > ./badges/release.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${current_pipeline_version}-green?style=plastic > ./badges/releasePipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${current_nextflow_version}-green?style=plastic > ./badges/releaseNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${master_pipeline_version}-green?style=plastic > ./badges/masterPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${master_nextflow_version}-green?style=plastic > ./badges/masterNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${develop_pipeline_version}-green?style=plastic > ./badges/developPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${develop_nextflow_version}-green?style=plastic > ./badges/developNextflow.svg