#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
echo "latest_release_tag =" ${latest_release_tag}
current_pipeline_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
echo "current_pipeline_version =" ${current_pipeline_version}
current_nextflow_version=$(git show ${latest_release_tag}:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
echo "current_nextflow_version =" ${current_nextflow_version}
master_pipeline_version=$(git show origin/master:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
echo "master_pipeline_version =" ${master_pipeline_version}
master_nextflow_version=$(git show origin/master:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
echo "master_nextflow_version =" ${master_nextflow_version}
develop_pipeline_version=$(git show origin/develop:workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
echo "develop_pipeline_version =" ${develop_pipeline_version}
develop_nextflow_version=$(git show origin/develop:workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
echo "develop_nextflow_version =" ${develop_nextflow_version}

echo "collecting badges"
mkdir badges
curl --request GET https://img.shields.io/badge/Latest%20Release-${latest_release_tag}-informational?style=flat > ./badges/release.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${current_pipeline_version}-informational?style=flat > ./badges/releasePipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${current_nextflow_version}-informational?style=flat > ./badges/releaseNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${master_pipeline_version}-informational?style=flat > ./badges/masterPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${master_nextflow_version}-informational?style=flat > ./badges/masterNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${develop_pipeline_version}-informational?style=flat > ./badges/developPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${develop_nextflow_version}-informational?style=flat > ./badges/developNextflow.svg