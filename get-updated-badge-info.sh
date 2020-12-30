#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git describe --tags `git rev-list --tags --max-count=1`)
current_pipeline_version=$(cat ./workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_nextflow_version=$(cat ./workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")

echo "collecting badges"
mkdir badges
curl --request GET https://img.shields.io/badge/Latest%20Release-${latest_release_tag}-green?style=plastic > ./badges/release.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${current_pipeline_version}-green?style=plastic > ./badges/pipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${current_nextflow_version}-green?style=plastic > ./badges/nextflow.svg