#!/bin/bash

echo "collecting stats for badges"

commits=`git rev-list --all --count`
latest_release_tag=$(git describe --tags `git rev-list --tags --max-count=1`)
current_pipeline_version=$(cat ./workflow/nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_nextflow_version=$(cat ./workflow/nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
echo "{\"commits\":\"$commits\", \"release_tag\":\"$latest_release_tag\", \"pipeline\":\"$current_pipeline_version\", \"nextflow\":\"$current_nextflow_version\"}" > badges.json