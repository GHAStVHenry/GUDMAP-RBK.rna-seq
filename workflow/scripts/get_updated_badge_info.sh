#!/bin/bash

echo "collecting stats for badges"
latest_release_tag=$(git tag --sort=-committerdate -l *.*.* | head -1)
current_pipeline_version=$(git show ${latest_release_tag}:nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
current_nextflow_version=$(git show ${latest_release_tag}:nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
master_pipeline_version=$(git show origin/master:nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
master_nextflow_version=$(git show origin/master:nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")
develop_pipeline_version=$(git show origin/develop:nextflow.config | grep -o version.* | grep -oP "(?<=').*(?=')")
develop_nextflow_version=$(git show origin/develop:nextflow.config | grep -o nextflowVersion.* | grep -oP "(?<=').*(?=')")

echo "collecting tool version for badges"
python_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o Python.* | grep -oP "(?<=d>).*(?=\<)")
deriva_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o DERIVA.* | grep -oP "(?<=d>).*(?=\<)")
bdbag_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o BDBag.* | grep -oP "(?<=d>).*(?=\<)")
trimgalore_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o 'Trim Galore!'.* | grep -oP "(?<=d>).*(?=\<)")
hisat2_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o HISAT2.* | grep -oP "(?<=d>).*(?=\<)")
samtools_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o Samtools.* | grep -oP "(?<=d>).*(?=\<)")
picard_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o 'picard (MarkDuplicates)'.* | grep -oP "(?<=d>).*(?=\<)")
featurecounts_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o featureCounts.* | grep -oP "(?<=d>).*(?=\<)")
deeptools_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o deepTools.* | grep -oP "(?<=d>).*(?=\<)")
seqtk_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o Seqtk.* | grep -oP "(?<=d>).*(?=\<)")
r_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o '>R<'.* | grep -oP "(?<=d>).*(?=\<)")
fastqc_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o FastQC.* | grep -oP "(?<=d>).*(?=\<)")
seqwho_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o SeqWho.* | grep -oP "(?<=d>).*(?=\<)")
rseqc_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o RSeQC.* | grep -oP "(?<=d>).*(?=\<)")
multiqc_version=$(git show ${latest_release_tag}:docs/software_versions_mqc.yaml | grep -o MultiQC.* | grep -oP "(?<=d>).*(?=\<)")

echo "collecting badges"
mkdir -p ./badges/tools
curl --request GET https://img.shields.io/badge/Latest%20Release-${latest_release_tag}-informational?style=flat > ./badges/release.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${current_pipeline_version}-informational?style=flat > ./badges/releasePipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${current_nextflow_version}-informational?style=flat > ./badges/releaseNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${master_pipeline_version}-informational?style=flat > ./badges/masterPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${master_nextflow_version}-informational?style=flat > ./badges/masterNextflow.svg
curl --request GET https://img.shields.io/badge/Pipeline%20Version-${develop_pipeline_version}-informational?style=flat > ./badges/developPipeline.svg
curl --request GET https://img.shields.io/badge/Nextflow%20Version-${develop_nextflow_version}-informational?style=flat > ./badges/developNextflow.svg

curl --request GET https://img.shields.io/badge/Python%20Version-${python_version}-blueviolet?style=flat > ./badges/tools/python.svg
curl --request GET https://img.shields.io/badge/DERIVA%20Version-${deriva_version}-blueviolet?style=flat > ./badges/tools/deriva.svg
curl --request GET https://img.shields.io/badge/BDBag%20Version-${bdbag_version}-blueviolet?style=flat > ./badges/tools/bdbag.svg
curl --request GET https://img.shields.io/badge/Trim%20Galore%20Version-${trimgalore_version}-blueviolet?style=flat > ./badges/tools/trimgalore.svg
curl --request GET https://img.shields.io/badge/HISAT2%20Version-${hisat2_version}-blueviolet?style=flat > ./badges/tools/hisat2.svg
curl --request GET https://img.shields.io/badge/Samtools%20Version-${samtools_version}-blueviolet?style=flat > ./badges/tools/samtools.svg
curl --request GET https://img.shields.io/badge/picard%20Version-${picard_version}-blueviolet?style=flat > ./badges/tools/picard.svg
curl --request GET https://img.shields.io/badge/featureCounts%20Version-${featurecounts_version}-blueviolet?style=flat > ./badges/tools/featurecounts.svg
curl --request GET https://img.shields.io/badge/deepTools%20Version-${deeptools_version}-blueviolet?style=flat > ./badges/tools/deeptools.svg
curl --request GET https://img.shields.io/badge/Seqtk%20Version-${seqtk_version}-blueviolet?style=flat > ./badges/tools/seqtk.svg
curl --request GET https://img.shields.io/badge/R%20Version-${r_version}-blueviolet?style=flat > ./badges/tools/r.svg
curl --request GET https://img.shields.io/badge/FastQC%20Version-${fastqc_version}-blueviolet?style=flat > ./badges/tools/fastqc.svg
curl --request GET https://img.shields.io/badge/SeqWho%20Version-${seqwho_version}-blueviolet?style=flat > ./badges/tools/seqwho.svg
curl --request GET https://img.shields.io/badge/RSeQC%20Version-${rseqc_version}-blueviolet?style=flat > ./badges/tools/rseqc.svg
curl --request GET https://img.shields.io/badge/MultiQC%20Version-${multiqc_version}-blueviolet?style=flat > ./badges/tools/multiqc.svg

echo "creating blank env badges if not tested"
if [ ! -f ./badges/env/dnanexus.svg ]
then
curl --request GET https://img.shields.io/badge/Envronment%3A%20DNAnexus-not%2Dtested-important?style=flat > ./badges/env/dnanexus.svg
fi
if [ ! -f ./badges/env/aws.svg ]
then
curl --request GET https://img.shields.io/badge/Envronment%3A%20AWS-not%2Dtested-important?style=flat > ./badges/env/aws.svg
fi
if [ ! -f ./badges/env/azure.svg ]
then
curl --request GET https://img.shields.io/badge/Envronment%3A%20Azure-not%2Dtested-important?style=flat > ./badges/env/azure.svg
fi
if [ ! -f ./badges/env/gcp.svg ]
then
curl --request GET https://img.shields.io/badge/Envronment%3A%20GCP-not%2Dtested-important?style=flat > ./badges/env/gcp.svg
fi