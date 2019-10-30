#!/usr/bin/env nextflow

// Define input variables
params.deriva = "/project/BICF/BICF_Core/shared/gudmap/cookies/deriva-cookies.txt"
params.bdbag = "${baseDir}/../test_data/Study_Q-Y4H0.zip"

params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = file(params.deriva, checkIfExists: 'true')
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "bdbag zip file not found: ${params.bdbag}" }

outDir = params.outDir

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
process splitData {
  tag "${bdbag.baseName}"

  input:
    file bdbag
    path cookies, stageAs: 'cookies.txt' from deriva

  output:
    file("Replicate_*.zip") into bdbagSplit mode flatten
    file("${bdbag.baseName}/data/File.csv") into fileMeta
    file("${bdbag.baseName}/data/Experiment Settings.csv") into experimentSettingsMeta
    file("${bdbag.baseName}/data/Experiment.csv") into experimentMeta

  script:
    """
    hostname
    ulimit -a
    ln -sf `readlink -e cookies.txt` ~/.bdbag/deriva-cookies.txt
    echo "LOG: deriva cookie linked"
    study=`echo "${bdbag}" | cut -d '.' -f1`
    echo "LOG: \${study}"
    unzip ${bdbag}
    echo "LOG: bdgag unzipped"
    python3 ${baseDir}/scripts/modifyFetch.py --fetchFile \${study}
    echo "LOG: fetch file filtered for only .fastq.gz"
    python3 ${baseDir}/scripts/splitFetch.py --fetchFile \${study}
    echo "LOG: fetch file split by replicates"
    sh ${baseDir}/scripts/splitBag.sh \${study}
    echo "LOG: bag recreated with replicate split fetch file"
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${rep.baseName}"
  publishDir "${outDir}/tempOut/fastqs", mode: "symlink"

  input:
    each rep from bdbagSplit

  output:
    path ("*.R*.fastq.gz", type: 'file', maxDepth: '0') into fastq

  script:
    """
    hostname
    ulimit -a
    export https_proxy=\${http_proxy}
    replicate=\$(basename "${rep}" | cut -d '.' -f1)
    echo "LOG: \${replicate}"
    unzip ${rep}
    echo "LOG: replicate bdbag unzipped"
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate}
    echo "LOG: replicate bdbag fetched"
    """
 }
