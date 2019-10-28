#!/usr/bin/env nextflow

// Define input variables
params.deriva = "${baseDir}/../test_data/deriva-cookies.txt"
params.bdbag = "${baseDir}/../test_data/Study_Q-Y4H0.zip"

params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = file(params.deriva)
deriva.copyTo('~/.bdbag/deriva-cookies.txt')
//deriva = Channel
//  .fromPath(params.deriva)
//  .ifEmpty { exit 1, "deriva cookie file not found: ${params.deriva}" }
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "bdbag zip file not found: ${params.bdbag}" }

outDir = params.outDir

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel
 */
process splitData {
  tag "${bdbag.baseName}"
  publishDir "${outDir}/temp/${task.process}", mode: "symlink"

  input:
    file bdbag
    file deriva

  output:
    file("Replicate_*.zip") into bdbagSplit mode flatten
    file("${bdbag.baseName}/data/File.csv") into fileMeta
    file("${bdbag.baseName}/data/Experiment Settings.csv") into experimentSettingsMeta
    file("${bdbag.baseName}/data/Experiment.csv") into experimentMeta

  script:
    """
    hostname
    ulimit -a
    ln -sf `readlink -e ${deriva}` ~/.bdbag/deriva-cookies.txt
    study=`echo "${bdbag}" | cut -d'.' -f1`
    echo LOG: \${study}
    unzip ${bdbag}
    echo LOG: bdgag unzipped
    python3 ${baseDir}/scripts/modifyFetch.py --fetchFile \${study}
    cd \${study}
    bash ${baseDir}/scripts/fixFetch.sh
    cd ..
    echo LOG: fetch file filtered for only .fastq.gz
    python3 ${baseDir}/scripts/splitFetch.py --fetchFile \${study}
    echo LOG: fetch file split by replicates
    sh ${baseDir}/scripts/splitBag.sh \${study}
    echo LOG: bag recreated with replicate split fetch file
    """
}

println {${http_proxy}}
println {${https_proxy}}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${rep.baseName}"
  publishDir "${outDir}/temp/${task.process}", mode: "symlink"

  input:
    file deriva
    each rep from bdbagSplit

  output:
    file("**/*.R*.fastq.gz") into fastq

  script:
    """
    hostname
    ulimit -a
    echo LOG:\${http_proxy}
    export https_proxy=\${http_proxy}
    ln -sf `readlink -e ${deriva}` ~/.bdbag/deriva-cookies.txt
    replicate=\$(echo "${rep}" | cut -d'.' -f1 | rev | cut -f1 -d '/' | rev)
    echo LOG: \${replicate}
    echo LOG: deriva cookie loaded
    unzip ${rep}
    echo LOG: replicate bdbag unzipped
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate}
    echo LOG: replicate bdbag fetched
    sh ${baseDir}/scripts/renameFastq.sh \${replicate}
    echo LOG: fastq.gz files renamed to replicate RID
    """
 }
