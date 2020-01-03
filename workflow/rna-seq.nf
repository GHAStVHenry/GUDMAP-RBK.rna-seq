#!/usr/bin/env nextflow

// Define input variables
params.deriva = "/project/BICF/BICF_Core/shared/gudmap/cookies/deriva-cookies.txt"
params.bdbag = "${baseDir}/../test_data/Study_Q-Y4H0.zip"
params.spikein = "false"
params.species = "human"

params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = file(params.deriva, checkIfExists: 'true')
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "bdbag zip file not found: ${params.bdbag}" }

outDir = params.outDir
logsDir = "${outDir}/Logs"

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
process splitData {
  tag "${bdbag.baseName}"
  executor 'local'
  publishDir "${logsDir}/splitData", mode: 'symlink', pattern: "${bdbag.baseName}.splitData.err"

  input:
    file bdbag
    path cookies, stageAs: 'cookies.txt' from deriva

  output:
    file("Replicate_*.zip") into bdbagSplit mode flatten
    file("${bdbag.baseName}/data/File.csv") into fileMeta
    file("${bdbag.baseName}/data/Experiment Settings.csv") into experimentSettingsMeta
    file("${bdbag.baseName}/data/Experiment.csv") into experimentMeta
    file ("${bdbag.baseName}.splitData.err")

  script:
    """
    hostname >> ${bdbag.baseName}.splitData.err
    ulimit -a >> ${bdbag.baseName}.splitData.err
    ln -sf `readlink -e cookies.txt` ~/.bdbag/deriva-cookies.txt 2>>${bdbag.baseName}.splitData.err
    echo "LOG: deriva cookie linked" >> ${bdbag.baseName}.splitData.err 
    study=`echo "${bdbag}" | cut -d '.' -f1` 2>>${bdbag.baseName}.splitData.err
    echo "LOG: \${study}" >> ${bdbag.baseName}.splitData.err
    unzip ${bdbag} 2>>${bdbag.baseName}.splitData.err
    echo "LOG: bdgag unzipped" >> ${bdbag.baseName}.splitData.err
    python3 ${baseDir}/scripts/modifyFetch.py --fetchFile \${study} 2>>${bdbag.baseName}.splitData.err
    echo "LOG: fetch file filtered for only .fastq.gz" >> ${bdbag.baseName}.splitData.err
    python3 ${baseDir}/scripts/splitFetch.py --fetchFile \${study} 2>>${bdbag.baseName}.splitData.err
    echo "LOG: fetch file split by replicates" >> ${bdbag.baseName}.splitData.err
    sh ${baseDir}/scripts/splitBag.sh \${study} 2>>${bdbag.baseName}.splitData.err
    echo "LOG: bag recreated with replicate split fetch file" >> ${bdbag.baseName}.splitData.err
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${rep.baseName}"
  publishDir "${logsDir}/getData", mode: 'symlink', pattern: "${rep.baseName}.getData.err"

  input:
    each rep from bdbagSplit

  output:
    set val ("${rep.baseName}"), file ("*.R{1,2}.fastq.gz") into trimming

  script:
    """
    hostname >>${rep.baseName}.getData.err
    ulimit -a >>${rep.baseName}.getData.err
    export https_proxy=\${http_proxy}
    replicate=\$(basename "${rep}" | cut -d '.' -f1)
    echo "LOG: \${replicate}" >>${rep.baseName}.getData.err
    unzip ${rep} 2>>${rep.baseName}.getData.err
    echo "LOG: replicate bdbag unzipped" >>${rep.baseName}.getData.err
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate} 2>>${rep.baseName}.getData.err
    echo "LOG: replicate bdbag fetched" >>${rep.baseName}.getData.err
    """
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "trim-${repID}"
  publishDir "${outDir}/tempOut/trimmed", mode: "symlink", pattern: "*_val_{1,2}.fq.gz"
  publishDir "${logsDir}/trimData", mode: 'symlink', pattern: "\${rep}.trimData.*"

  input:
    set repID, reads from trimming

  output:
    set repID, path ("*_val_{1,2}.fq.gz", type: 'file', maxDepth: '0') into aligning

  script:
    """
    rep=`echo ${repID} | cut -f2- -d '_'`;
    trim_galore --gzip --max_n 1 --paired --basename \${rep} -j `nproc` ${reads[0]} ${reads[1]} 1>>\${rep}.trimData.log 2>>\${rep}.trimData.err;
    """
}

/*
 * alignReads: aligns the reads to a reference database
*/
process alignReads {
  tag "align-${repID}"
  publishDir "${outDir}/tempOut/aligned", mode: "symlink"

  input:
    set repID, fqs from aligning

  
