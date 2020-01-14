#!/usr/bin/env nextflow

// Define input variables
params.deriva = "${baseDir}/../test_data/credential.json"
params.bdbag = "${baseDir}/../test_data/cookies.txt"
//params.repRID = "16-1ZX4"
params.repRID = "Q-Y5JA"
//params.bdbag = "${baseDir}/../test_data/Study_Q-Y4H0.zip"

params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "deriva cookie file for bdbag not found: ${params.bdbag}" }
//bdbag = Channel
//  .fromPath(params.bdbag)
//  .ifEmpty { exit 1, "bdbag zip file not found: ${params.bdbag}" }

Channel.from(params.repRID)
  .into {
    repRID_getBag
    repRID_getData
    repRID_trimData
  }

outDir = params.outDir
logsDir = "${outDir}/Logs"

derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
 /*
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
*/


/*
 * getData: get bagit file from consortium
 */
process getBag {
  executor 'local'
  tag "${repRID_getBag}"
  publishDir "${logsDir}/getBag", mode: 'symlink', pattern: "${repRID_getBag}.getBag.err"

  input:
    val repRID_getBag
    path credential, stageAs: 'credential.json' from deriva
    path derivaConfig

  output:
    path ("Replicate_*.zip") into bagit

  script:
    """
    hostname >>${repRID_getBag}.getBag.err
    ulimit -a >>${repRID_getBag}.getBag.err
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json 2>>${repRID_getBag}.getBag.err
    echo "LOG: deriva credentials linked" >>${repRID_getBag}.getBag.err
    deriva-download-cli dev.gudmap.org --catalog 2 ${derivaConfig} . rid=${repRID_getBag} 2>>${repRID_getBag}.getBag.err
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID_getData}"
  publishDir "${logsDir}/getData", mode: 'symlink', pattern: "${repRID_getData}.getData.err"

  input:
    val repRID_getData
    path cookies, stageAs: 'deriva-cookies.txt' from bdbag
    path bagit

  output:
    path ("**/*.R{1,2}.fastq.gz") into fastqs

  script:
    """
    hostname >>${repRID_getData}.getData.err
    ulimit -a >>${repRID_getData}.getData.err
    export https_proxy=\${http_proxy}
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt >>${repRID_getData}.getData.err
    echo "LOG: deriva cookie linked" >>${repRID_getData}.getData.err
    replicate=\$(basename "${bagit}" | cut -d '.' -f1)
    echo "LOG: \${replicate}" >>${repRID_getData}.getData.err
    unzip ${bagit} 2>>${repRID_getData}.getData.err
    echo "LOG: replicate bdbag unzipped" >>${repRID_getData}.getData.err
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate} 2>>${repRID_getData}.getData.err
    echo "LOG: replicate bdbag fetched" >>${repRID_getData}.getData.err
    """
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID_trimData}"
  publishDir "${outDir}/tempOut/trimmed", mode: "symlink", pattern: "*_val_{1,2}.fq.gz"
  publishDir "${logsDir}/trimData", mode: 'symlink', pattern: "\${repRID_trimData}.trimData.*"

  input:
    val repRID_trimData
    file(fastq) from fastqs

  output:
    path ("*_val_{1,2}.fq.gz", type: 'file', maxDepth: '0')

  script:
    """
    if [ `nproc` -gt 8 ]
    then
      ncore=8
    else
      ncore=`nproc`
    fi
    if [ -z ${fastq[1]} ]
    then
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID_trimData} -j \${ncore} ${fastq[0]} 1>>${repRID_trimData}.trimData.log 2>>${repRID_trimData}.trimData.err;
    else
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID_trimData} -j \${ncore} ${fastq[0]} ${fastq[1]} 1>>${repRID_trimData}.trimData.log 2>>${repRID_trimData}.trimData.err;
    fi
    """
}