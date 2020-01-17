#!/usr/bin/env nextflow

// Define input variables
params.deriva = "${baseDir}/../test_data/credential.json"
params.bdbag = "${baseDir}/../test_data/cookies.txt"
//params.repRID = "16-1ZX4"
params.repRID = "Q-Y5JA"

params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "deriva cookie file for bdbag not found: ${params.bdbag}" }

Channel.from(params.repRID)
  .into {
    repRID_getBag
    repRID_getData
    repRID_trimData
  }

outDir = params.outDir
logsDir = "${outDir}/Logs"

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")

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
    file ("${repRID_getBag}.getBag.err")

  script:
    """
    hostname >>${repRID_getBag}.getBag.err
    ulimit -a >>${repRID_getBag}.getBag.err
    export https_proxy=\${http_proxy}
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
    executor 'local'
    path cookies, stageAs: 'deriva-cookies.txt' from bdbag
    path bagit

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    file("**/File.csv") into fileMeta
    file("**/Experiment Settings.csv") into experimentSettingsMeta
    file("**/Experiment.csv") into experimentMeta
    file ("${repRID_getData}.getData.err")


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
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate} ${repRID_getData} 2>>${repRID_getData}.getData.err
    echo "LOG: replicate bdbag fetched" >>${repRID_getData}.getData.err
    """
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID_trimData}"
  publishDir "${logsDir}/trimData", mode: 'symlink', pattern: "\${repRID_trimData}.trimData.*"

  input:
    val repRID_trimData
    file(fastq) from fastqs

  output:
    path ("*.fq.gz") into fastqs_trimmed
    val ends
    file ("${repRID_trimData}.trimData.log")
    file ("${repRID_trimData}.trimData.err")

  script:
    """
    if [ `nproc` -gt 8 ]
    then
      ncore=8
    else
      ncore=`nproc`
    fi
    if [ '${fastq[1]}' == 'null' ]
    then
      ends='se'
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID_trimData} -j \${ncore} ${fastq[0]} 1>>${repRID_trimData}.trimData.log 2>>${repRID_trimData}.trimData.err;
    else
      ends='pe'
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID_trimData} -j \${ncore} ${fastq[0]} ${fastq[1]} 1>>${repRID_trimData}.trimData.log 2>>${repRID_trimData}.trimData.err;
    fi
    """
}