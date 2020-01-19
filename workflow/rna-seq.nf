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
repRID = params.repRID

outDir = params.outDir
logsDir = "${outDir}/Logs"

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbagFetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parseMeta.py")

/*
 * getData: get bagit file from consortium
 */
process getBag {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.getBag.err"

  input:
    path credential, stageAs: 'credential.json' from deriva
    path derivaConfig

  output:
    path ("Replicate_*.zip") into bagit
    file ("${repRID}.getBag.err")

  script:
    """
    hostname >>${repRID}.getBag.err
    ulimit -a >>${repRID}.getBag.err
    export https_proxy=\${http_proxy}

    # link credential file for authentication
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json 2>>${repRID}.getBag.err
    echo "LOG: deriva credentials linked" >>${repRID}.getBag.err

    # deriva-download replicate RID
    deriva-download-cli dev.gudmap.org --catalog 2 ${derivaConfig} . rid=${repRID} 2>>${repRID}.getBag.err
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.getData.err"

  input:
    path script_bdbagFetch
    path cookies, stageAs: 'deriva-cookies.txt' from bdbag
    path bagit

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    file("**/File.csv") into fileMeta
    file("**/Experiment Settings.csv") into experimentSettingsMeta
    file("**/Experiment.csv") into experimentMeta
    file ("${repRID}.getData.err")

  script:
    """
    hostname >>${repRID}.getData.err
    ulimit -a >>${repRID}.getData.err
    export https_proxy=\${http_proxy}
    
    # link deriva cookie for authentication
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt >>${repRID}.getData.err
    echo "LOG: deriva cookie linked" >>${repRID}.getData.err
    
    # get bagit basename
    replicate=\$(basename "${bagit}" | cut -d '.' -f1)
    echo "LOG: \${replicate}" >>${repRID}.getData.err
    
    # unzip bagit
    unzip ${bagit} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag unzipped" >>${repRID}.getData.err
    
    # bagit fetch fastq's only and rename by repRID
    sh ${script_bdbagFetch} \${replicate} ${repRID} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag fetched" >>${repRID}.getData.err
    """
}

/*
 * parseMetadata: parses metadata to extract experiment parameters
*/
process parseMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.parseMetadata.err"

  input:
    path script_parseMeta
    val repRID
    path fileMeta
    path experimentSettingsMeta
    path experimentMeta

  output:
    val ends
    val stranded
    val spike
    val specie

  script:
    """
    hostname >>${repRID}.parseMetadata.err
    ulimit -a >>${repRID}.parseMetadata.err

    # Check replicate RID metadata
    rep=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p repRID)
    echo "LOG: replicate RID metadata parsed: \${rep}" >>${repRID}.parseMetadata.err
    
    # Get endedness metadata
    ends=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p ends)
    echo "LOG: endedness metadata parsed: \${ends}" >>${repRID}.parseMetadata.err
    
    # Get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p stranded)
    echo "LOG: strandedness metadata parsed: \${stranded}" >>${repRID}.parseMetadata.err
    
    # Get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p spike)
    echo "LOG: spike-in metadata parsed: \${spike}" >>${repRID}.parseMetadata.err
    
    # Get species metadata
    specie=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentMeta}" -p specie)
    echo "LOG: species metadata parsed: \${specie}" >>${repRID}.parseMetadata.err
    """
}

ends.set {
  ends_trimData
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "\${repRID}.trimData.*"

  input:
    file(fastq) from fastqs
    val ends_trimData

  output:
    path ("*.fq.gz") into fastqs_trimmed
    val ends
    file ("${repRID}.trimData.log")
    file ("${repRID}.trimData.err")

  script:
    """
    hostname >>${repRID}.trimData.err
    ulimit -a >>${repRID}.trimData.err

    # trim fastqs
    if [ '${ends_trimData}' == 'se' ]
    then
      ends='se'
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    else
      ends='pe'
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    fi
    """
}