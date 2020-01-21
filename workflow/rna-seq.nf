#!/usr/bin/env nextflow

// Define input variables
params.deriva = "${baseDir}/../test_data/credential.json"
params.bdbag = "${baseDir}/../test_data/cookies.txt"
//params.repRID = "16-1ZX4"
params.repRID = "Q-Y5JA"
params.spikein = "false"
params.reference = "mouse"

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
reference = params.reference

if (params.spikein) {
  if (params.reference == "human") {
    reference = file ("/project/BICF/BICF_Core/s181706/github/gudmap/rna-seq/References/GRCh38.p12-S/hisat2")
  } else if (params.reference == "mouse") {
    reference = file ("/project/BICF/BICF_Core/s181706/github/gudmap/rna-seq/References/GRCm38.P6-S/hisat2")
  } 
} else if (params.reference == "mouse") {
  reference = file ("/project/BICF/BICF_Core/s181706/github/gudmap/rna-seq/References/GRCm38.P6/hisat2")
} else if (params.reference == "human") {
  reference = file ("/project/BICF/BICF_Core/s181706/github/gudmap/rna-seq/References/GRCh38.p12/hisat2")
} else {
  print ("Warning: Reference genome not specified, defaulting to GRCm38.P6 with NO spike-in")
  reference = file ("/project/BICF/BICF_Core/s181706/github/gudmap/rna-seq/References/GRCm38.P6/hisat2")
}

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")

/*
 * getData: get bagit file from consortium
 */
process getBag {
  executor 'local'
  tag "${repRID}"
  publishDir "${logsDir}/getBag", mode: 'copy', pattern: "${repRID}.getBag.err"

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
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json 2>>${repRID}.getBag.err
    echo "LOG: deriva credentials linked" >>${repRID}.getBag.err
    deriva-download-cli dev.gudmap.org --catalog 2 ${derivaConfig} . rid=${repRID} 2>>${repRID}.getBag.err
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID}"
  publishDir "${logsDir}/getData", mode: 'copy', pattern: "${repRID}.getData.err"

  input:
    executor 'local'
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
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt >>${repRID}.getData.err
    echo "LOG: deriva cookie linked" >>${repRID}.getData.err
    replicate=\$(basename "${bagit}" | cut -d '.' -f1)
    echo "LOG: \${replicate}" >>${repRID}.getData.err
    unzip ${bagit} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag unzipped" >>${repRID}.getData.err
    sh ${baseDir}/scripts/bdbagFetch.sh \${replicate} ${repRID} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag fetched" >>${repRID}.getData.err
    """
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"
  publishDir "${logsDir}/trimData", mode: 'copy', pattern: "\${repRID}.trimData.*"

  input:
    file(fastq) from fastqs

  output:
    tuple file ("params.csv"), file ("*.fq.gz") into fastqs_trimmed
    file ("${repRID}.trimData.log")
    file ("${repRID}.trimData.err")

  script:
    """
    if [ '${fastq[1]}' == 'null' ]
    then
      ends='se'
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    else
      ends='pe'
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    fi
    echo \${ends} > params.csv;
    """
}

/*
 * alignReads: aligns the reads to a reference database
*/
process alignReads {
  tag "align-${repRID}"
  publishDir "${outDir}/tempOut/aligned", mode: "copy"

  input:
    set ends, fqs from fastqs_trimmed
    file reference

  output:
    set repRID, file ("${repRID}.unal.gz"), file ("${repRID}.bam"), file ("${repRID}.bai")

  script:
    """
    ends=`cat ${ends}`;
    if [ "\${ends}" == 'pe' ]; then
    hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x ${reference}/genome -1 ${fqs[0]} -2 ${fqs[1]} 1>${repRID}.align.out 2> ${repRID}.align.err;
    else hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x ${reference}/genome -U ${fqs[0]} 1>${repRID}.align.out 2> ${repRID}.align.err;
    fi;
    samtools view -1 --threads `nproc` -o ${repRID}.bam ${repRID}.sam 1>${repRID}.align.out 2> ${repRID}.align.err;
    samtools sort -@ `nproc` -O BAM ${repRID}.bam 1>${repRID}.align.out 2> ${repRID}.align.err;
    """
}

