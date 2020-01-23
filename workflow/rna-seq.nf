#!/usr/bin/env nextflow

// Define input variables
params.deriva = "${baseDir}/../test_data/auth/credential.json"
params.bdbag = "${baseDir}/../test_data/auth/cookies.txt"
//params.repRID = "16-1ZX4"
params.repRID = "Q-Y5JA"
params.outDir = "${baseDir}/../output"
params.reference = "/project/BICF/BICF_Core/shared/gudmap/references"

// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "deriva cookie file for bdbag not found: ${params.bdbag}" }
repRID = params.repRID
outDir = params.outDir
referenceBase = file ("${params.reference}")
logsDir = "${outDir}/Logs"

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbagFetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parseMeta.py")

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
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
    path 'design.csv' into metadata

  script:
    """
    hostname >>${repRID}.parseMetadata.err
    ulimit -a >>${repRID}.parseMetadata.err

    # Check replicate RID metadata
    rep=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p repRID)
    echo "LOG: replicate RID metadata parsed: \${rep}" >>${repRID}.parseMetadata.err
    
    # Get endedness metadata
    endsMeta=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p endsMeta)
    echo "LOG: endedness metadata parsed: \${endsMeta}" >>${repRID}.parseMetadata.err
    
    # Manually get endness
    endsManual=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p endsManual)
    echo "LOG: endedness manually detected: \${endsManual}" >>${repRID}.parseMetadata.err

    # Get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p stranded)
    echo "LOG: strandedness metadata parsed: \${stranded}" >>${repRID}.parseMetadata.err
    
    # Get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p spike)
    echo "LOG: spike-in metadata parsed: \${spike}" >>${repRID}.parseMetadata.err
    
    # Get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentMeta}" -p species)
    echo "LOG: species metadata parsed: \${species}" >>${repRID}.parseMetadata.err

    #Set the correct Reference Directory
    if [ "\${spike}" == 'yes' ]; then
      if [ "\${species}" == 'Homo sapiens' ]; then
        referenceDir="/project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12-S/hisat2/genome";
        echo "LOG: Referece set to Homo sapiens with spike-in." >>${repRID}.parseMetadata.err;
      elif [ "\${species}" == 'Mus musculus' ]; then
        referenceDir="/project/BICF/BICF_Core/shared/gudmap/references/GRCm38.P6-S/hisat2/genome";
        echo "LOG: Referece set to Mus musculus with spike-in." >>${repRID}.parseMetadata.err;
      fi;
    else
      if [ "\${species}" == 'Homo sapiens' ]; then
        referenceDir="/project/BICF/BICF_Core/shared/gudmap/references/GRCh38.p12/hisat2/genome";
        echo "LOG: Referece set to Homo sapiens without spike-in." >>${repRID}.parseMetadata.err;
      elif [ "\${species}" == 'Mus musculus' ]; then
        referenceDir="/project/BICF/BICF_Core/shared/gudmap/references/GRCm38.P6/hisat2/genome";
        echo "LOG: Referece set to Mus musculus without spike-in." >>${repRID}.parseMetadata.err;
      fi;
    fi;
    
    # Save design file
    echo "\${rep},\${endsMeta},\${endsManual},\${stranded},\${spike},\${species},\${referenceDir}" > design.csv
    """
}

// Split metadata into separate channels
rep = Channel.create()
endsMeta = Channel.create()
endsManual = Channel.create()
stranded = Channel.create()
spike = Channel.create()
species = Channel.create()
referenceDir = Channel.create()
metadata.splitCsv(sep: ',', header: false).separate(
  rep,
  endsMeta,
  endsManual,
  stranded,
  spike,
  species,
  referenceDir
)
endsManual.into {
  endsManual_trimData
  endsManual_alignReads
}
stranded.into {
  stranded_alignReads
}
referenceDir.into {
  referenceDir_alignReads
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.trimData.*"

  input:
    val endsManual_trimData
    file(fastq) from fastqs

  output:
    path ("*.fq.gz") into fastqs_trimmed
    file ("${repRID}.trimData.log")
    file ("${repRID}.trimData.err")

  script:
    """
    hostname >>${repRID}.trimData.err
    ulimit -a >>${repRID}.trimData.err

    # trim fastqs
    if [ '${endsManual_trimData}' == 'se' ]
    then
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    else
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err;
    fi
    """
}

/*
 * alignReads: aligns the reads to a reference database
*/
process alignReads {
  tag "align-${repRID}"
  publishDir "${outDir}/aligned", mode: "copy", pattern: "${repRID}.{unal,sorted}.{gz,bam,bai}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.align.{out,err}"

  input:
    val endsManual_alignReads
    val stranded_alignReads
    path fqs from fastqs_trimmed
    val referenceDir_alignReads

  output:
    set repRID, file ("${repRID}.unal.gz"), file ("${repRID}.sorted.bam"), file ("${repRID}.sorted.bai")
    set repRID, file ("${repRID}.align.out"), file ("${repRID}.align.err")

  script:
    """
    if [ "${endsManual_alignReads}" == 'pe' ]; then
    hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x ${referenceDir_alignReads} ${stranded_alignReads} --no-mixed --no-discordant -1 ${fqs[0]} -2 ${fqs[1]} 1>${repRID}.align.out 2>${repRID}.align.err;
    else hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x ${referenceDir_alignReads} ${stranded_alignReads} -U ${fqs[0]} 1>${repRID}.align.out 2>${repRID}.align.err;
    fi;
    samtools view -1 --threads `nproc` -o ${repRID}.bam ${repRID}.sam 1>>${repRID}.align.out 2>>${repRID}.align.err;
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.bam ${repRID}.bam 1>>${repRID}.align.out 2>>${repRID}.align.err;
    samtools index -@ `nproc` -b ${repRID}.sorted.bam ${repRID}.sorted.bai 1>>${repRID}.align.out 2>>${repRID}.align.err;
    """
}
