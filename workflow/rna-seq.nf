#!/usr/bin/env nextflow

//  ########  ####  ######  ########
//  ##     ##  ##  ##    ## ##
//  ##     ##  ##  ##       ##
//  ########   ##  ##       ######
//  ##     ##  ##  ##       ##
//  ##     ##  ##  ##    ## ##
//  ########  ####  ######  ##

// Define input variables
params.deriva = "${baseDir}/../test_data/auth/credential.json"
params.bdbag = "${baseDir}/../test_data/auth/cookies.txt"
//params.repRID = "16-1ZX4"
params.repRID = "Q-Y5F6"
params.source = "dev"
params.refMoVersion = "38.p6.vM25"
params.refHuVersion = "38.p13.v36"
params.refERCCVersion = "92"
params.outDir = "${baseDir}/../output"
params.upload = false
params.email = ""
params.track = false


// Define override input variable
params.refSource = "biohpc"
params.inputBagForce = ""
params.fastqsForce = ""
params.speciesForce = ""
params.strandedForce = ""
params.spikeForce = ""

// Define tracking input variables
params.ci = false
params.dev = true


// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
deriva.into {
  deriva_getBag
  deriva_getRefInfer
  deriva_getRef
  deriva_uploadInputBag
  deriva_uploadExecutionRun
  deriva_uploadQC
  deriva_uploadQC_fail
  deriva_uploadProcessedFile
  deriva_uploadOutputBag
  deriva_finalizeExecutionRun
  deriva_failPreExecutionRun
  deriva_failExecutionRun
}
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "deriva cookie file for bdbag not found: ${params.bdbag}" }
repRID = params.repRID
refMoVersion = params.refMoVersion
refHuVersion = params.refHuVersion
refERCCVersion = params.refERCCVersion
outDir = params.outDir
logsDir = "${outDir}/Logs"
upload = params.upload
inputBagForce = params.inputBagForce
fastqsForce = params.fastqsForce
speciesForce = params.speciesForce
strandedForce = params.strandedForce
spikeForce = params.spikeForce
email = params.email

// Define fixed files and variables
replicateExportConfig = Channel.fromPath("${baseDir}/conf/Replicate_For_Input_Bag.json")
executionRunExportConfig = Channel.fromPath("${baseDir}/conf/Execution_Run_For_Output_Bag.json")
if (params.source == "dev") {
  source = "dev.gudmap.org"
} else if (params.source == "staging") {
  source = "staging.gudmap.org"
} else if (params.source == "production") {
  source = "www.gudmap.org"
}
if (params.refSource == "biohpc") {
  referenceBase = "/project/BICF/BICF_Core/shared/gudmap/references/new"
} else if (params.refSource == "datahub") {
  referenceBase = "www.gudmap.org"
}
referenceInfer = Channel.fromList(["ERCC","GRCh","GRCm"])
multiqcConfig = Channel.fromPath("${baseDir}/conf/multiqc_config.yaml")
bicfLogo = Channel.fromPath("${baseDir}/../docs/bicf_logo.png")
softwareReferences = Channel.fromPath("${baseDir}/../docs/software_references_mqc.yaml")
softwareVersions = Channel.fromPath("${baseDir}/../docs/software_versions_mqc.yaml")

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbag_fetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parse_meta.py")
script_inferMeta = Channel.fromPath("${baseDir}/scripts/infer_meta.sh")
script_refDataInfer = Channel.fromPath("${baseDir}/scripts/extract_ref_data.py")
script_refData = Channel.fromPath("${baseDir}/scripts/extract_ref_data.py")
script_calculateTPM = Channel.fromPath("${baseDir}/scripts/calculateTPM.R")
script_convertGeneSymbols = Channel.fromPath("${baseDir}/scripts/convertGeneSymbols.R")
script_tinHist = Channel.fromPath("${baseDir}/scripts/tin_hist.py")
script_uploadInputBag = Channel.fromPath("${baseDir}/scripts/upload_input_bag.py")
script_uploadExecutionRun_uploadExecutionRun = Channel.fromPath("${baseDir}/scripts/upload_execution_run.py")
script_uploadExecutionRun_finalizeExecutionRun = Channel.fromPath("${baseDir}/scripts/upload_execution_run.py")
script_uploadExecutionRun_failPreExecutionRun = Channel.fromPath("${baseDir}/scripts/upload_execution_run.py")
script_uploadExecutionRun_failExecutionRun = Channel.fromPath("${baseDir}/scripts/upload_execution_run.py")
script_uploadQC = Channel.fromPath("${baseDir}/scripts/upload_qc.py")
script_uploadQC_fail = Channel.fromPath("${baseDir}/scripts/upload_qc.py")
script_uploadOutputBag = Channel.fromPath("${baseDir}/scripts/upload_output_bag.py")
script_deleteEntry_uploadQC = Channel.fromPath("${baseDir}/scripts/delete_entry.py")
script_deleteEntry_uploadQC_fail = Channel.fromPath("${baseDir}/scripts/delete_entry.py")
script_deleteEntry_uploadProcessedFile = Channel.fromPath("${baseDir}/scripts/delete_entry.py")

/*
 * trackStart: track start of pipeline
 */
process trackStart {
  container 'docker://gudmaprbk/gudmap-rbk_base:1.0.0'
  script:
    """
    hostname
    ulimit -a

    curl -H 'Content-Type: application/json' -X PUT -d \
      '{ \
        "sessionId": "${workflow.sessionId}", \
        "pipeline": "gudmap.rbk_rnaseq", \
        "start": "${workflow.start}", \
        "repRID": "${repRID}", \
        "astrocyte": false, \
        "status": "started", \
        "nextflowVersion": "${workflow.nextflow.version}", \
        "pipelineVersion": "${workflow.manifest.version}", \
        "ci": ${params.ci}, \
        "dev": ${params.dev} \
      }' \
      "https://xku43pcwnf.execute-api.us-east-1.amazonaws.com/ProdDeploy/pipeline-tracking"

    if [ ${params.track} == true ]
    then
      curl -H 'Content-Type: application/json' -X PUT -d \
        '{ \
          "ID": "${workflow.sessionId}", \
          "repRID": "${repRID}", \
          "PipelineVersion": "${workflow.manifest.version}", \
          "Server": "${params.source}", \
          "Queued": "NA", \
          "CheckedOut": "NA", \
          "Started": "${workflow.start}" \
        }' \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/track"
    fi
    """
}

log.info """\
====================================
BICF RNA-seq Pipeline for GUDMAP/RBK
====================================
Replicate RID          : ${params.repRID}
Source Server          : ${params.source}
Mouse Reference Version: ${params.refMoVersion}
Human Reference Version: ${params.refHuVersion}
ERCC Reference Version : ${params.refERCCVersion}
Reference source       : ${params.refSource}
Output Directory       : ${params.outDir}
Upload                 : ${upload}
Track                  : ${params.track}
------------------------------------
Nextflow Version       : ${workflow.nextflow.version}
Pipeline Version       : ${workflow.manifest.version}
Session ID             : ${workflow.sessionId}
------------------------------------
CI                     : ${params.ci}
Development            : ${params.dev}
------------------------------------
"""

/*
 * getBag: download input bag
 */
process getBag {
  tag "${repRID}"
  publishDir "${outDir}/inputBag", mode: 'copy', pattern: "*_inputBag_*.zip"

  input:
    path credential, stageAs: "credential.json" from deriva_getBag
    path replicateExportConfig

  output:
    path ("*.zip") into bag

  when:
    inputBagForce == ""

  script:
    """
    hostname > ${repRID}.getBag.log
    ulimit -a >> ${repRID}.getBag.log

    # link credential file for authentication
    echo -e "LOG: linking deriva credentials" >> ${repRID}.getBag.log
    mkdir -p ~/.deriva
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json
    echo -e "LOG: linked" >> ${repRID}.getBag.log

    # deriva-download replicate RID
    echo -e "LOG: fetching bag for ${repRID} in GUDMAP" >> ${repRID}.getBag.log
    deriva-download-cli ${source} --catalog 2 ${replicateExportConfig} . rid=${repRID}
    echo -e "LOG: fetched" >> ${repRID}.getBag.log

    name=\$(ls *.zip)
    name=\$(basename \${name} | cut -d "." -f1)
    yr=\$(date +'%Y')
    mn=\$(date +'%m')
    dy=\$(date +'%d')
    mv \${name}.zip \${name}_\${yr}\${mn}\${dy}.zip
    """
}

// Set inputBag to downloaded or forced input
if (inputBagForce != "") {
  inputBag = Channel
    .fromPath(inputBagForce)
    .ifEmpty { exit 1, "override inputBag file not found: ${inputBagForce}" }
} else {
  inputBag = bag
}
inputBag.into {
  inputBag_getData
  inputBag_uploadInputBag
}

/*
 * getData: fetch replicate files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID}"

  input:
    path script_bdbagFetch
    path cookies, stageAs: "deriva-cookies.txt" from bdbag
    path inputBag from inputBag_getData

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    path ("**/File.csv") into fileMeta
    path ("**/Experiment Settings.csv") into experimentSettingsMeta
    path ("**/Experiment.csv") into experimentMeta
    path "fastqCount.csv" into fastqCount_fl

  script:
    """
    hostname > ${repRID}.getData.log
    ulimit -a >> ${repRID}.getData.log

    # link deriva cookie for authentication
    echo -e "LOG: linking deriva cookie" >> ${repRID}.getData.log
    mkdir -p ~/.bdbag
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt
    echo -e "LOG: linked" >> ${repRID}.getData.log

    # get bag basename
    replicate=\$(basename "${inputBag}")
    echo -e "LOG: bag replicate name \${replicate}" >> ${repRID}.getData.log

    # unzip bag
    echo -e "LOG: unzipping replicate bag" >> ${repRID}.getData.log
    unzip ${inputBag}
    echo -e "LOG: unzipped" >> ${repRID}.getData.log

    # bag fetch fastq's only and rename by repRID
    echo -e "LOG: fetching replicate bdbag" >> ${repRID}.getData.log
    sh ${script_bdbagFetch} \${replicate::-13} ${repRID}
    echo -e "LOG: fetched" >> ${repRID}.getData.log
    
    fastqCount=\$(ls *.fastq.gz | wc -l)
    if [ "\${fastqCount}" == "0" ]
    then
      touch dummy.R1.fastq.gz
    fi
    echo "\${fastqCount}" > fastqCount.csv
    """
}

// Split fastq count into channel
fastqCount = Channel.create()
fastqCount_fl.splitCsv(sep: ",", header: false).separate(
  fastqCount
)

// Set raw fastq to downloaded or forced input and replicate them for multiple process inputs
if (fastqsForce != "") {
  Channel
    .fromPath(fastqsForce)
    .ifEmpty { exit 1, "override inputBag file not found: ${fastqsForce}" }
    .collect().into {
      fastqs_parseMetadata
      fastqs_fastqc
    }
} else {
  fastqs.collect().into {
    fastqs_parseMetadata
    fastqs_fastqc
  }
}

/*
 * parseMetadata: parses metadata to extract experiment parameters
*/
process parseMetadata {
  tag "${repRID}"

  input:
    path script_parseMeta
    path file from fileMeta
    path experimentSettings, stageAs: "ExperimentSettings.csv" from experimentSettingsMeta
    path experiment from experimentMeta
    path (fastq) from fastqs_parseMetadata.collect()
    val fastqCount

  output:
    path "design.csv" into metadata_fl
    path "fastqError.csv" into fastqError_fl

  script:
    """
    hostname > ${repRID}.parseMetadata.log
    ulimit -a >> ${repRID}.parseMetadata.log

    # check replicate RID metadata
    rep=\$(python3 ${script_parseMeta} -r ${repRID} -m "${file}" -p repRID)
    echo -e "LOG: replicate RID metadata parsed: \${rep}" >> ${repRID}.parseMetadata.log

    # get experiment RID metadata
    exp=\$(python3 ${script_parseMeta} -r ${repRID} -m "${file}" -p expRID)
    echo -e "LOG: experiment RID metadata parsed: \${exp}" >> ${repRID}.parseMetadata.log

    # get study RID metadata
    study=\$(python3 ${script_parseMeta} -r ${repRID} -m "${file}" -p studyRID)
    echo -e "LOG: study RID metadata parsed: \${study}" >> ${repRID}.parseMetadata.log

    # get endedness metadata
    endsRaw=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p endsMeta)
    echo -e "LOG: endedness metadata parsed: \${endsRaw}" >> ${repRID}.parseMetadata.log
    if [ "\${endsRaw}" == "Single End" ]
    then
      endsMeta="se"
    elif [ "\${endsRaw}" == "Paired End" ]
    then
      endsMeta="pe"
    elif [ "\${endsRaw}" == "Single Read" ]
    # "Single Read" depreciated as of Jan 2021, this option is present for backwards compatibility
    then
      endsMeta="se"
    elif [ "\${endsRaw}" == "nan" ]
    then
      endsRaw="_No value_"
      endsMeta="NA"
    fi

    # ganually get endness
    if [ "${fastqCount}" == "1" ]
    then
      endsManual="se"
    else
      endsManual="pe"
    fi
    echo -e "LOG: endedness manually detected: ${fastqCount}" >> ${repRID}.parseMetadata.log

    # get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p stranded)
    echo -e "LOG: strandedness metadata parsed: \${stranded}" >> ${repRID}.parseMetadata.log
    if [ "\${stranded}" == "nan" ]
    then
      stranded="_No value_"
    fi

    # get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p spike)
    echo -e "LOG: spike-in metadata parsed: \${spike}" >> ${repRID}.parseMetadata.log
    if [ "\${spike}" == "nan" ]
    then
      spike="_No value_"
    fi
    if [ "\${spike}" == "f" ]
    then
      spike="false"
    elif [ "\${spike}" == "t" ]
    then
      spike="true"
    elif [ "\${spike}" == "no" ]
    # "yes"/"no" depreciated as of Jan 2021, this option is present for backwards compatibility
    then
      spike="false"
    elif [ "\${spike}" == "yes" ]
    # "yes"/"no" depreciated as of Jan 2021, this option is present for backwards compatibility
    then
      spike="true"
    elif [ "\${spike}" == "nan" ]
    then
      endsRaw="_No value_"
      endsMeta="NA"
    fi

    # get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experiment}" -p species)
    echo -e "LOG: species metadata parsed: \${species}" >> ${repRID}.parseMetadata.log
    if [ "\${species}" == "nan" ]
    then
      species="_No value_"
    fi

    # get read length metadata
    readLength=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p readLength)
    if [ "\${readLength}" = "nan" ]
    then
      readLength="NA"
    fi
    echo -e "LOG: read length metadata parsed: \${readLength}" >> ${repRID}.parseMetadata.log

    # check not incorrect number of fastqs
    fastqCountError=false
    fastqCountError_details=""
    if [ "${fastqCount}" -gt "2" ]
    then
      fastqCountError=true
      fastqCountError_details="**Too many fastqs detected (>2)**"
    elif [ "${fastqCount}" -eq "0" ]
    then
      fastqCountError=true
      fastqCountError_details="**No valid fastqs detected \\(may not match {_.}R{12}.fastq.gz convention\\)**"
    elif [ "\${endsMeta}" == "se" ] && [ "${fastqCount}" -ne "1" ]
    then
      fastqCountError=true
      fastqCountError_details="**Number of fastqs detected does not match submitted endness**"
    elif [ "\${endsMeta}" == "pe" ] && [ "${fastqCount}" -ne "2" ]
    then
      fastqCountError=true
      fastqCountError_details="**Number of fastqs detected does not match submitted endness**"
    fi

    # check read counts match for fastqs
    fastqReadError=false
    fastqReadError_details=""
    if [ "\${endsManual}" == "pe" ]
    then
      r1Count=\$(zcat ${fastq[0]} | wc -l)
      r2Count=\$(zcat ${fastq[1]} | wc -l)
      if [ "\${r1Count}" -ne "\${r2Count}" ]
      then
        fastqReadError=true
        fastqReadError_details="**Number of reads do not match for R1 and R2:** there may be a trunkation or mismatch of fastq files"
      fi
    fi

    # save design file
    echo "\${endsMeta},\${endsRaw},\${endsManual},\${stranded},\${spike},\${species},\${readLength},\${exp},\${study}" > design.csv

    # save fastq error file
    echo "\${fastqCountError},\${fastqCountError_details},\${fastqReadError},\${fastqReadError_details}" > fastqError.csv
    """
}

// Split metadata into separate channels
endsMeta = Channel.create()
endsRaw = Channel.create()
endsManual = Channel.create()
strandedMeta = Channel.create()
spikeMeta = Channel.create()
speciesMeta = Channel.create()
readLengthMeta = Channel.create()
expRID = Channel.create()
studyRID = Channel.create()
metadata_fl.splitCsv(sep: ",", header: false).separate(
  endsMeta,
  endsRaw,
  endsManual,
  strandedMeta,
  spikeMeta,
  speciesMeta,
  readLengthMeta,
  expRID,
  studyRID
)

// Replicate metadata for multiple process inputs
endsMeta.into {
  endsMeta_checkMetadata
  endsMeta_aggrQC
  endsMeta_failExecutionRun
}
endsManual.into {
  endsManual_trimData
  endsManual_downsampleData
  endsManual_alignSampleData
  endsManual_aggrQC
}
strandedMeta.into {
  strandedMeta_checkMetadata
  strandedMeta_aggrQC
  strandedMeta_failExecutionRun
}
spikeMeta.into {
  spikeMeta_checkMetadata
  spikeMeta_aggrQC
  spikeMeta_failPreExecutionRun
  spikeMeta_failExecutionRun
}
speciesMeta.into {
  speciesMeta_checkMetadata
  speciesMeta_aggrQC
  speciesMeta_failPreExecutionRun
  speciesMeta_failExecutionRun
}
studyRID.into {
  studyRID_aggrQC
  studyRID_uploadInputBag
  studyRID_uploadProcessedFile
  studyRID_uploadOutputBag
}
expRID.into {
  expRID_aggrQC
  expRID_uploadProcessedFile
}

// Split fastq count error into separate channel
fastqCountError = Channel.create()
fastqCountError_details = Channel.create()
fastqReadError = Channel.create()
fastqReadError_details = Channel.create()
fastqError_fl.splitCsv(sep: ",", header: false).separate(
  fastqCountError,
  fastqCountError_details,
  fastqReadError,
  fastqReadError_details
)

//  Replicate errors for multiple process inputs
fastqCountError.into {
  fastqCountError_fastqc
  fastqCountError_trimData
  fastqCountError_getRefInfer
  fastqCountError_downsampleData
  fastqCountError_alignSampleData
  fastqCountError_inferMetadata
  fastqCountError_checkMetadata
  fastqCountError_uploadExecutionRun
  fastqCountError_getRef
  fastqCountError_alignData
  fastqCountError_dedupData
  fastqCountError_makeBigWig
  fastqCountError_countData
  fastqCountError_dataQC
  fastqCountError_aggrQC
  fastqCountError_uploadQC
  fastqCountError_uploadQC_fail
  fastqCountError_uploadProcessedFile
  fastqCountError_uploadOutputBag
  fastqCountError_failPreExecutionRun_fastq
}
fastqReadError.into {
  fastqReadError_fastqc
  fastqReadError_trimData
  fastqReadError_getRefInfer
  fastqReadError_downsampleData
  fastqReadError_alignSampleData
  fastqReadError_inferMetadata
  fastqReadError_checkMetadata
  fastqReadError_uploadExecutionRun
  fastqReadError_getRef
  fastqReadError_alignData
  fastqReadError_dedupData
  fastqReadError_makeBigWig
  fastqReadError_countData
  fastqReadError_dataQC
  fastqReadError_aggrQC
  fastqReadError_uploadQC
  fastqReadError_uploadQC_fail
  fastqReadError_uploadProcessedFile
  fastqReadError_uploadOutputBag
  fastqReadError_failPreExecutionRun_fastq
}

/*
 *fastqc: run fastqc on untrimmed fastq's
*/
process fastqc {
  tag "${repRID}"

  input:
    path (fastq) from fastqs_fastqc.collect()
    val fastqCountError_fastqc
    val fastqReadError_fastqc

  output:
    path ("*.R{1,2}.fastq.gz", includeInputs:true) into fastqs_trimData
    path ("*_fastqc.zip") into fastqc
    path ("rawReads.csv") into rawReadsInfer_fl
    path "fastqFileError.csv" into fastqFileError_fl

  when:
    fastqCountError_fastqc == 'false' && fastqReadError_fastqc == 'false'

  script:
    """
    hostname > ${repRID}.fastqc.log
    ulimit -a >> ${repRID}.fastqc.log

    # run fastqc
    echo -e "LOG: running fastq on raw fastqs" >> ${repRID}.fastqc.log
    fastqc *.fastq.gz -o . &> fastqc.out || true
    fastqcErrorOut=\$(cat fastqc.out | grep -c 'Failed to process file') || fastqcErrorOut=0
    fastqFileError=false
    fastqFileError_details=""
    if [ "\${fastqcErrorOut}" -ne "0" ]
    then
      fastqFileError=true
      fastqFileError_details="**There is an error with the structure of the fastq**"
      echo -e "LOG: There is an error with the structure of the fastq" >> ${repRID}.fastqc.log
      touch dummy_fastqc.zip
    else
      echo -e "LOG: The structure of the fastq is correct" >> ${repRID}.fastqc.log
    fi

    # count raw reads
    zcat *.R1.fastq.gz | echo \$((`wc -l`/4)) > rawReads.csv

    # save fastq error file
    echo "\${fastqFileError},\${fastqFileError_details}" > fastqFileError.csv
    """
}

// Extract number of raw reads metadata into channel
rawReadsInfer = Channel.create()
rawReadsInfer_fl.splitCsv(sep: ",", header: false).separate(
  rawReadsInfer
)

// Replicate inferred raw reads for multiple process inputs
rawReadsInfer.into {
  rawReadsInfer_aggrQC
  rawReadsInfer_uploadQC
}

// Split fastq count error into separate channel
fastqFileError = Channel.create()
fastqFileError_details = Channel.create()
fastqFileError_fl.splitCsv(sep: ",", header: false).separate(
  fastqFileError,
  fastqFileError_details
)

//  Replicate errors for multiple process inputs
fastqFileError.into {
  fastqFileError_fastqc
  fastqFileError_trimData
  fastqFileError_getRefInfer
  fastqFileError_downsampleData
  fastqFileError_alignSampleData
  fastqFileError_inferMetadata
  fastqFileError_checkMetadata
  fastqFileError_uploadExecutionRun
  fastqFileError_getRef
  fastqFileError_alignData
  fastqFileError_dedupData
  fastqFileError_makeBigWig
  fastqFileError_countData
  fastqFileError_dataQC
  fastqFileError_aggrQC
  fastqFileError_uploadQC
  fastqFileError_uploadQC_fail
  fastqFileError_uploadProcessedFile
  fastqFileError_uploadOutputBag
  fastqFileError_failPreExecutionRun_fastqFile
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"

  input:
    path (fastq) from fastqs_trimData
    val ends from endsManual_trimData
    val fastqCountError_trimData
    val fastqReadError_trimData
    val fastqFileError_trimData

  output:
    path ("*.fq.gz") into fastqsTrim
    path ("*_trimming_report.txt") into trimQC
    path ("readLength.csv") into readLengthInfer_fl

  when:
    fastqCountError_trimData == "false"
    fastqReadError_trimData == "false"
    fastqFileError_trimData == "false"

  script:
    """
    hostname > ${repRID}.trimData.log
    ulimit -a >> ${repRID}.trimData.log

    # trim fastq's using trim_galore and extract median read length
    echo -e "LOG: trimming ${ends}" >> ${repRID}.trimData.log
    if [ "${ends}" == "se" ]
    then
      trim_galore --gzip -q 25 --length 35 --basename ${repRID} ${fastq[0]}
      readLength=\$(zcat *_trimmed.fq.gz | awk '{if(NR%4==2) print length(\$1)}' | sort -n | awk '{a[NR]=\$0}END{print(NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}')
    elif [ "${ends}" == "pe" ]
    then
      trim_galore --gzip -q 25 --length 35 --paired --basename ${repRID} ${fastq[0]} ${fastq[1]}
      readLength=\$(zcat *_1.fq.gz | awk '{if(NR%4==2) print length(\$1)}' | sort -n | awk '{a[NR]=\$0}END{print(NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2}')
    fi
    echo -e "LOG: trimmed" >> ${repRID}.trimData.log
    echo -e "LOG: average trimmed read length: \${readLength}" >> ${repRID}.trimData.log

    # save read length file
    echo "\${readLength}" > readLength.csv
    """
}

// Extract calculated read length metadata into channel
readLengthInfer = Channel.create()
readLengthInfer_fl.splitCsv(sep: ",", header: false).separate(
  readLengthInfer
)

// Replicate inferred read length for multiple process inputs
readLengthInfer.into {
  readLengthInfer_aggrQC
  readLengthInfer_uploadQC
}
// Replicate trimmed fastq's for multiple process inputs
fastqsTrim.into {
  fastqsTrim_alignData
  fastqsTrim_downsampleData
}

// Combine inputs of getRefInfer
getRefInferInput = referenceInfer.combine(deriva_getRefInfer.combine(script_refDataInfer.combine(fastqCountError_getRefInfer.combine(fastqReadError_getRefInfer.combine(fastqFileError_getRefInfer)))))

/*
  * getRefInfer: dowloads appropriate reference for metadata inference
*/
process getRefInfer {
  tag "${refName}"

  input:
    tuple val (refName), path (credential, stageAs: "credential.json"), path (script_refDataInfer), val (fastqCountError), val (fastqReadError), val (fastqFileError) from getRefInferInput

  output:
    tuple val (refName), path ("hisat2", type: 'dir'), path ("*.fna"), path ("*.gtf")  into refInfer
    path ("${refName}", type: 'dir') into bedInfer

  when:
    fastqCountError == "false"
    fastqReadError == "false"
    fastqFileError == "false"

  script:
    """
    hostname > ${repRID}.${refName}.getRefInfer.log
    ulimit -a >> ${repRID}.${refName}.getRefInfer.log

    # link credential file for authentication
    echo -e "LOG: linking deriva credentials" >> ${repRID}.${refName}.getRefInfer.log
    mkdir -p ~/.deriva
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json
    echo -e "LOG: linked" >> ${repRID}.${refName}.getRefInfer.log

    # set the reference name
    if [ "${refName}" == "ERCC" ]
    then
      references=\$(echo ${referenceBase}/ERCC${refERCCVersion})
    elif [ "${refName}" == "GRCm" ]
    then
      references=\$(echo ${referenceBase}/GRCm${refMoVersion})
    elif [ '${refName}' == "GRCh" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
    else
      echo -e "LOG: ERROR - References could not be set!\nReference found: ${referenceBase}" >> ${repRID}.${refName}.getRefInfer.log
      exit 1
    fi

    # retreive appropriate reference appropriate location
    echo -e "LOG: fetching ${refName} reference files from ${referenceBase}" >> ${repRID}.${refName}.getRefInfer.log
    if [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references/new" ]
    then
      unzip \${references}.zip
      mv \$(basename \${references})/data/* .
    elif [ params.refSource == "datahub" ]
    then
      GRCv=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f1)
      GRCp=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f2)
      GENCODE=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f3)
      if [ "${refName}" != "ERCC" ]
      then
        query=\$(echo 'https://${referenceBase}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Reference_Version='\${GRCv}'.'\${GRCp}'/Annotation_Version=GENCODE%20'\${GENCODE})
      else
        query=\$(echo 'https://${referenceBase}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Reference_Version=${refName}${refERCCVersion}/Annotation_Version=${refName}${refERCCVersion}')
      fi
      curl --request GET \${query} > refQuery.json
      refURL=\$(python ${script_refDataInfer} --returnParam URL)
      loc=\$(dirname \${refURL})
      fName=\$(python ${script_refDataInfer} --returnParam fName)
      fName=\${fName%.*}
      if [ "\${loc}" = "/hatrac/*" ]; then echo "LOG: Reference not present in hatrac"; exit 1; fi
      filename=\$(echo \$(basename \${refURL}) | grep -oP '.*(?=:)')
      deriva-hatrac-cli --host ${referenceBase} get \${refURL}
      unzip \$(basename \${refURL})
      mv \${fName}/data/* .
    fi
    mv ./annotation/genome.gtf .
    mv ./sequence/genome.fna .
    mkdir ${refName}
    if [ "${refName}" != "ERCC" ]
    then
      mv ./annotation/genome.bed ./${refName}
    fi
    echo -e "LOG: fetched" >> ${repRID}.${refName}.getRefInfer.log
    """
}

/*
 * downsampleData: downsample fastq's for metadata inference
 */
process downsampleData {
  tag "${repRID}"

  input:
    path fastq from fastqsTrim_downsampleData
    val ends from endsManual_downsampleData
    val fastqCountError_downsampleData
    val fastqReadError_downsampleData
    val fastqFileError_downsampleData

  output:
    path ("sampled.1.fq") into fastqs1Sample
    path ("sampled.2.fq") into fastqs2Sample

  when:
    fastqCountError_downsampleData == "false"
    fastqReadError_downsampleData == "false"
    fastqFileError_downsampleData == "false"

  script:
    """
    hostname > ${repRID}.downsampleData.log
    ulimit -a >> ${repRID}.downsampleData.log

    if [ "${ends}" == "se" ]
    then
      echo -e "LOG: downsampling SE trimmed fastq" >> ${repRID}.downsampleData.log
      seqtk sample -s100 *trimmed.fq.gz 100000 1> sampled.1.fq
      touch sampled.2.fq
    elif [ "${ends}" == "pe" ]
    then
      echo -e "LOG: downsampling R1 of PE trimmed fastq" >> ${repRID}.downsampleData.log
      seqtk sample -s100 *1.fq.gz 1000000 1> sampled.1.fq
      echo -e "LOG: downsampling R2 of PE trimmed fastq" >> ${repRID}.downsampleData.log
      seqtk sample -s100 *2.fq.gz 1000000 1> sampled.2.fq
    fi
    echo -e "LOG: downsampled" >> ${repRID}.downsampleData.log
    """
}

// Replicate the dowsampled fastq's and attatched to the references
inferInput = endsManual_alignSampleData.combine(refInfer.combine(fastqs1Sample.collect().combine(fastqs2Sample.collect().combine(fastqCountError_alignSampleData.combine(fastqReadError_alignSampleData.combine(fastqFileError_alignSampleData))))))

/*
 * alignSampleData: aligns the downsampled reads to a reference database
*/
process alignSampleData {
  tag "${ref}"

  input:
    tuple val (ends), val (ref), path (hisat2), path (fna), path (gtf), path (fastq1), path (fastq2), val (fastqCountError), val (fastqReadError), val (fastqFileError) from inferInput

  output:
    path ("${ref}.sampled.sorted.bam") into sampleBam
    path ("${ref}.sampled.sorted.bam.bai") into sampleBai
    path ("${ref}.alignSampleSummary.txt") into alignSampleQC

  when:
    fastqCountError == "false"
    fastqReadError == "false"
    fastqFileError == "false"

  script:
    """
    hostname > ${repRID}.${ref}.alignSampleData.log
    ulimit -a >> ${repRID}.${ref}.alignSampleData.log

    # align the reads with Hisat2
    echo -e "LOG: aligning ${ends}" >> ${repRID}.${ref}.alignSampleData.log
    if [ "${ends}" == "se" ]
    then

      hisat2 -p `nproc` --add-chrname -S ${ref}.sampled.sam -x hisat2/genome -U ${fastq1} --summary-file ${ref}.alignSampleSummary.txt --new-summary
    elif [ "${ends}" == "pe" ]
    then
      hisat2 -p `nproc` --add-chrname -S ${ref}.sampled.sam -x hisat2/genome --no-mixed --no-discordant -1 ${fastq1} -2 ${fastq2} --summary-file ${ref}.alignSampleSummary.txt --new-summary
    fi
    echo -e "LOG: aliged" >> ${repRID}.${ref}.alignSampleData.log

    # convert the output sam file to a sorted bam file using Samtools
    echo -e "LOG: converting from sam to bam" >> ${repRID}.${ref}.alignSampleData.log
    samtools view -1 -@ `nproc` -F 4 -F 8 -F 256 -o ${ref}.sampled.bam ${ref}.sampled.sam

    # sort the bam file using Samtools
    echo -e "LOG: sorting the bam file" >> ${repRID}.${ref}.alignSampleData.log
    proc=\$(expr `nproc` - 1)
    mem=\$(vmstat -s -S K | grep 'total memory' | grep -o '[0-9]*')
    mem=\$(expr \${mem} / \${proc} \\* 85 / 100)
    samtools sort -@ \${proc} -m \${mem}K -O BAM -o ${ref}.sampled.sorted.bam ${ref}.sampled.bam

    # index the sorted bam using Samtools
    echo -e "LOG: indexing sorted bam file" >> ${repRID}.${ref}.alignSampleData.log
    samtools index -@ `nproc` -b ${ref}.sampled.sorted.bam ${ref}.sampled.sorted.bam.bai
    """
}

alignSampleQC.into {
  alignSampleQC_inferMetadata
  alignSampleQC_aggrQC
}

process inferMetadata {
  tag "${repRID}"

  input:
    path script_inferMeta
    path beds from bedInfer.collect()
    path bam from sampleBam.collect()
    path bai from sampleBai.collect()
    path alignSummary from alignSampleQC_inferMetadata.collect()
    val strandedForce
    val spikeForce
    val fastqCountError_inferMetadata
    val fastqReadError_inferMetadata
    val fastqFileError_inferMetadata

  output:
    path "infer.csv" into inferMetadata_fl
    path "${repRID}.infer_experiment.txt" into inferExperiment
    path "speciesError.csv" into speciesError_fl

  when:
    fastqCountError_inferMetadata == "false"
    fastqReadError_inferMetadata == "false"
    fastqFileError_inferMetadata == "false"

  script:
    """
    hostname > ${repRID}.inferMetadata.log
    ulimit -a >> ${repRID}.inferMetadata.log

    # collect alignment rates (round down to integers)
    align_ercc=\$(echo \$(grep "Overall alignment rate" ERCC.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_ercc=\$(echo \${align_ercc%.*})
    echo -e "LOG: alignment rate to ERCC: \${align_ercc}" >> ${repRID}.inferMetadata.log
    align_hu=\$(echo \$(grep "Overall alignment rate" GRCh.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_hu=\$(echo \${align_hu%.*})
    echo -e "LOG: alignment rate to GRCh: \${align_hu}" >> ${repRID}.inferMetadata.log
    align_mo=\$(echo \$(grep "Overall alignment rate" GRCm.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_mo=\$(echo \${align_mo%.*})
    echo -e "LOG: alignment rate to GRCm: \${align_mo}" >> ${repRID}.inferMetadata.log

    # determine spike-in
    if [ 1 -eq \$(echo \$(expr \${align_ercc} ">=" 10)) ]
    then
      spike="true"
    else
      spike="false"
    fi
    echo -e "LOG: inference of strandedness results is: \${spike}" >> ${repRID}.inferMetadata.log
    if [ "${spikeForce}" != "" ]
    then
      spike=${spikeForce}
      echo -e "LOG: spike-in metadata forced: \${spike}" >> ${repRID}.parseMetadata.log
    fi

    speciesError=false
    speciesError_details=""
    # determine species
    if [ 1 -eq \$(echo \$(expr \${align_hu} ">=" 40)) ] && [ 1 -eq \$(echo \$(expr \${align_mo} "<" 40)) ]
    then
      species="Homo sapiens"
      bam="GRCh.sampled.sorted.bam"
      bed="./GRCh/genome.bed"
      echo -e "LOG: inference of species results in: \${species}" >> ${repRID}.inferMetadata.log
    elif [ 1 -eq \$(echo \$(expr \${align_mo} ">=" 40)) ] && [ 1 -eq \$(echo \$(expr \${align_hu} "<" 40)) ]
    then
      species="Mus musculus"
      bam="GRCm.sampled.sorted.bam"
      bed="./GRCm/genome.bed"
      echo -e "LOG: inference of species results in: \${species}" >> ${repRID}.inferMetadata.log
    else
      echo -e "LOG: ERROR - inference of species returns an ambiguous result: hu=\${align_hu} mo=\${align_mo}" >> ${repRID}.inferMetadata.log
      if [ "${speciesForce}" == "" ]
      then
        speciesError=true
        speciesError_details="**Inference of species returns an ambiguous result:** Percent aligned to human = \${align_hu} and percent aligned to mouse = \${align_mo}"
      fi
    fi
    if [ "${speciesForce}" != "" ]
    then
      speciesError=false
      echo -e "LOG: species overridden to: ${speciesForce}"
      species="${speciesForce}"
      if [ "${speciesForce}" == "Homo sapiens" ]
      then
        bam="GRCh.sampled.sorted.bam"
        bed="./GRCh/genome.bed"
      elif [ "${speciesForce}" == "Mus musculus" ]
      then
        bam="GRCm.sampled.sorted.bam"
        bed="./GRCm/genome.bed"
      fi
    fi

    if [ "\${speciesError}" == false ]
    then
      # infer experimental setting from dedup bam
      echo -e "LOG: infer experimental setting from dedup bam" >> ${repRID}.inferMetadata.log
      infer_experiment.py -r "\${bed}" -i "\${bam}" 1>> ${repRID}.infer_experiment.txt
      echo -e "LOG: inferred" >> ${repRID}.inferMetadata.log

      ended=`bash ${script_inferMeta} endness ${repRID}.infer_experiment.txt`
      fail=`bash ${script_inferMeta} fail ${repRID}.infer_experiment.txt`
      if [ \${ended} == "PairEnd" ]
      then
        ends="pe"
        percentF=`bash ${script_inferMeta} pef ${repRID}.infer_experiment.txt`
        percentR=`bash ${script_inferMeta} per ${repRID}.infer_experiment.txt`
      elif [ \${ended} == "SingleEnd" ]
      then
        ends="se"
        percentF=`bash ${script_inferMeta} sef ${repRID}.infer_experiment.txt`
        percentR=`bash ${script_inferMeta} ser ${repRID}.infer_experiment.txt`
      fi
      echo -e "LOG: percentage reads in the same direction as gene: \${percentF}" >> ${repRID}.inferMetadata.log
      echo -e "LOG: percentage reads in the opposite direction as gene: \${percentR}" >> ${repRID}.inferMetadata.log
      if [ 1 -eq \$(echo \$(expr \${percentF#*.} ">" 2500)) ] && [ 1 -eq \$(echo \$(expr \${percentR#*.} "<" 2500)) ]
      then
        stranded="forward"
      elif [ 1 -eq \$(echo \$(expr \${percentR#*.} ">" 2500)) ] && [ 1 -eq \$(echo \$(expr \${percentF#*.} "<" 2500)) ]
      then
        stranded="reverse"
      else
        stranded="unstranded"
      fi
      echo -e "LOG: stradedness set to: \${stranded}" >> ${repRID}.inferMetadata.log
      if [ "${strandedForce}" != "" ]
      then
        stranded=${strandedForce}
        echo -e "LOG: spike-in metadata forced: \${stranded}" >> ${repRID}.inferMetadata.log
      fi
    else
      ends=""
      stranded=""
      spike=""
      species=""
      percentF=""
      percentR=""
      fail=""
      touch ${repRID}.infer_experiment.txt
    fi

    # write inferred metadata to file
    echo "\${ends},\${stranded},\${spike},\${species},\${align_ercc},\${align_hu},\${align_mo},\${percentF},\${percentR},\${fail}" > infer.csv

    # save species error file
    echo "\${speciesError},\${speciesError_details}" > speciesError.csv
    """
}

// Split metadata into separate channels
endsInfer = Channel.create()
strandedInfer = Channel.create()
spikeInfer = Channel.create()
speciesInfer = Channel.create()
align_erccInfer = Channel.create()
align_huInfer = Channel.create()
align_moInfer = Channel.create()
percentFInfer = Channel.create()
percentRInfer = Channel.create()
failInfer = Channel.create()
inferMetadata_fl.splitCsv(sep: ",", header: false).separate(
  endsInfer,
  strandedInfer,
  spikeInfer,
  speciesInfer,
  align_erccInfer,
  align_huInfer,
  align_moInfer,
  percentFInfer,
  percentRInfer,
  failInfer
)

// Replicate metadata for multiple process inputs
endsInfer.into {
  endsInfer_checkMetadata
  endsInfer_alignData
  endsInfer_countData
  endsInfer_dataQC
  endsInfer_aggrQC
  endsInfer_uploadQC
  endsInfer_failExecutionRun
}
strandedInfer.into {
  strandedInfer_checkMetadata
  strandedInfer_alignData
  strandedInfer_countData
  strandedInfer_aggrQC
  strandedInfer_uploadQC
  strandedInfer_failExecutionRun
}
spikeInfer.into{
  spikeInfer_checkMetadata
  spikeInfer_getRef
  spikeInfer_aggrQC
  spikeInfer_uploadExecutionRun
  spikeInfer_failExecutionRun
}
speciesInfer.into {
  speciesInfer_checkMetadata
  speciesInfer_getRef
  speciesInfer_aggrQC
  speciesInfer_uploadExecutionRun
  speciesInfer_uploadProcessedFile
  speciesInfer_failExecutionRun
}

// Split species count error into separate channel
speciesError = Channel.create()
speciesError_details = Channel.create()
speciesError_fl.splitCsv(sep: ",", header: false).separate(
  speciesError,
  speciesError_details
)

//  Replicate errors for multiple process inputs
speciesError.into {
  speciesError_checkMetadata
  speciesError_uploadExecutionRun
  speciesError_getRef
  speciesError_alignData
  speciesError_dedupData
  speciesError_makeBigWig
  speciesError_countData
  speciesError_fastqc
  speciesError_dataQC
  speciesError_aggrQC
  speciesError_uploadQC
  speciesError_uploadQC_fail
  speciesError_uploadProcessedFile
  speciesError_uploadOutputBag
  speciesError_failPreExecutionRun_species
}

/* 
 * checkMetadata: checks the submitted metada against inferred
*/
process checkMetadata {
  tag "${repRID}"

  input:
    val endsMeta from endsMeta_checkMetadata
    val strandedMeta from strandedMeta_checkMetadata
    val spikeMeta from spikeMeta_checkMetadata
    val speciesMeta from speciesMeta_checkMetadata
    val endsInfer from endsInfer_checkMetadata
    val strandedInfer from strandedInfer_checkMetadata
    val spikeInfer from spikeInfer_checkMetadata
    val speciesInfer from speciesInfer_checkMetadata
    val fastqCountError_checkMetadata
    val fastqReadError_checkMetadata
    val fastqFileError_checkMetadata
    val speciesError_checkMetadata

  output:
    path ("check.csv") into checkMetadata_fl
    path ("outputBagRID.csv") optional true into outputBagRID_fl_dummy

  when:
    fastqCountError_checkMetadata == "false"
    fastqReadError_checkMetadata == "false"
    fastqFileError_checkMetadata == "false"
    speciesError_checkMetadata == "false"

  script:
    """
    hostname > ${repRID}.checkMetadata.log
    ulimit -a >> ${repRID}.checkMetadata.log

    pipelineError=false
    pipelineError_ends=false
    pipelineError_stranded=false
    pipelineError_spike=false
    pipelineError_species=false
    # check if submitted metadata matches inferred
    if [ "${strandedMeta}" != "${strandedInfer}" ]
    then
      if [ "${params.strandedForce}" != "" ]
      then
        pipelineError=false
        pipelineError_stranded=false
        echo -e "LOG: stranded forced: Submitted=${strandedMeta}; Inferred=${strandedInfer}" >> ${repRID}.checkMetadata.log
      else
        pipelineError=true
        pipelineError_stranded=true
        if [ "${strandedMeta}" == "stranded" ]
        then
          if [[ "${strandedInfer}" == "forward" ]] || [[ "${strandedInfer}" == "reverse" ]]
          then
            pipelineError=false
            pipelineError_stranded=false
            echo -e "LOG: stranded matches: Submitted=${strandedMeta}; Inferred=${strandedInfer}" >> ${repRID}.checkMetadata.log
          else
            echo -e "LOG: stranded does not match: Submitted=${strandedMeta}; Inferred=${strandedInfer}" >> ${repRID}.checkMetadata.log
          fi
        else
          echo -e "LOG: stranded does not match: Submitted=${strandedMeta}; Inferred=${strandedInfer}" >> ${repRID}.checkMetadata.log
        fi
      fi
    else
      pipelineError=false
      pipelineError_stranded=false
      echo -e "LOG: stranded matches: Submitted=${strandedMeta}; Inferred=${strandedInfer}" >> ${repRID}.checkMetadata.log
    fi
    if [ "${endsMeta}" != "${endsInfer}" ]
    then
      pipelineError=true
      pipelineError_ends=true
      echo -e "LOG: ends do not match: Submitted=${endsMeta}; Inferred=${endsInfer}" >> ${repRID}.checkMetadata.log
    else
      pipelineError_ends=false
      echo -e "LOG: ends matches: Submitted=${endsMeta}; Inferred=${endsInfer}" >> ${repRID}.checkMetadata.log
    fi
    if [ "${spikeMeta}" != "${spikeInfer}" ]
    then
      if [[ "${params.spikeForce}" != "" ]]
      then
        pipelineError_spike=false
        echo -e "LOG: spike forced: Submitted=${spikeMeta}; Inferred=${spikeInfer}" >> ${repRID}.checkMetadata.log
      else
        pipelineError=true
        pipelineError_spike=true
        echo -e "LOG: spike does not match: Submitted=${spikeMeta}; Inferred=${spikeInfer}" >> ${repRID}.checkMetadata.log
      fi
    else
      pipelineError_spike=false
      echo -e "LOG: spike matches: Submitted=${spikeMeta}; Inferred=${spikeInfer}" >> ${repRID}.checkMetadata.log
    fi
    if [ "${speciesMeta}" != "${speciesInfer}" ]
    then
    if [[ "${params.speciesForce}" != "" ]]
      then
        pipelineError_species=false
        echo -e "LOG: species forced: Submitted=${speciesMeta}; Inferred=${speciesInfer}" >> ${repRID}.checkMetadata.log
      else
        pipelineError=true
        pipelineError_species=true
        echo -e "LOG: species does not match: Submitted=${speciesMeta}; Inferred=${speciesInfer}" >> ${repRID}.checkMetadata.log
      fi
    else
      pipelineError_species=false
      echo -e "LOG: species matches: Submitted=${speciesMeta}; Inferred=${speciesInfer}" >> ${repRID}.checkMetadata.log
    fi

    # create dummy output bag rid if failure
    if [ \${pipelineError} == true ]
    then
      echo "fail" > outputBagRID.csv
    fi

    # write checks to file
    echo "\${pipelineError},\${pipelineError_ends},\${pipelineError_stranded},\${pipelineError_spike},\${pipelineError_species}" > check.csv
    """
}

// Split errors into separate channels
pipelineError = Channel.create()
pipelineError_ends = Channel.create()
pipelineError_stranded = Channel.create()
pipelineError_spike = Channel.create()
pipelineError_species = Channel.create()
checkMetadata_fl.splitCsv(sep: ",", header: false).separate(
  pipelineError,
  pipelineError_ends,
  pipelineError_stranded,
  pipelineError_spike,
  pipelineError_species
)

// Replicate errors for multiple process inputs
pipelineError.into {
  pipelineError_getRef
  pipelineError_alignData
  pipelineError_dedupData
  pipelineError_makeBigWig
  pipelineError_countData
  pipelineError_fastqc
  pipelineError_dataQC
  pipelineError_aggrQC
  pipelineError_uploadQC
  pipelineError_uploadQC_fail
  pipelineError_uploadProcessedFile
  pipelineError_uploadOutputBag
  pipelineError_failExecutionRun
}

/* 
 * uploadInputBag: uploads the input bag
*/
process uploadInputBag {
  tag "${repRID}"

  input:
    path script_uploadInputBag
    path credential, stageAs: "credential.json" from deriva_uploadInputBag
    path inputBag from inputBag_uploadInputBag
    val studyRID from studyRID_uploadInputBag

  output:
    path ("inputBagRID.csv") into inputBagRID_fl

  when:
    upload

  script:
    """
    hostname > ${repRID}.uploadInputBag.log
    ulimit -a >> ${repRID}.uploadInputBag.log

    yr=\$(date +'%Y')
    mn=\$(date +'%m')
    dy=\$(date +'%d')

    file=\$(basename -a ${inputBag})
    md5=\$(md5sum ./\${file} | awk '{ print \$1 }')
    echo LOG: ${repRID} input bag md5 sum - \${md5} >> ${repRID}.uploadInputBag.log
    size=\$(wc -c < ./\${file})
    echo LOG: ${repRID} input bag size - \${size} bytes >> ${repRID}.uploadInputBag.log
    
    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Input_Bag/File_MD5=\${md5})
    if [ "\${exist}" == "[]" ]
    then
        cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
        cookie=\${cookie:11:-1}

        loc=\$(deriva-hatrac-cli --host ${source} put ./\${file} /hatrac/resources/rnaseq/pipeline/input_bag/study/${studyRID}/replicate/${repRID}/\${file} --parents)
        inputBag_rid=\$(python3 ${script_uploadInputBag} -f \${file} -l \${loc} -s \${md5} -b \${size} -o ${source} -c \${cookie})
        echo LOG: input bag RID uploaded - \${inputBag_rid} >> ${repRID}.uploadInputBag.log
        rid=\${inputBag_rid}
    else
        exist=\$(echo \${exist} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
        exist=\${exist:7:-6}
        echo LOG: input bag RID already exists - \${exist} >> ${repRID}.uploadInputBag.log
        rid=\${exist}
    fi

    echo "\${rid}" > inputBagRID.csv
    """
}

// Extract input bag RID into channel
inputBagRID = Channel.create()
inputBagRID_fl.splitCsv(sep: ",", header: false).separate(
  inputBagRID
)

// Replicate input bag RID for multiple process inputs
inputBagRID.into {
  inputBagRID_uploadExecutionRun
  inputBagRID_finalizeExecutionRun
  inputBagRID_failPreExecutionRun
  inputBagRID_failExecutionRun
}

/* 
 * uploadExecutionRun: uploads the execution run
*/
process uploadExecutionRun {
  tag "${repRID}"

  input:
    path script_uploadExecutionRun_uploadExecutionRun
    path credential, stageAs: "credential.json" from deriva_uploadExecutionRun
    val spike from spikeInfer_uploadExecutionRun
    val species from speciesInfer_uploadExecutionRun
    val inputBagRID from inputBagRID_uploadExecutionRun
    val fastqCountError_uploadExecutionRun
    val fastqReadError_uploadExecutionRun
    val fastqFileError_uploadExecutionRun
    val speciesError_uploadExecutionRun
    
  output:
    path ("executionRunRID.csv") into executionRunRID_fl

  when:
    upload
    fastqCountError_uploadExecutionRun == "false"
    fastqReadError_uploadExecutionRun == "false"
    fastqFileError_uploadExecutionRun == "false"
    speciesError_uploadExecutionRun == "false"

  script:
    """
    hostname > ${repRID}.uploadExecutionRun.log
    ulimit -a >> ${repRID}.uploadExecutionRun.log

    echo LOG: searching for workflow RID - BICF mRNA ${workflow.manifest.version} >> ${repRID}.uploadExecutionRun.log
    workflow=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Workflow/Name=BICF%20mRNA%20Replicate/Version=${workflow.manifest.version})
    workflow=\$(echo \${workflow} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
    workflow=\${workflow:7:-6}
    echo LOG: workflow RID extracted - \${workflow} >> ${repRID}.uploadExecutionRun.log

    if [ "${species}" == "Homo sapiens" ]
    then
      genomeName=\$(echo GRCh${refHuVersion})
    elif [ "${species}" == "Mus musculus" ]
    then
      genomeName=\$(echo GRCm${refMoVersion})
    fi
    if [ "${spike}" == "true" ]
    then
      genomeName=\$(echo \${genomeName}-S)
    fi
    echo LOG: searching for genome name - \${genomeName} >> ${repRID}.uploadExecutionRun.log
    genome=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Name=\${genomeName})
    genome=\$(echo \${genome} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
    genome=\${genome:7:-6}
    echo LOG: genome RID extracted - \${genome} >> ${repRID}.uploadExecutionRun.log

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Execution_Run/Workflow=\${workflow}/Replicate=${repRID}/Input_Bag=${inputBagRID})
    echo \${exist} >> ${repRID}.uploadExecutionRun.log
    if [ "\${exist}" == "[]" ]
    then
      executionRun_rid=\$(python3 ${script_uploadExecutionRun_uploadExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s In-progress -d 'Run in process' -o ${source} -c \${cookie} -u F)
      echo LOG: execution run RID uploaded - \${executionRun_rid} >> ${repRID}.uploadExecutionRun.log
    else
      rid=\$(echo \${exist} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
      rid=\${rid:7:-6}
      echo \${rid} >> ${repRID}.uploadExecutionRun.log
      executionRun_rid=\$(python3 ${script_uploadExecutionRun_uploadExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s In-progress -d 'Run in process' -o ${source} -c \${cookie} -u \${rid})
      echo LOG: execution run RID updated - \${executionRun_rid} >> ${repRID}.uploadExecutionRun.log
    fi

    echo "\${executionRun_rid}" > executionRunRID.csv

    if [ ${params.track} == true ]
    then
      curl -H 'Content-Type: application/json' -X PUT -d \
        '{ \
          "ID": "${workflow.sessionId}", \
          "ExecutionRunRID": "'\${executionRun_rid}'" \
        }' \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/track"
    fi
    """
}

// Extract execution run RID into channel
executionRunRID = Channel.create()
executionRunRID_fl.splitCsv(sep: ",", header: false).separate(
  executionRunRID
)

// Replicate execution run RID for multiple process inputs
executionRunRID.into {
  executionRunRID_uploadQC
  executionRunRID_uploadProcessedFile
  executionRunRID_uploadOutputBag
  executionRunRID_finalizeExecutionRun
  executionRunRID_failExecutionRun
  executionRunRID_fail
}

/*
  * getRef: downloads appropriate reference
*/
process getRef {
  tag "${species}"

  input:
    path script_refData
    path credential, stageAs: "credential.json" from deriva_getRef
    val spike from spikeInfer_getRef
    val species from speciesInfer_getRef
    val fastqCountError_getRef
    val fastqReadError_getRef
    val fastqFileError_getRef
    val speciesError_getRef
    val pipelineError_getRef

  output:
    tuple path ("hisat2", type: 'dir'), path ("*.bed"), path ("*.fna"), path ("*.gtf"), path ("geneID.tsv"), path ("Entrez.tsv")  into reference

  when:
    fastqCountError_getRef == "false"
    fastqReadError_getRef == "false"
    fastqFileError_getRef == "false"
    speciesError_getRef == "false"
    pipelineError_getRef == "false"

  script:
    """
    hostname > ${repRID}.getRef.log
    ulimit -a >> ${repRID}.getRef.log

    # link credential file for authentication
    echo -e "LOG: linking deriva credentials" >> ${repRID}.getRef.log
    mkdir -p ~/.deriva
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json
    echo -e "LOG: linked" >> ${repRID}.getRef.log

    # set the reference name
    if [ "${species}" == "Mus musculus" ]
    then
      references=\$(echo ${referenceBase}/GRCm${refMoVersion})
      refName=GRCm
    elif [ '${species}' == "Homo sapiens" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
      refName=GRCh
    else
      echo -e "LOG: ERROR - References could not be set!\nSpecies reference found: ${species}" >> ${repRID}.getRef.log
      exit 1
    fi
    if [ "${spike}" == "true" ]
    then
      references=\$(echo \${reference}-S)
    elif [ "${spike}" == "false" ]
    then
      reference=\$(echo \${references})
    fi
    echo -e "LOG: species set to \${references}" >> ${repRID}.getRef.log

    # retreive appropriate reference appropriate location
    echo -e "LOG: fetching ${species} reference files from ${referenceBase}" >> ${repRID}.getRef.log
    if [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references/new" ]
    then
      echo -e "LOG: grabbing reference files from local (BioHPC)" >> ${repRID}.getRef.log
      unzip \${reference}.zip
      mv \$(basename \${reference})/data/* .
    elif [ arams.refSource == "datahub" ]
    then
      echo -e "LOG: grabbing reference files from datahub" >> ${repRID}.getRef.log
      GRCv=\$(echo \${references} | grep -o \${refName}.* | cut -d '.' -f1)
      GRCp=\$(echo \${references} | grep -o \${refName}.* | cut -d '.' -f2)
      GENCODE=\$(echo \${references} | grep -o \${refName}.* | cut -d '.' -f3)
      query=\$(echo 'https://${referenceBase}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Reference_Version='\${GRCv}'.'\${GRCp}'/Annotation_Version=GENCODE%20'\${GENCODE})
      curl --request GET \${query} > refQuery.json
      refURL=\$(python ${script_refData} --returnParam URL)
      loc=\$(dirname \${refURL})
      fName=\$(python ${script_refData} --returnParam fName)
      fName=\${fName%.*}
      if [ "\${loc}" = "/hatrac/*" ]; then echo "LOG: Reference not present in hatrac"; exit 1; fi
      filename=\$(echo \$(basename \${refURL}) | grep -oP '.*(?=:)')
      deriva-hatrac-cli --host ${referenceBase} get \${refURL}
      unzip \$(basename \${refURL})
      mv \${fName}/data/* .
    fi
    echo -e "LOG: fetched" >> ${repRID}.getRef.log

    mv ./annotation/genome.gtf .
    mv ./sequence/genome.fna .
    mv ./annotation/genome.bed .
    mv ./metadata/Entrez.tsv .
    mv ./metadata/geneID.tsv .
    """
}

// Replicate reference for multiple process inputs
reference.into {
  reference_alignData
  reference_countData
  reference_dataQC
}

/*
 * alignData: aligns the reads to a reference database
*/
process alignData {
  tag "${repRID}"

  input:
    path fastq from fastqsTrim_alignData
    path reference_alignData
    val ends from endsInfer_alignData
    val stranded from strandedInfer_alignData
    val fastqCountError_alignData
    val fastqReadError_alignData
    val fastqFileError_alignData
    val speciesError_alignData
    val pipelineError_alignData

  output:
    tuple path ("${repRID}.sorted.bam"), path ("${repRID}.sorted.bam.bai") into rawBam
    path ("*.alignSummary.txt") into alignQC

  when:
    fastqCountError_alignData == "false"
    fastqReadError_alignData == "false"
    fastqFileError_alignData == "false"
    speciesError_alignData == "false"
    pipelineError_alignData == "false"

  script:
    """
    hostname > ${repRID}.align.log
    ulimit -a >> ${repRID}.align.log

    # set stranded param for hisat2
    if [ "${stranded}"=="unstranded" ]
    then
      strandedParam=""
    elif [ "${stranded}" == "forward" ] && [ "${ends}" == "se" ]
    then
        strandedParam="--rna-strandness F"
    elif [ "${stranded}" == "forward" ] && [ "${ends}" == "pe" ]
    then
      strandedParam="--rna-strandness FR"
    elif [ "${stranded}" == "reverse" ] && [ "${ends}" == "se" ]
    then
        strandedParam="--rna-strandness R"
    elif [ "${stranded}" == "reverse" ] && [ "${ends}" == "pe" ]
    then
      strandedParam="--rna-strandness RF"
    fi

    # align the reads with Hisat2
    echo -e "LOG: aligning ${ends}" >> ${repRID}.align.log
    if [ "${ends}" == "se" ]
    then
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome \${strandedParam} -U ${fastq[0]} --summary-file ${repRID}.alignSummary.txt --new-summary
    elif [ "${ends}" == "pe" ]
    then
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome \${strandedParam} --no-mixed --no-discordant -1 ${fastq[0]} -2 ${fastq[1]} --summary-file ${repRID}.alignSummary.txt --new-summary
    fi
    echo -e "LOG: alignined" >> ${repRID}.align.log

    # convert the output sam file to a sorted bam file using Samtools
    echo -e "LOG: converting from sam to bam" >> ${repRID}.align.log
    samtools view -1 -@ `nproc` -F 4 -F 8 -F 256 -o ${repRID}.bam ${repRID}.sam

    # sort the bam file using Samtools
    echo -e "LOG: sorting the bam file" >> ${repRID}.align.log
    proc=\$(expr `nproc` - 1)
    mem=\$(vmstat -s -S K | grep 'total memory' | grep -o '[0-9]*')
    mem=\$(expr \${mem} / \${proc} \\* 75 / 100)
    samtools sort -@ \${proc} -m \${mem}K -O BAM -o ${repRID}.sorted.bam ${repRID}.bam

    # index the sorted bam using Samtools
    echo -e "LOG: indexing sorted bam file" >> ${repRID}.align.log
    samtools index -@ `nproc` -b ${repRID}.sorted.bam ${repRID}.sorted.bam.bai
    """
}

// Replicate rawBam for multiple process inputs
rawBam.set {
  rawBam_dedupData
}

/*
 *dedupData: mark the duplicate reads, specifically focused on PCR or optical duplicates
*/
process dedupData {
  tag "${repRID}"
  publishDir "${outDir}/bam", mode: 'copy', pattern: "*.deduped.bam"

  input:
    tuple path (bam), path (bai) from rawBam_dedupData
    val fastqCountError_dedupData
    val fastqReadError_dedupData
    val fastqFileError_dedupData
    val speciesError_dedupData
    val pipelineError_dedupData

  output:
    tuple path ("${repRID}_sorted.deduped.bam"), path ("${repRID}_sorted.deduped.bam.bai") into dedupBam
    tuple path ("${repRID}_sorted.deduped.*.bam"), path ("${repRID}_sorted.deduped.*.bam.bai") into dedupChrBam
    path ("*.deduped.Metrics.txt") into dedupQC

  when:
    fastqCountError_dedupData == 'false'
    fastqReadError_dedupData == 'false'
    fastqFileError_dedupData == 'false'
    speciesError_dedupData == 'false'
    pipelineError_dedupData == 'false'

  script:
    """
    hostname > ${repRID}.dedup.log
    ulimit -a >> ${repRID}.dedup.log

    # remove duplicated reads using Picard's MarkDuplicates
    echo -e "LOG: deduplicating reads" >> ${repRID}.dedup.log
    java -jar /picard/build/libs/picard.jar MarkDuplicates I=${bam} O=${repRID}.deduped.bam M=${repRID}.deduped.Metrics.txt REMOVE_DUPLICATES=true
    echo -e "LOG: deduplicated" >> ${repRID}.dedup.log

    # sort the bam file using Samtools
    echo -e "LOG: sorting the bam file" >> ${repRID}.dedup.log
    samtools sort -@ `nproc` -O BAM -o ${repRID}_sorted.deduped.bam ${repRID}.deduped.bam

    # index the sorted bam using Samtools
    echo -e "LOG: indexing sorted bam file" >> ${repRID}.dedup.log
    samtools index -@ `nproc` -b ${repRID}_sorted.deduped.bam ${repRID}_sorted.deduped.bam.bai

    # split the deduped BAM file for multi-threaded tin calculation
    for i in `samtools view ${repRID}_sorted.deduped.bam | cut -f3 | grep -o chr.[0-9]* | sort | uniq`;
      do
      echo "echo \"LOG: splitting each chromosome into its own BAM and BAI files with Samtools\"; samtools view -b ${repRID}_sorted.deduped.bam \${i} 1>> ${repRID}_sorted.deduped.\${i}.bam; samtools index -@ `nproc` -b ${repRID}_sorted.deduped.\${i}.bam ${repRID}_sorted.deduped.\${i}.bam.bai"
    done | parallel -j `nproc` -k
    """
}

// Replicate dedup bam/bai for multiple process inputs
dedupBam.into {
  dedupBam_countData
  dedupBam_makeBigWig
  dedupBam_dataQC
  dedupBam_uploadProcessedFile
}

/*
 *makeBigWig: make BigWig files for output
*/
process makeBigWig {
  tag "${repRID}"
  publishDir "${outDir}/bigwig", mode: 'copy', pattern: "${repRID}.bw"

  input:
    tuple path (bam), path (bai) from dedupBam_makeBigWig
    val fastqCountError_makeBigWig
    val fastqReadError_makeBigWig
    val fastqFileError_makeBigWig
    val speciesError_makeBigWig
    val pipelineError_makeBigWig

  output:
    path ("${repRID}_sorted.deduped.bw") into bigwig

  when:
    fastqCountError_makeBigWig == 'false'
    fastqReadError_makeBigWig == 'false'
    fastqFileError_makeBigWig == 'false'
    speciesError_makeBigWig == 'false'
    pipelineError_makeBigWig == 'false'

  script:
    """
    hostname > ${repRID}.makeBigWig.log
    ulimit -a >> ${repRID}.makeBigWig.log

    # create bigwig
    echo -e "LOG: creating bibWig" >> ${repRID}.makeBigWig.log
    bamCoverage -p `nproc` -b ${bam} -o ${repRID}_sorted.deduped.bw
    echo -e "LOG: created" >> ${repRID}.makeBigWig.log
    """
}

/*
 *countData: count data and calculate tpm
*/
process countData {
  tag "${repRID}"
  publishDir "${outDir}/count", mode: 'copy', pattern: "${repRID}*_tpmTable.csv"

  input:
    path script_calculateTPM
    path script_convertGeneSymbols
    tuple path (bam), path (bai) from dedupBam_countData
    path ref from reference_countData
    val ends from endsInfer_countData
    val stranded from strandedInfer_countData
    val fastqCountError_countData
    val fastqReadError_countData
    val fastqFileError_countData
    val speciesError_countData
    val pipelineError_countData

  output:
    path ("*_tpmTable.csv") into counts
    path ("*_countData.summary") into countsQC
    path ("assignedReads.csv") into assignedReadsInfer_fl

  when:
    fastqCountError_countData == 'false'
    fastqReadError_countData == 'false'
    fastqFileError_countData == 'false'
    speciesError_countData == 'false'
    pipelineError_countData == 'false'

  script:
    """
    hostname > ${repRID}.countData.log
    ulimit -a >> ${repRID}.countData.log

    # determine strandedness and setup strandig for countData
    stranding=0;
    if [ "${stranded}" == "unstranded" ]
    then
      stranding=0
      echo -e "LOG: strandedness set to unstranded [0]" >> ${repRID}.countData.log
    elif [ "${stranded}" == "forward" ]
    then
      stranding=1
      echo -e "LOG: strandedness set to forward stranded [1]" >> ${repRID}.countData.log
    elif [ "${stranded}" == "reverse" ]
    then
      stranding=2
      echo -e "LOG: strandedness set to reverse stranded [2]" >> ${repRID}.countData.log
    fi

    # run featureCounts
    echo -e "LOG: counting ${ends} features" >> ${repRID}.countData.log
    if [ "${ends}" == "se" ]
    then
      featureCounts -T `nproc` -a ./genome.gtf -G ./genome.fna -g 'gene_name' --extraAttributes 'gene_id' -o ${repRID}_countData -s \${stranding} -R SAM --primary --ignoreDup ${repRID}_sorted.deduped.bam
    elif [ "${ends}" == "pe" ]
    then
      featureCounts -T `nproc` -a ./genome.gtf -G ./genome.fna -g 'gene_name' --extraAttributes 'gene_id' -o ${repRID}_countData -s \${stranding} -p -B -R SAM --primary --ignoreDup ${repRID}_sorted.deduped.bam
    fi
    echo -e "LOG: counted" >> ${repRID}.countData.log

    # extract assigned reads
    grep -m 1 'Assigned' *_countData.summary | grep -oe '\\([0-9.]*\\)' > assignedReads.csv

    # calculate TPM from the resulting countData table
    echo -e "LOG: calculating TPM with R" >> ${repRID}.countData.log
    Rscript ${script_calculateTPM} --count "${repRID}_countData"

    # convert gene symbols to Entrez id's
    echo -e "LOG: convert gene symbols to Entrez id's" >> ${repRID}.countData.log
    Rscript ${script_convertGeneSymbols} --repRID "${repRID}"
    """
}

// Extract number of assigned reads metadata into channel
assignedReadsInfer = Channel.create()
assignedReadsInfer_fl.splitCsv(sep: ",", header: false).separate(
  assignedReadsInfer
)

// Replicate inferred assigned reads for multiple process inputs
assignedReadsInfer.into {
  assignedReadsInfer_aggrQC
  assignedReadsInfer_uploadQC
}

/*
 *dataQC: calculate transcript integrity numbers (TIN) and bin as well as calculate innerdistance of PE replicates
*/
process dataQC {
  tag "${repRID}"

  input:
    path script_tinHist
    path ref from reference_dataQC
    tuple path (bam), path (bai) from dedupBam_dataQC
    tuple path (chrBam), path (chrBai) from dedupChrBam
    val ends from endsInfer_dataQC
    val fastqCountError_dataQC
    val fastqReadError_dataQC
    val fastqFileError_dataQC
    val speciesError_dataQC
    val pipelineError_dataQC

  output:
    path "${repRID}_tin.hist.tsv" into tinHist
    path "${repRID}_tin.med.csv" into  tinMedInfer_fl
    path "${repRID}_insertSize.inner_distance_freq.txt" into innerDistance

  when:
    fastqCountError_dataQC == 'false'
    fastqReadError_dataQC == 'false'
    fastqFileError_dataQC == 'false'
    speciesError_dataQC == 'false'
    pipelineError_dataQC == 'false'

  script:
    """
    hostname > ${repRID}.dataQC.log
    ulimit -a >> ${repRID}.dataQC.log

    # calcualte TIN values per feature on each chromosome
    echo -e  "geneID\tchrom\ttx_start\ttx_end\tTIN" > ${repRID}_sorted.deduped.tin.xls
    for i in `cat ./genome.bed | cut -f1 | grep -o chr.[0-9]* | sort | uniq`; do
      echo "echo \"LOG: running tin.py on \${i}\" >> ${repRID}.dataQC.log; tin.py -i ${repRID}_sorted.deduped.\${i}.bam  -r ./genome.bed; cat ${repRID}_sorted.deduped.\${i}.tin.xls | tr -s \"\\w\" \"\\t\" | grep -P \\\"\\\\t\${i}\\\\t\\\";";
    done | parallel -j `nproc` -k 1>> ${repRID}_sorted.deduped.tin.xls

    # bin TIN values
    echo -e "LOG: binning TINs" >> ${repRID}.dataQC.log
    python3 ${script_tinHist} -r ${repRID}
    echo -e "LOG: binned" >> ${repRID}.dataQC.log

    # calculate inner-distances for PE data
    if [ "${ends}" == "pe" ]
    then
      echo -e "LOG: calculating inner distances for ${ends}" >> ${repRID}.dataQC.log
      inner_distance.py -i "${bam}" -o ${repRID}_insertSize -r ./genome.bed
      echo -e "LOG: calculated" >> ${repRID}.dataQC.log
    elif [ "${ends}" == "se" ]
    then
      echo -e "LOG: creating dummy inner distance file for ${ends}" >> ${repRID}.dataQC.log
      touch ${repRID}_insertSize.inner_distance_freq.txt
    fi
    """
}

// Extract median TIN metadata into channel
tinMedInfer = Channel.create()
tinMedInfer_fl.splitCsv(sep: ",", header: false).separate(
  tinMedInfer
)

// Replicate inferred median TIN for multiple process inputs
tinMedInfer.into {
  tinMedInfer_aggrQC
  tinMedInfer_uploadQC
}

/*
 *aggrQC: aggregate QC from processes as well as metadata and run MultiQC
*/
process aggrQC {
  tag "${repRID}"
  publishDir "${outDir}/report", mode: 'copy', pattern: "${repRID}.multiqc.html"
  publishDir "${outDir}/qc", mode: 'copy', pattern: "${repRID}.multiqc_data.json"

  input:
    path multiqcConfig
    path bicfLogo
    path softwareReferences
    path softwareVersions
    path fastqc
    path trimQC
    path alignQC
    path dedupQC
    path countsQC
    path innerDistance
    path tinHist
    path alignSampleQCs from alignSampleQC_aggrQC.collect()
    path inferExperiment
    val endsManual from endsManual_aggrQC
    val endsM from endsMeta_aggrQC
    val strandedM from strandedMeta_aggrQC
    val spikeM from spikeMeta_aggrQC
    val speciesM from speciesMeta_aggrQC
    val endsI from endsInfer_aggrQC
    val strandedI from strandedInfer_aggrQC
    val spikeI from spikeInfer_aggrQC
    val speciesI from speciesInfer_aggrQC
    val readLengthM from readLengthMeta
    val readLengthI from readLengthInfer_aggrQC
    val rawReadsI from rawReadsInfer_aggrQC
    val assignedReadsI from assignedReadsInfer_aggrQC
    val tinMedI from tinMedInfer_aggrQC
    val studyRID from studyRID_aggrQC
    val expRID from expRID_aggrQC
    val fastqCountError_aggrQC
    val fastqReadError_aggrQC
    val fastqFileError_aggrQC
    val speciesError_aggrQC
    val pipelineError_aggrQC

  output:
    path "${repRID}.multiqc.html" into multiqc
    path "${repRID}.multiqc_data.json" into multiqcJSON

  when:
    fastqCountError_aggrQC == 'false'
    fastqReadError_aggrQC == 'false'
    fastqFileError_aggrQC == 'false'
    speciesError_aggrQC == 'false'
    pipelineError_aggrQC == 'false'

  script:
    """
    hostname > ${repRID}.aggrQC.log
    ulimit -a >> ${repRID}.aggrQC.log

    # make run table
    if [ "${params.inputBagForce}" == "" ] && [ "${params.fastqsForce}" == "" ] && [ "${params.speciesForce}" == "" ] && [ "${params.strandedForce}" == "" ] && [ "${params.spikeForce}" == "" ]
    then
      input="default"
    else
      input="override:"
      if [ "${params.inputBagForce}" != "" ]
      then
        input=\$(echo \${input} inputBag)
      fi
      if [ "${params.fastqsForce}" != "" ]
      then
        input=\$(echo \${input} fastq)
      fi
      if [ "${params.speciesForce}" != "" ]
      then
        input=\$(echo \${input} species)
      fi
      if [ "${params.strandedForce}" != "" ]
      then
        input=\$(echo \${input} stranded)
      fi
      if [ "${params.spikeForce}" != "" ]
      then
        input=\$(echo \${input} spike)
      fi
    fi
    echo -e "LOG: creating run table" >> ${repRID}.aggrQC.log
    echo -e "Session\tSession ID\tStart Time\tPipeline Version\tInput" > run.tsv
    echo -e "Session\t${workflow.sessionId}\t${workflow.start}\t${workflow.manifest.version}\t\${input}" >> run.tsv


    # make RID table
    echo -e "LOG: creating RID table" >> ${repRID}.aggrQC.log
    echo -e "Replicate\tReplicate RID\tExperiment RID\tStudy RID" > rid.tsv
    echo -e "Replicate\t${repRID}\t${expRID}\t${studyRID}" >> rid.tsv

    # make metadata table
    echo -e "LOG: creating metadata table" >> ${repRID}.aggrQC.log
    echo -e "Source\tSpecies\tEnds\tStranded\tSpike-in\tRaw Reads\tAssigned Reads\tMedian Read Length\tMedian TIN" > metadata.tsv
    echo -e "Submitter\t${speciesM}\t${endsM}\t${strandedM}\t${spikeM}\t-\t-\t'${readLengthM}'\t-" >> metadata.tsv
    if [ "${params.speciesForce}" == "" ]
    then
      input=\$(echo "Inferred\\t${speciesI}\\t")
    else
      input=\$(echo "Inferred\\t${speciesI} (FORCED)\\t")
    fi
    input=\$(echo \${input}"${endsI}\\t")
    if [ "${params.strandedForce}" == "" ]
    then
      input=\$(echo \${input}"${strandedI}\\t")
    else
      input=\$(echo \${input}"${strandedI} (FORCED)\\t")
    fi
    if [ "${params.spikeForce}" == "" ]
    then
      input=\$(echo \${input}"${spikeI}\\t-\\t-\\t-\\t-")
    else
      input=\$(echo \${input}"${spikeI} (FORCED)\\t-\\t-\\t-\\t-")
    fi
    echo -e \${input} >> metadata.tsv
    echo -e "Measured\t-\t${endsManual}\t-\t-\t'${rawReadsI}'\t'${assignedReadsI}'\t'${readLengthI}'\t'${tinMedI}'" >> metadata.tsv

    # make reference table
    echo -e "LOG: creating referencerun table" >> ${repRID}.aggrQC.log
    echo -e "Species\tGenome Reference Consortium Build\tGenome Reference Consortium Patch\tGENCODE Annotation Release" > reference.tsv
    echo -e "Human\tGRCh\$(echo `echo ${params.refHuVersion} | cut -d "." -f 1`)\t\$(echo `echo ${params.refHuVersion} | cut -d "." -f 2`)\t'\$(echo `echo ${params.refHuVersion} | cut -d "." -f 3 | sed "s/^v//"`)'" >> reference.tsv
    echo -e "Mouse\tGRCm\$(echo `echo ${params.refMoVersion} | cut -d "." -f 1`)\t\$(echo `echo ${params.refMoVersion} | cut -d "." -f 2`)\t'\$(echo `echo ${params.refMoVersion} | cut -d "." -f 3 | sed "s/^v//"`)'" >> reference.tsv

    # remove inner distance report if it is empty (SE repRID)
    echo -e "LOG: removing dummy inner distance file" >> ${repRID}.aggrQC.log
    if [ "${endsM}" == "se" ]
    then
      rm -f ${innerDistance}
    fi

    # run MultiQC
    echo -e "LOG: running multiqc" >> ${repRID}.aggrQC.log
    multiqc -c ${multiqcConfig} . -n ${repRID}.multiqc.html
    cp ${repRID}.multiqc_data/multiqc_data.json ${repRID}.multiqc_data.json

    if [ ${params.track} == true ]
    then
      curl -H 'Content-Type: application/json' -X PUT -d \
        @./${repRID}.multiqc_data.json \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/qc"
    fi
    """
}

/* 
 * uploadQC: uploads the mRNA QC
*/
process uploadQC {
  tag "${repRID}"

  input:
    path script_deleteEntry_uploadQC
    path script_uploadQC
    path credential, stageAs: "credential.json" from deriva_uploadQC
    val executionRunRID from executionRunRID_uploadQC
    val ends from endsInfer_uploadQC
    val stranded from strandedInfer_uploadQC
    val length from readLengthInfer_uploadQC
    val rawCount from rawReadsInfer_uploadQC
    val finalCount from assignedReadsInfer_uploadQC
    val tinMed from tinMedInfer_uploadQC
    val fastqCountError_uploadQC
    val fastqReadError_uploadQC
    val fastqFileError_uploadQC
    val speciesError_uploadQC
    val pipelineError_uploadQC

  output:
    path ("qcRID.csv") into qcRID_fl

  when:
    upload
    fastqCountError_uploadQC == 'false'
    fastqReadError_uploadQC == 'false'
    fastqFileError_uploadQC == 'false'
    speciesError_uploadQC == 'false'
    pipelineError_uploadQC == 'false'

  script:
    """
    hostname > ${repRID}.uploadQC.log
    ulimit -a >> ${repRID}.uploadQC.log

    if [ "${ends}" == "pe" ]
    then
      end="Paired End"
    elif [ "${ends}" == "se" ]
    then
      end="Single End"
    fi

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:mRNA_QC/Replicate=${repRID})
    if [ "\${exist}" != "[]" ]
    then
      rids=\$(echo \${exist} | grep -o '\\"RID\\":\\".\\{7\\}' | sed 's/^.\\{7\\}//')
      for rid in \${rids}
      do
        python3 ${script_deleteEntry_uploadQC} -r \${rid} -t mRNA_QC -o ${source} -c \${cookie}
        echo LOG: old mRNA QC RID deleted - \${rid} >> ${repRID}.uploadQC.log
      done
      echo LOG: all old mRNA QC RIDs deleted >> ${repRID}.uploadQC.log
    fi

    qc_rid=\$(python3 ${script_uploadQC} -r ${repRID} -e ${executionRunRID} -p "\${end}" -s ${stranded} -l ${length} -w ${rawCount} -f ${finalCount} -t ${tinMed} -o ${source} -c \${cookie} -u F)
    echo LOG: mRNA QC RID uploaded - \${qc_rid} >> ${repRID}.uploadQC.log

    echo "\${qc_rid}" > qcRID.csv
    """
}

/*
 *uploadProcessedFile: uploads the processed files
*/
process uploadProcessedFile {
  tag "${repRID}"
  publishDir "${outDir}/outputBag", mode: 'copy', pattern: "Replicate_${repRID}.outputBag.zip"

  input:
    path script_deleteEntry_uploadProcessedFile
    path credential, stageAs: "credential.json" from deriva_uploadProcessedFile
    path executionRunExportConfig
    path multiqc
    path multiqcJSON
    tuple path (bam),path (bai) from dedupBam_uploadProcessedFile
    path bigwig
    path counts
    val species from speciesInfer_uploadProcessedFile
    val studyRID from studyRID_uploadProcessedFile
    val expRID from expRID_uploadProcessedFile
    val executionRunRID from executionRunRID_uploadProcessedFile
    val fastqCountError_uploadProcessedFile
    val fastqReadError_uploadProcessedFile
    val fastqFileError_uploadProcessedFile
    val speciesError_uploadProcessedFile
    val pipelineError_uploadProcessedFile

  output:
    path ("${repRID}_Output_Bag.zip") into outputBag

  when:
    upload
    fastqCountError_uploadProcessedFile == 'false'
    fastqReadError_uploadProcessedFile == 'false'
    fastqFileError_uploadProcessedFile == 'false'
    speciesError_uploadProcessedFile == 'false'
    pipelineError_uploadProcessedFile == 'false'

  script:
    """
    hostname > ${repRID}.outputBag.log
    ulimit -a >> ${repRID}.outputBag.log

    mkdir -p ./deriva/Seq/pipeline/${studyRID}/${executionRunRID}/
    cp ${bam} ./deriva/Seq/pipeline/${studyRID}/${executionRunRID}/
    cp ${bai} ./deriva/Seq/pipeline/${studyRID}/${executionRunRID}/
    cp ${bigwig} ./deriva/Seq/pipeline/${studyRID}/${executionRunRID}/
    cp ${counts} ./deriva/Seq/pipeline/${studyRID}/${executionRunRID}/

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Processed_File/Replicate=${repRID})
    if [ "\${exist}" != "[]" ]
    then
      rids=\$(echo \${exist} | grep -o '\\"RID\\":\\".\\{7\\}' | sed 's/^.\\{7\\}//')
      for rid in \${rids}
      do
        python3 ${script_deleteEntry_uploadProcessedFile} -r \${rid} -t Processed_File -o ${source} -c \${cookie}
      done
      echo LOG: all old processed file RIDs deleted >> ${repRID}.uploadQC.log
    fi

    deriva-upload-cli --catalog 2 --token \${cookie:9} ${source} ./deriva
    echo LOG: processed files uploaded >> ${repRID}.outputBag.log

    deriva-download-cli --catalog 2 --token \${cookie:9} ${source} ${executionRunExportConfig} . rid=${executionRunRID}
    echo LOG: execution run bag downloaded >> ${repRID}.outputBag.log

    echo -e "### Run Details" >> runDetails.md
    echo -e "**Workflow URL:** https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq" >> runDetails.md
    echo -e "**Workflow Version:** ${workflow.manifest.version}" >> runDetails.md
    echo -e "**Description:** ${workflow.manifest.description}" >> runDetails.md
    if [ "${species}" == "Mus musculus" ]; then
      genome=\$(echo GRCm${refMoVersion} | cut -d '.' -f1)
      patch=\$(echo ${refMoVersion} | cut -d '.' -f2)
      annotation=\$(echo ${refMoVersion} | cut -d '.' -f3 | tr -d 'v')
    elif [ "${species}" == "Homo sapiens" ]; then
      genome=\$(echo GRCh${refHuVersion} | cut -d '.' -f1)
      patch=\$(echo ${refHuVersion} | cut -d '.' -f2)
      annotation=\$(echo ${refHuVersion} | cut -d '.' -f3 | tr -d 'v')
    fi
    echo -e "**Genome Assembly Version:** \${genome} patch \${patch}" >> runDetails.md
    echo -e "**Annotation Version:** GENCODE release \${annotation}" >> runDetails.md
    echo -e "**Run ID:** ${repRID}" >> runDetails.md
    echo LOG: runDetails.md created >> ${repRID}.outputBag.log

    unzip Execution_Run_${executionRunRID}.zip
    yr=\$(date +'%Y')
    mn=\$(date +'%m')
    dy=\$(date +'%d')
    mv Execution_Run_${executionRunRID} ${repRID}_Output_Bag_\${yr}\${mn}\${dy}
    loc=./${repRID}_Output_Bag/data/assets/Study/${studyRID}/Experiment/${expRID}/Replicate/${repRID}/Execution_Run/${executionRunRID}/Output_Files/
    mkdir -p \${loc}
    cp runDetails.md \${loc}
    cp ${multiqc} \${loc}
    cp ${multiqcJSON} \${loc}

    bdbag ./${repRID}_Output_Bag/ --update --archiver zip --debug
    echo LOG: output bag created >> ${repRID}.outputBag.log
    """
}

/* 
 * uploadOutputBag: uploads the output bag
*/
process uploadOutputBag {
  tag "${repRID}"

  input:
    path script_uploadOutputBag
    path credential, stageAs: "credential.json" from deriva_uploadOutputBag
    path outputBag
    val studyRID from studyRID_uploadOutputBag
    val executionRunRID from executionRunRID_uploadOutputBag
    val fastqCountError_uploadOutputBag
    val fastqReadError_uploadOutputBag
    val fastqFileError_uploadOutputBag
    val speciesError_uploadOutputBag
    val pipelineError_uploadOutputBag

  output:
    path ("outputBagRID.csv") into outputBagRID_fl

  when:
    upload
    fastqCountError_uploadOutputBag == 'false'
    fastqReadError_uploadOutputBag == 'false'
    fastqFileError_uploadOutputBag == 'false'
    speciesError_uploadOutputBag == 'false'
    pipelineError_uploadOutputBag == 'false'

  script:
    """
    hostname > ${repRID}.uploadOutputBag.log
    ulimit -a >> ${repRID}.uploadOutputBag.log

    yr=\$(date +'%Y')
    mn=\$(date +'%m')
    dy=\$(date +'%d')

    file=\$(basename -a ${outputBag})
    md5=\$(md5sum ./\${file} | awk '{ print \$1 }')
    echo LOG: ${repRID} output bag md5 sum - \${md5} >> ${repRID}.uploadOutputBag.log
    size=\$(wc -c < ./\${file})
    echo LOG: ${repRID} output bag size - \${size} bytes >> ${repRID}.uploadOutputBag.log
    
    loc=\$(deriva-hatrac-cli --host ${source} put ./\${file} /hatrac/resources/rnaseq/pipeline/output_bag/study/${studyRID}/replicate/${repRID}/\${file} --parents)
    echo LOG: output bag uploaded - \${loc} >> ${repRID}.uploadOutputBag.log
    # url-ify the location
    loc=\${loc//\\//%2F}
    loc=\${loc//:/%3A}
    loc=\${loc// /@20}

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Output_Bag/File_URL=\${loc})
    if [ "\${exist}" == "[]" ]
    then
      outputBag_rid=\$(python3 ${script_uploadOutputBag} -e ${executionRunRID} -f \${file} -l \${loc} -s \${md5} -b \${size} -o ${source} -c \${cookie} -u F)
      echo LOG: output bag RID uploaded - \${outputBag_rid} >> ${repRID}.uploadOutputBag.log
      rid=\${outputBag_rid}
    else
      exist=\$(echo \${exist} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
      exist=\${exist:8:-6}
      outputBag_rid=\$(python3 ${script_uploadOutputBag} -e ${executionRunRID} -o ${source} -c \${cookie} -u \${exist})
      echo LOG: output bag RID already exists - \${exist} >> ${repRID}.uploadOutputBag.log
      rid=\${exist}
    fi

    echo "\${rid}" > outputBagRID.csv
    """
}

// Extract output bag RID into channel
outputBagRID = Channel.create()
outputBagRID_fl.splitCsv(sep: ",", header: false).separate(
  outputBagRID
)

/* 
 * finalizeExecutionRun: finalizes the execution run
*/
process finalizeExecutionRun {
  tag "${repRID}"

  input:
    path script_uploadExecutionRun_finalizeExecutionRun
    path credential, stageAs: "credential.json" from deriva_finalizeExecutionRun
    val executionRunRID from executionRunRID_finalizeExecutionRun
    val inputBagRID from inputBagRID_finalizeExecutionRun
    val outputBagRID

  when:
    upload

  script:
    """
    hostname > ${repRID}.finalizeExecutionRun.log
    ulimit -a >> ${repRID}.finalizeExecutionRun.log

    executionRun=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Execution_Run/RID=${executionRunRID})
    workflow=\$(echo \${executionRun} | grep -o '\\"Workflow\\":.*\\"Reference' | grep -oP '(?<=\\"Workflow\\":\\").*(?=\\",\\"Reference)')
    genome=\$(echo \${executionRun} | grep -o '\\"Reference_Genome\\":.*\\"Input_Bag' | grep -oP '(?<=\\"Reference_Genome\\":\\").*(?=\\",\\"Input_Bag)')

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    rid=\$(python3 ${script_uploadExecutionRun_finalizeExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s Success -d 'Run Successful' -o ${source} -c \${cookie} -u ${executionRunRID})
    echo LOG: execution run RID marked as successful - \${rid} >> ${repRID}.finalizeExecutionRun.log

    if [ ${params.track} == true ]
    then
    dt=`date +%FT%T.%3N%:z`
      curl -H 'Content-Type: application/json' -X PUT -d \
        '{ \
          "ID": "${workflow.sessionId}", \
          "Complete": "'\${dt}'" \
        }' \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/track"
    fi
    """
}

// Combine errors
error_meta = fastqCountError_uploadQC_fail.ifEmpty(false).combine(fastqReadError_uploadQC_fail.ifEmpty(false).combine(fastqFileError_uploadQC_fail.ifEmpty(false).combine(speciesError_uploadQC_fail.ifEmpty(false).combine(pipelineError_uploadQC_fail.ifEmpty(false)))))
error_meta. into{
  error_failPreExecutionRun
  error_uploadQC_fail
}
errorDetails = fastqCountError_details.ifEmpty("").combine(fastqReadError_details.ifEmpty("").combine(fastqFileError_details.ifEmpty("").combine(speciesError_details.ifEmpty(""))))

/* 
 * failPreExecutionRun_fastq: fail the execution run prematurely for fastq errors
*/
process failPreExecutionRun {
  tag "${repRID}"

  input:
    path script_uploadExecutionRun from script_uploadExecutionRun_failPreExecutionRun
    path credential, stageAs: "credential.json" from deriva_failPreExecutionRun
    val spike from spikeMeta_failPreExecutionRun
    val species from speciesMeta_failPreExecutionRun
    val inputBagRID from inputBagRID_failPreExecutionRun
    tuple val (fastqCountError), val (fastqReadError), val (fastqFileError), val (speciesError), val (pipelineError) from error_failPreExecutionRun
    tuple val (fastqCountError_details), val (fastqReadError_details), val (fastqFileError_details), val (speciesError_details) from errorDetails

  output:
    path ("executionRunRID.csv") into executionRunRID_preFail_fl

  when:
    upload
    fastqCountError == 'true' || fastqReadError == 'true' || fastqFileError == 'true' || speciesError == 'true'

  script:
    """
    hostname > ${repRID}.failPreExecutionRun.log
    ulimit -a >> ${repRID}.failPreExecutionRun.log

    errorDetails=""
    if [ ${fastqCountError} == true ]
    then
      errorDetails=\$(echo ${fastqCountError_details}"\\n")
    elif [ ${fastqReadError} == true ]
    then
      errorDetails=\$(echo \$(errorDetails)${fastqReadError_details}"\\n")
    elif [ ${fastqFileError} == true ]
    then
      errorDetails=\$(echo \$(errorDetails)${fastqReadError_details}"\\n")
    elif [ ${speciesError} == true ]
    then
      errorDetails=\$(echo \$(errorDetails)${fastqReadError_details}"\\n")
    fi

    echo LOG: searching for workflow RID - BICF mRNA ${workflow.manifest.version} >> ${repRID}.failPreExecutionRun.log
    workflow=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Workflow/Name=BICF%20mRNA%20Replicate/Version=${workflow.manifest.version})
    workflow=\$(echo \${workflow} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
    workflow=\${workflow:7:-6}
    echo LOG: workflow RID extracted - \${workflow} >> ${repRID}.failPreExecutionRun.log

    if [ "${species}" == "Homo sapiens" ]
    then
      genomeName=\$(echo GRCh${refHuVersion})
    elif [ "${species}" == "Mus musculus" ]
    then
      genomeName=\$(echo GRCm${refMoVersion})
    fi
    if [ "${spike}" == "true" ]
    then
      genomeName=\$(echo \${genomeName}-S)
    fi
    echo LOG: searching for genome name - \${genomeName} >> ${repRID}.failPreExecutionRun.log
    genome=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Name=\${genomeName})
    genome=\$(echo \${genome} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
    genome=\${genome:7:-6}
    echo LOG: genome RID extracted - \${genome} >> ${repRID}.failPreExecutionRun.log

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Execution_Run/Workflow=\${workflow}/Replicate=${repRID}/Input_Bag=${inputBagRID})
    echo \${exist} >> ${repRID}.failPreExecutionRun.log
    if [ "\${exist}" == "[]" ]
    then
      rid=\$(python3 ${script_uploadExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s Error -d "\${errorDetails}" -o ${source} -c \${cookie} -u F)
      echo LOG: execution run RID uploaded - \${rid} >> ${repRID}.failPreExecutionRun.log
    else
      rid=\$(echo \${exist} | grep -o '\\"RID\\":\\".*\\",\\"RCT')
      rid=\${rid:7:-6}
      echo \${rid} >> ${repRID}.failPreExecutionRun.log
      executionRun_rid=\$(python3 ${script_uploadExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s Error -d "\${errorDetails}" -o ${source} -c \${cookie} -u \${rid})
      echo LOG: execution run RID updated - \${executionRun_rid} >> ${repRID}.failPreExecutionRun.log
    fi

    echo "\${rid}" > executionRunRID.csv

    if [ ${params.track} == true ]
    then
    dt=`date +%FT%T.%3N%:z`
      curl -H 'Content-Type: application/json' -X PUT -d \
        '{ \
          "ID": "${workflow.sessionId}", \
          "ExecutionRunRID": "'\${rid}'", \
          "Failure": "'\${dt}'" \
        }' \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/track"
    fi
  """
}
// Extract execution run RID into channel
executionRunRID_preFail = Channel.create()
executionRunRID_preFail_fl.splitCsv(sep: ",", header: false).separate(
  executionRunRID_preFail
)

failExecutionRunRID = executionRunRID_fail.ifEmpty('').mix(executionRunRID_preFail.ifEmpty('')).filter { it != "" }

/* 
 * failExecutionRun: fail the execution run
*/
process failExecutionRun {
  tag "${repRID}"

  input:
    path script_uploadExecutionRun_failExecutionRun
    path credential, stageAs: "credential.json" from deriva_failExecutionRun
    val executionRunRID from executionRunRID_failExecutionRun
    val inputBagRID from inputBagRID_failExecutionRun
    val endsMeta from endsMeta_failExecutionRun
    val endsRaw
    val strandedMeta from strandedMeta_failExecutionRun
    val spikeMeta from spikeMeta_failExecutionRun
    val speciesMeta from speciesMeta_failExecutionRun
    val endsInfer from endsInfer_failExecutionRun
    val strandedInfer from strandedInfer_failExecutionRun
    val spikeInfer from spikeInfer_failExecutionRun
    val speciesInfer from speciesInfer_failExecutionRun
    val pipelineError from pipelineError_failExecutionRun
    val pipelineError_ends
    val pipelineError_stranded
    val pipelineError_spike
    val pipelineError_species

  when:
    upload
    pipelineError == 'true'

  script:
    """
    hostname > ${repRID}.failExecutionRun.log
    ulimit -a >> ${repRID}.failExecutionRun.log

    executionRun=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:Execution_Run/RID=${executionRunRID})
    workflow=\$(echo \${executionRun} | grep -o '\\"Workflow\\":.*\\"Reference' | grep -oP '(?<=\\"Workflow\\":\\").*(?=\\",\\"Reference)')
    genome=\$(echo \${executionRun} | grep -o '\\"Reference_Genome\\":.*\\"Input_Bag' | grep -oP '(?<=\\"Reference_Genome\\":\\").*(?=\\",\\"Input_Bag)')

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    errorDetails=""
    if [ ${pipelineError} == false ]
    then
      rid=\$(python3 ${script_uploadExecutionRun_failExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s Success -d 'Run Successful' -o ${source} -c \${cookie} -u ${executionRunRID})
      echo LOG: execution run RID marked as successful - \${rid} >> ${repRID}.failExecutionRun.log
    else
      pipelineError_details=\$(echo "**Submitted metadata does not match inferred:**\\n")
      pipelineError_details=\$(echo \${pipelineError_details}"|Metadata|Submitted value|Inferred value|\\n")
      pipelineError_details=\$(echo \${pipelineError_details}"|:-:|-:|-:|\\n")
      if ${pipelineError_ends}
      then
        if [ "${endsInfer}" == "se" ]
        then
          endInfer="Single End"
        elif [ "${endsInfer}" == "pe" ]
        then
          endInfer="Paired End"
        else
          endInfer="unknown"
        fi
        pipelineError_details=\$(echo \${pipelineError_details}"|Paired End|${endsRaw}|"\${endInfer}"|\\n")
      fi
      if ${pipelineError_stranded}
      then
        pipelineError_details=\$(echo \${pipelineError_details}"|Strandedness|${strandedMeta}|${strandedInfer}|\\n")
      fi
      if ${pipelineError_spike}
      then
        pipelineError_details=\$(echo \${pipelineError_details}"|Used Spike Ins|${spikeMeta}|${spikeInfer}|\\n")
      fi
      if ${pipelineError_species}
      then
        pipelineError_details=\$(echo \${pipelineError_details}"|Species|${speciesMeta}|${speciesInfer}|\\n")
      fi
      pipelineError_details=\${pipelineError_details::-2}
      rid=\$(python3 ${script_uploadExecutionRun_failExecutionRun} -r ${repRID} -w \${workflow} -g \${genome} -i ${inputBagRID} -s Error -d "\${pipelineError_details}" -o ${source} -c \${cookie} -u ${executionRunRID})
      echo LOG: execution run RID marked as error - \${rid} >> ${repRID}.failExecutionRun.log
    fi
    
    if [ ${params.track} == true ]
    then
      dt=`date +%FT%T.%3N%:z`
      curl -H 'Content-Type: application/json' -X PUT -d \
        '{ \
          "ID": "${workflow.sessionId}", \
          "ExecutionRunRID": "'\${rid}'", \
          "Failure": "'\${dt}'" \
        }' \
        "https://9ouc12dkwb.execute-api.us-east-2.amazonaws.com/prod/db/track"
    fi
  """
}

/* 
 * uploadQC_fail: uploads the mRNA QC on failed execution run
*/
process uploadQC_fail {
  tag "${repRID}"

  input:
    path script_deleteEntry_uploadQC_fail
    path script_uploadQC_fail
    path credential, stageAs: "credential.json" from deriva_uploadQC_fail
    val executionRunRID from failExecutionRunRID
    tuple val (fastqCountError), val (fastqReadError), val (fastqFileError), val (speciesError), val (pipelineError) from error_uploadQC_fail

  when:
    upload
    fastqCountError == 'true' || fastqReadError == 'true' || fastqFileError == 'true' || speciesError == 'true' || pipelineError == 'true'

  script:
    """
    hostname > ${repRID}.uploadQC.log
    ulimit -a >> ${repRID}.uploadQC.log

    cookie=\$(cat credential.json | grep -A 1 '\\"${source}\\": {' | grep -o '\\"cookie\\": \\".*\\"')
    cookie=\${cookie:11:-1}

    exist=\$(curl -s https://${source}/ermrest/catalog/2/entity/RNASeq:mRNA_QC/Replicate=${repRID})
    if [ "\${exist}" != "[]" ]
    then
      rids=\$(echo \${exist} | grep -o '\\"RID\\":\\".\\{7\\}' | sed 's/^.\\{7\\}//')
      for rid in \${rids}
      do
        python3 ${script_deleteEntry_uploadQC_fail} -r \${rid} -t mRNA_QC -o ${source} -c \${cookie}
        echo LOG: old mRNA QC RID deleted - \${rid} >> ${repRID}.uploadQC.log
      done
      echo LOG: all old mRNA QC RIDs deleted >> ${repRID}.uploadQC.log
    fi

    qc_rid=\$(python3 ${script_uploadQC_fail} -r ${repRID} -e ${executionRunRID} -o ${source} -c \${cookie} -u E)
    echo LOG: mRNA QC RID uploaded - \${qc_rid} >> ${repRID}.uploadQC.log

    echo "\${qc_rid}" > qcRID.csv
    """
}


workflow.onError = {
  subject = "$workflow.manifest.name FAILED: $params.repRID"

  def msg = """\

      Pipeline error summary
      ---------------------------
      RID         : ${params.repRID}
      Version     : ${workflow.manifest.version}
      Duration    : ${workflow.duration}
      Nf Version  : ${workflow.nextflow.version}
      Message     : ${workflow.errorMessage}
      exit status : ${workflow.exitStatus}
      """
      .stripIndent()
  if (email != '') {
    sendMail(to: email, subject: subject , body: msg)
  }
}
