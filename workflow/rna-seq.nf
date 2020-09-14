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
params.repRID = "Q-Y5JA"
params.source = "dev"
params.refMoVersion = "38.p6.vM22"
params.refHuVersion = "38.p12.v31"
params.refERCCVersion = "92"
params.outDir = "${baseDir}/../output"

// Define override input variable
params.refSource = "biohpc"
params.inputBagForce = ""
params.fastqsForce = ""
params.speciesForce = ""

// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
deriva.into {
  deriva_getBag
  deriva_getRefInfer
  deriva_getRef
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
inputBagForce = params.inputBagForce
fastqsForce = params.fastqsForce
speciesForce = params.speciesForce

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")
if (params.source == "dev") {
  source = "dev.gudmap.org"
} else if (params.source == "staging") {
  source = "staging.gudmap.org"
} else if (params.source == "production") {
  source = "www.gudmap.org"
}
if (params.refSource == "biohpc") {
  referenceBase = "/project/BICF/BICF_Core/shared/gudmap/references"
//} else if (params.refSource == "aws") {
//  referenceBase = "s3://bicf-references"
} else if (params.refSource == "datahub") {
  referenceBase = "dev.gudmap.org"
}
referenceInfer = Channel.fromList(["ERCC","GRCh","GRCm"])
multiqcConfig = Channel.fromPath("${baseDir}/conf/multiqc_config.yaml")
bicfLogo = Channel.fromPath("${baseDir}/../docs/bicf_logo.png")

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbagFetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parseMeta.py")
script_inferMeta = Channel.fromPath("${baseDir}/scripts/inferMeta.sh")
script_refDataInfer = Channel.fromPath("${baseDir}/scripts/extractRefData.py")
script_refData = Channel.fromPath("${baseDir}/scripts/extractRefData.py")
script_calculateTPM = Channel.fromPath("${baseDir}/scripts/calculateTPM.R")
script_convertGeneSymbols = Channel.fromPath("${baseDir}/scripts/convertGeneSymbols.R")
script_tinHist = Channel.fromPath("${baseDir}/scripts/tinHist.py")


/*
 * trackStart: track start of pipeline
 */
params.ci = false
params.dev = false
process trackStart {
  container 'docker://bicf/bicfbase:2.1.0'
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
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
process getBag {
  tag "${repRID}"
  publishDir "${outDir}/inputBag", mode: 'copy', pattern: "Replicate_*.zip"

  input:
    path credential, stageAs: "credential.json" from deriva_getBag
    path derivaConfig

  output:
    path ("Replicate_*.zip") into bag

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
    deriva-download-cli ${source} --catalog 2 ${derivaConfig} . rid=${repRID}
    echo -e "LOG: fetched" >> ${repRID}.getBag.log
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

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID}"

  input:
    path script_bdbagFetch
    path cookies, stageAs: "deriva-cookies.txt" from bdbag
    path inputBag

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    path ("**/File.csv") into fileMeta
    path ("**/Experiment Settings.csv") into experimentSettingsMeta
    path ("**/Experiment.csv") into experimentMeta

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
    replicate=\$(basename "${inputBag}" | cut -d "." -f1)
    echo -e "LOG: bag replicate name \${replicate}" >> ${repRID}.getData.log
    
    # unzip bag
    echo -e "LOG: unzipping replicate bag" >> ${repRID}.getData.log
    unzip ${inputBag}
    echo -e "LOG: unzipped" >> ${repRID}.getData.log
    
    # bag fetch fastq's only and rename by repRID
    echo -e "LOG: fetching replicate bdbag" >> ${repRID}.getData.log
    sh ${script_bdbagFetch} \${replicate} ${repRID}
    echo -e "LOG: fetched" >> ${repRID}.getData.log
    """
}

// Set raw fastq to downloaded or forced input and replicate them for multiple process inputs
if (fastqsForce != "") {
  Channel
    .fromPath(fastqsForce)
    .ifEmpty { exit 1, "override inputBag file not found: ${fastqsForce}" }
    .collect().into {
      fastqs_trimData
      fastqs_fastqc
    }
} else {
  fastqs.into {
    fastqs_trimData
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

  output:
    path "design.csv" into metadata

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
    endsMeta=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p endsMeta)
    echo -e "LOG: endedness metadata parsed: \${endsMeta}" >> ${repRID}.parseMetadata.log
    
    # ganually get endness
    endsManual=\$(python3 ${script_parseMeta} -r ${repRID} -m "${file}" -p endsManual)
    echo -e "LOG: endedness manually detected: \${endsManual}" >> ${repRID}.parseMetadata.log

    # get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p stranded)
    echo -e "LOG: strandedness metadata parsed: \${stranded}" >> ${repRID}.parseMetadata.log
    
    # get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p spike)
    echo -e "LOG: spike-in metadata parsed: \${spike}" >> ${repRID}.parseMetadata.log
    
    # get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experiment}" -p species)
    echo -e "LOG: species metadata parsed: \${species}" >> ${repRID}.parseMetadata.log

    # get read length metadata
    readLength=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettings}" -p readLength)
    if [ "\${readLength}" = "nan"]
    then
      readLength="NA"
    fi
    echo -e "LOG: read length metadata parsed: \${readLength}" >> ${repRID}.parseMetadata.log

    # save design file
    echo -e "\${endsMeta},\${endsManual},\${stranded},\${spike},\${species},\${readLength},\${exp},\${study}" > design.csv
    """
}

// Split metadata into separate channels
endsMeta = Channel.create()
endsManual = Channel.create()
strandedMeta = Channel.create()
spikeMeta = Channel.create()
speciesMeta = Channel.create()
readLengthMeta = Channel.create()
expRID = Channel.create()
studyRID = Channel.create()
metadata.splitCsv(sep: ",", header: false).separate(
  endsMeta,
  endsManual,
  strandedMeta,
  spikeMeta,
  speciesMeta,
  readLengthMeta,
  expRID,
  studyRID
)
// Replicate metadata for multiple process inputs
endsManual.into {
  endsManual_trimData
  endsManual_downsampleData
  endsManual_alignSampleData
  endsManual_aggrQC
}


/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"

  input:
    val ends from endsManual_trimData
    path (fastq) from fastqs_trimData

  output:
    path ("*.fq.gz") into fastqsTrim
    path ("*_trimming_report.txt") into trimQC
    path ("readLength.csv") into inferMetadata_readLength

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
    echo -e "\${readLength}" > readLength.csv
    """
}

// Extract calculated read length metadata into channel
readLengthInfer = Channel.create()
inferMetadata_readLength.splitCsv(sep: ",", header: false).separate(
  readLengthInfer
)

// Replicate trimmed fastq's
fastqsTrim.into {
  fastqsTrim_alignData
  fastqsTrim_downsampleData
}

/*
  * getRefInfer: dowloads appropriate reference for metadata inference
*/  
process getRefInfer {
  tag "${refName}"

  input:
    path credential, stageAs: "credential.json" from deriva_getRefInfer
    path script_refDataInfer
    val refName from referenceInfer

  output:
    tuple val (refName), path ("hisat2", type: 'dir'), path ("*.fna"), path ("*.gtf")  into refInfer
    path ("${refName}", type: 'dir') into bedInfer
 
  script:
    """
    hostname > ${repRID}.${refName}.getRefInfer.log
    ulimit -a >> ${repRID}.${refName}.getRefInfer.log

    # link credential file for authentication
    echo -e "LOG: linking deriva credentials" >> ${repRID}.getRefInfer.log
    mkdir -p ~/.deriva
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json
    echo -e "LOG: linked" >> ${repRID}.getRefInfer.log

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
    mkdir ${refName}

    # retreive appropriate reference appropriate location
    echo -e "LOG: fetching ${refName} reference files from ${referenceBase}" >> ${repRID}.${refName}.getRefInfer.log
    if [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references" ]
    then
      ln -s "\${references}"/hisat2
      ln -s "\${references}"/bed ${refName}/bed
      ln -s "\${references}"/genome.fna
      ln -s "\${references}"/genome.gtf
    #elif [ ${referenceBase} == "s3://bicf-references" ]
    #then
    #  aws s3 cp "\${references}"/hisat2 ./hisat2 --recursive
    #  aws s3 cp "\${references}"/bed ./${refName}/bed --recursive
    #  aws s3 cp "\${references}"/genome.fna ./
    #  aws s3 cp "\${references}"/genome.gtf ./
    elif [ ${referenceBase} == "dev.gudmap.org" ]
    then
      GRCv=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f1)
      GRCp=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f2)
      GENCODE=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f3)
      query=\$(echo 'https://${referenceBase}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Reference_Version='\${GRCv}'.'\${GRCp}'/Annotation_Version=GENCODE%20'\${GENCODE})
      curl --request GET \${query} > refQuery.json
      refURL=\$(python extractRefData.py --returnParam URL)
      loc=\$(dirname \${refURL})
      fName=\$(python extractRefData.py --returnParam fName)
      fName=\${fName%.*}
      if [ "\${loc}" = "/hatrac/*" ]; then echo "LOG: Reference not present in hatrac"; exit 1; fi
      filename=\$(echo \$(basename \${refURL}) | grep -oP '.*(?=:)')
      deriva-hatrac-cli --host ${referenceBase} get \${refURL}
      unzip \$(basename \${refURL})
      mv \${fName}/* .
      mv bed ${refName}/
    fi
    echo -e "LOG: fetched" >> ${repRID}.${refName}.getRefInfer.log

    # make blank bed folder for ERCC
    echo -e "LOG: making dummy bed folder for ERCC" >> ${repRID}.${refName}.getRefInfer.log
    if [ "${refName}" == "ERCC" ]
    then
      rm -rf ${refName}/bed
      mkdir ${refName}/bed
      touch ${refName}/bed/temp
    fi
    """
}

/*
 * downsampleData: downsample fastq's for metadata inference
 */
process downsampleData {
  tag "${repRID}"

  input:
    val ends from endsManual_downsampleData
    path fastq from fastqsTrim_downsampleData

  output:
    path ("sampled.1.fq") into fastqs1Sample
    path ("sampled.2.fq") into fastqs2Sample

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
inferInput = endsManual_alignSampleData.combine(refInfer.combine(fastqs1Sample.collect().combine(fastqs2Sample.collect())))

/*
 * alignSampleData: aligns the downsampled reads to a reference database
*/
process alignSampleData {
  tag "${ref}"

  input:
    tuple val (ends), val (ref), path (hisat2), path (fna), path (gtf), path (fastq1), path (fastq2) from inferInput

  output:
    path ("${ref}.sampled.sorted.bam") into sampleBam
    path ("${ref}.sampled.sorted.bam.bai") into sampleBai
    path ("${ref}.alignSampleSummary.txt") into alignSampleQC

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
    samtools sort -@ `nproc` -O BAM -o ${ref}.sampled.sorted.bam ${ref}.sampled.bam

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

  output:
    path "infer.csv" into inferMetadata
    path "${repRID}.infer_experiment.txt" into inferExperiment

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
      spike="yes"
    else
      spike="no"
    fi
    echo -e "LOG: inference of strandedness results is: \${spike}" >> ${repRID}.inferMetadata.log

    # determine species
    if [ 1 -eq \$(echo \$(expr \${align_hu} ">=" 25)) ] && [ 1 -eq \$(echo \$(expr \${align_mo} "<" 25)) ]
    then
      species="Homo sapiens"
      bam="GRCh.sampled.sorted.bam"
      bed="./GRCh/bed/genome.bed"
    elif [ 1 -eq \$(echo \$(expr \${align_mo} ">=" 25)) ] && [ 1 -eq \$(echo \$(expr \${align_hu} "<" 25)) ]
    then
      species="Mus musculus"
      bam="GRCm.sampled.sorted.bam"
      bed="./GRCm/bed/genome.bed"
    else
      echo -e "LOG: ERROR - inference of species returns an ambiguous result: hu=\${align_hu} mo=\${align_mo}" >> ${repRID}.inferMetadata.log
      if [ "${speciesForce}" == "" ]
      then
        exit 1
      fi
    fi
    if [ "${speciesForce}" != "" ]
    then
      echo -e "LOG: species overridden to: ${speciesForce}"
      species="${speciesForce}"
      if [ "${speciesForce}" == "Homo sapiens" ]
      then
        bam="GRCh.sampled.sorted.bam"
        bed="./GRCh/bed/genome.bed"
      elif [ "${speciesForce}" == "Mus musculus" ]
      then
        bam="GRCm.sampled.sorted.bam"
        bed="./GRCm/bed/genome.bed"
      fi
    fi
    echo -e "LOG: inference of species results in: \${species}" >> ${repRID}.inferMetadata.log

    # infer experimental setting from dedup bam
    echo -e "LOG: infer experimental setting from dedup bam" >> ${repRID}.inferMetadata.log
    infer_experiment.py -r "\${bed}" -i "\${bam}" 1>> ${repRID}.infer_experiment.txt
    echo -e "LOG: infered" >> ${repRID}.inferMetadata.log

    ended=`bash inferMeta.sh endness ${repRID}.infer_experiment.txt`
    fail=`bash inferMeta.sh fail ${repRID}.infer_experiment.txt`
    if [ \${ended} == "PairEnd" ] 
    then
      ends="pe"
      percentF=`bash inferMeta.sh pef ${repRID}.infer_experiment.txt`
      percentR=`bash inferMeta.sh per ${repRID}.infer_experiment.txt`
    elif [ \${ended} == "SingleEnd" ]
    then
      ends="se"
      percentF=`bash inferMeta.sh sef ${repRID}.infer_experiment.txt`
      percentR=`bash inferMeta.sh ser ${repRID}.infer_experiment.txt`
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

    # write infered metadata to file
    echo "\${ends},\${stranded},\${spike},\${species},\${align_ercc},\${align_hu},\${align_mo},\${percentF},\${percentR},\${fail}" 1>> infer.csv
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
inferMetadata.splitCsv(sep: ",", header: false).separate(
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
  endsInfer_alignData
  endsInfer_countData
  endsInfer_dataQC
  endsInfer_aggrQC
}
strandedInfer.into {
  strandedInfer_alignData
  strandedInfer_countData
  strandedInfer_aggrQC
}
spikeInfer.into{
  spikeInfer_getRef
  spikeInfer_aggrQC
}
speciesInfer.into {
  speciesInfer_getRef
  speciesInfer_aggrQC
}


/*
  * getRef: downloads appropriate reference
*/
process getRef {
  tag "${species}"

  input:
    path credential, stageAs: "credential.json" from deriva_getRef
    path script_refData
    val spike from spikeInfer_getRef
    val species from speciesInfer_getRef

  output:
    tuple path ("hisat2", type: 'dir'), path ("bed", type: 'dir'), path ("*.fna"), path ("*.gtf"), path ("geneID.tsv"), path ("Entrez.tsv")  into reference
 
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
    elif [ '${species}' == "Homo sapiens" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
    else
      echo -e "LOG: ERROR - References could not be set!\nSpecies reference found: ${species}" >> ${repRID}.getRef.log
      exit 1
    fi
    if [ "${spike}" == "yes" ]
    then
      references=\$(echo \${reference}-S/)
    elif [ "${spike}" == "no" ]
    then
      reference=\$(echo \${references}/)
    fi
    echo -e "LOG: species set to \${references}" >> ${repRID}.getRef.log

    # retreive appropriate reference appropriate location
    echo -e "LOG: fetching ${species} reference files from ${referenceBase}" >> ${repRID}.getRef.log
    if [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references" ]
    then
      echo -e "LOG: grabbing reference files from local (BioHPC)" >> ${repRID}.getRef.log
      ln -s "\${references}"/hisat2
      ln -s "\${references}"/bed
      ln -s "\${references}"/genome.fna
      ln -s "\${references}"/genome.gtf
      ln -s "\${references}"/geneID.tsv
      ln -s "\${references}"/Entrez.tsv
    #elif [ ${referenceBase} == "s3://bicf-references" ]
    #then
    #  echo -e "LOG: grabbing reference files from S3" >> ${repRID}.getRef.log
    #  aws s3 cp "\${references}"/hisat2 ./hisat2 --recursive
    #  aws s3 cp "\${references}"/bed ./bed --recursive
    #  aws s3 cp "\${references}"/genome.fna ./
    #  aws s3 cp "\${references}"/genome.gtf ./
    #  aws s3 cp "\${references}"/geneID.tsv ./
    #  aws s3 cp "\${references}"/Entrez.tsv ./
    elif [ ${referenceBase} == "dev.gudmap.org" ]
    then
      GRCv=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f1)
      GRCp=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f2)
      GENCODE=\$(echo \${references} | grep -o ${refName}.* | cut -d '.' -f3)
      query=\$(echo 'https://${referenceBase}/ermrest/catalog/2/entity/RNASeq:Reference_Genome/Reference_Version='\${GRCv}'.'\${GRCp}'/Annotation_Version=GENCODE%20'\${GENCODE})
      curl --request GET \${query} > refQuery.json
      refURL=\$(python extractRefData.py --returnParam URL)
      loc=\$(dirname \${refURL})
      fName=\$(python extractRefData.py --returnParam fName)
      fName=\${fName%.*}
      if [ "\${loc}" = "/hatrac/*" ]; then echo "LOG: Reference not present in hatrac"; exit 1; fi
      filename=\$(echo \$(basename \${refURL}) | grep -oP '.*(?=:)')
      deriva-hatrac-cli --host ${referenceBase} get \${refURL}
      unzip \$(basename \${refURL})
      mv \${fName}/* .
      mv bed ${refName}/
    fi
    echo -e "LOG: fetched" >> ${repRID}.getRef.log
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
    val ends from endsInfer_alignData
    val stranded from strandedInfer_alignData
    path fastq from fastqsTrim_alignData
    path reference_alignData

  output:
    tuple path ("${repRID}.sorted.bam"), path ("${repRID}.sorted.bam.bai") into rawBam
    path ("*.alignSummary.txt") into alignQC

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
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.bam ${repRID}.bam

    # index the sorted bam using Samtools
    echo -e "LOG: indexing sorted bam file" >> ${repRID}.align.log
    samtools index -@ `nproc` -b ${repRID}.sorted.bam ${repRID}.sorted.bam.bai
    """
}

// Replicate rawBam for multiple process inputs
rawBam.into {
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

  output:
    tuple path ("${repRID}.sorted.deduped.bam"), path ("${repRID}.sorted.deduped.bam.bai") into dedupBam
    tuple path ("${repRID}.sorted.deduped.*.bam"), path ("${repRID}.sorted.deduped.*.bam.bai") into dedupChrBam 
    path ("*.deduped.Metrics.txt") into dedupQC

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
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.deduped.bam ${repRID}.deduped.bam
    
    # index the sorted bam using Samtools
    echo -e "LOG: indexing sorted bam file" >> ${repRID}.dedup.log
    samtools index -@ `nproc` -b ${repRID}.sorted.deduped.bam ${repRID}.sorted.deduped.bam.bai

    # split the deduped BAM file for multi-threaded tin calculation
    for i in `samtools view ${repRID}.sorted.deduped.bam | cut -f3 | sort | uniq`;
      do
      echo "echo \"LOG: splitting each chromosome into its own BAM and BAI files with Samtools\"; samtools view -b ${repRID}.sorted.deduped.bam \${i} 1>> ${repRID}.sorted.deduped.\${i}.bam; samtools index -@ `nproc` -b ${repRID}.sorted.deduped.\${i}.bam ${repRID}.sorted.deduped.\${i}.bam.bai"
    done | parallel -j `nproc` -k
    """
}

// Replicate dedup bam/bai for multiple process inputs
dedupBam.into {
  dedupBam_countData
  dedupBam_makeBigWig
  dedupBam_dataQC
}

/*
 *makeBigWig: make BigWig files for output
*/
process makeBigWig {
  tag "${repRID}"
  publishDir "${outDir}/bigwig", mode: 'copy', pattern: "${repRID}.bw"

  input:
    tuple path (bam), path (bai) from dedupBam_makeBigWig

  output:
    path ("${repRID}.bw")

  script:
    """
    hostname > ${repRID}.makeBigWig.log
    ulimit -a >> ${repRID}.makeBigWig.log

    # create bigwig
    echo -e "LOG: creating bibWig" >> ${repRID}.makeBigWig.log
    bamCoverage -p `nproc` -b ${bam} -o ${repRID}.bw
    echo -e "LOG: created" >> ${repRID}.makeBigWig.log
    """
}

/*
 *countData: count data and calculate tpm
*/
process countData {
  tag "${repRID}"
  publishDir "${outDir}/count", mode: 'copy', pattern: "${repRID}*.tpmTable.csv"

  input:
    path script_calculateTPM
    path script_convertGeneSymbols
    tuple path (bam), path (bai) from dedupBam_countData
    path ref from reference_countData
    val ends from endsInfer_countData
    val stranded from strandedInfer_countData

  output:
    path ("*.tpmTable.csv") into counts
    path ("*.countData.summary") into countsQC
    path ("assignedReads.csv") into inferMetadata_assignedReads

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
      featureCounts -T `nproc` -a ./genome.gtf -G ./genome.fna -g 'gene_name' --extraAttributes 'gene_id' -o ${repRID}.countData -s \${stranding} -R SAM --primary --ignoreDup ${repRID}.sorted.deduped.bam
    elif [ "${ends}" == "pe" ]
    then
      featureCounts -T `nproc` -a ./genome.gtf -G ./genome.fna -g 'gene_name' --extraAttributes 'gene_id' -o ${repRID}.countData -s \${stranding} -p -B -R SAM --primary --ignoreDup ${repRID}.sorted.deduped.bam
    fi
    echo -e "LOG: counted" >> ${repRID}.countData.log
    
    # extract assigned reads
    grep -m 1 'Assigned' *.countData.summary | grep -oe '\\([0-9.]*\\)' > assignedReads.csv

    # calculate TPM from the resulting countData table
    echo -e "LOG: calculating TPM with R" >> ${repRID}.countData.log
    Rscript calculateTPM.R --count "${repRID}.countData"

    # convert gene symbols to Entrez id's
    echo -e "LOG: convert gene symbols to Entrez id's" >> ${repRID}.countData.log
    Rscript convertGeneSymbols.R --repRID "${repRID}"
    """
}

// Extract number of assigned reads metadata into channel
assignedReadsInfer = Channel.create()
inferMetadata_assignedReads.splitCsv(sep: ",", header: false).separate(
  assignedReadsInfer
)

/*
 *fastqc: run fastqc on untrimmed fastq's
*/
process fastqc {
  tag "${repRID}"

  input:
    path (fastq) from fastqs_fastqc

  output:
    path ("*_fastqc.zip") into fastqc
    path ("rawReads.csv") into inferMetadata_rawReads

  script:
    """
    hostname > ${repRID}.fastqc.log
    ulimit -a >> ${repRID}.fastqc.log

    # run fastqc
    echo -e "LOG: running fastq on raw fastqs" >> ${repRID}.fastqc.log
    fastqc *.fastq.gz -o .

    # count raw reads
    zcat *.R1.fastq.gz | echo \$((`wc -l`/4)) > rawReads.csv
    """
}

// Extract number of raw reads metadata into channel
rawReadsInfer = Channel.create()
inferMetadata_rawReads.splitCsv(sep: ",", header: false).separate(
  rawReadsInfer
)

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
    
  output:
    path "${repRID}.tin.hist.tsv" into tinHist
    path "${repRID}.tin.med.csv" into inferMetadata_tinMed
    path "${repRID}.insertSize.inner_distance_freq.txt" into innerDistance
  
  script:
    """
    hostname > ${repRID}.dataQC.log
    ulimit -a >> ${repRID}.dataQC.log

    # calcualte TIN values per feature on each chromosome
    echo -e  "geneID\tchrom\ttx_start\ttx_end\tTIN" > ${repRID}.sorted.deduped.tin.xls
    for i in `cat ./bed/genome.bed | cut -f1 | sort | uniq`; do
      echo "echo \"LOG: running tin.py on \${i}\" >> ${repRID}.dataQC.log; tin.py -i ${repRID}.sorted.deduped.\${i}.bam  -r ./bed/genome.bed; cat ${repRID}.sorted.deduped.\${i}.tin.xls | tr -s \"\\w\" \"\\t\" | grep -P \\\"\\\\t\${i}\\\\t\\\";";
    done | parallel -j `nproc` -k 1>> ${repRID}.sorted.deduped.tin.xls

    # bin TIN values
    echo -e "LOG: binning TINs" >> ${repRID}.dataQC.log
    python3 ${script_tinHist} -r ${repRID}
    echo -e "LOG: binned" >> ${repRID}.dataQC.log

    # calculate inner-distances for PE data
    if [ "${ends}" == "pe" ]
    then
      echo -e "LOG: calculating inner distances for ${ends}" >> ${repRID}.dataQC.log
      inner_distance.py -i "${bam}" -o ${repRID}.insertSize -r ./bed/genome.bed
      echo -e "LOG: calculated" >> ${repRID}.dataQC.log
    elif [ "${ends}" == "se" ]
    then
      echo -e "LOG: creating dummy inner distance file for ${ends}" >> ${repRID}.dataQC.log
      touch ${repRID}.insertSize.inner_distance_freq.txt
    fi
    """
}

// Extract median TIN metadata into channel
tinMedInfer = Channel.create()
inferMetadata_tinMed.splitCsv(sep: ",", header: false).separate(
  tinMedInfer
)

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
    val endsM from endsMeta
    val strandedM from strandedMeta
    val spikeM from spikeMeta
    val speciesM from speciesMeta
    val endsI from endsInfer_aggrQC
    val strandedI from strandedInfer_aggrQC
    val spikeI from spikeInfer_aggrQC
    val speciesI from speciesInfer_aggrQC
    val readLengthM from readLengthMeta
    val readLengthI from readLengthInfer
    val rawReadsI from rawReadsInfer
    val assignedReadsI from assignedReadsInfer
    val tinMedI from tinMedInfer
    val expRID
    val studyRID

  output:
    path "${repRID}.multiqc.html" into multiqc
    path "${repRID}.multiqc_data.json" into multiqcJSON

  script:
    """
    hostname > ${repRID}.aggrQC.log
    ulimit -a >> ${repRID}.aggrQC.log

    # make run table
    if [ "${params.inputBagForce}" == "" ] && [ "${params.fastqsForce}" == "" ] && [ "${params.speciesForce}" == "" ]
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
      echo -e "Infered\t${speciesI}\t${endsI}\t${strandedI}\t${spikeI}\t-\t-\t-\t-" >> metadata.tsv
    else
      echo -e "Infered\t${speciesI} (FORCED)\t${endsI}\t${strandedI}\t${spikeI}\t-\t-\t-\t-" >> metadata.tsv
    fi
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
    """
}

/*
 *ouputBag: create ouputBag
*/
process outputBag {
  tag "${repRID}"
  publishDir "${outDir}/outputBag", mode: 'copy', pattern: "Replicate_${repRID}.outputBag.zip"
  
  input:
    path multiqc
    path multiqcJSON
  
  output:
    path ("Replicate_*.zip") into outputBag

  script:
  """
  mkdir Replicate_${repRID}.outputBag
  cp ${multiqc} Replicate_${repRID}.outputBag
  cp ${multiqcJSON} Replicate_${repRID}.outputBag
  bdbag Replicate_${repRID}.outputBag --archiver zip
  """
}