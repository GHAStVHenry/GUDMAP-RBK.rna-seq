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
params.refMoVersion = "38.p6.vM22"
params.refHuVersion = "38.p12.v31"
params.refERCCVersion = "92"
params.outDir = "${baseDir}/../output"

// Parse input variables
deriva = Channel
  .fromPath(params.deriva)
  .ifEmpty { exit 1, "deriva credential file not found: ${params.deriva}" }
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "deriva cookie file for bdbag not found: ${params.bdbag}" }
repRID = params.repRID
refMoVersion = params.refMoVersion
refHuVersion = params.refHuVersion
refERCCVersion = params.refERCCVersion
outDir = params.outDir
logsDir = "${outDir}/Logs"

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")
//referenceBase = "s3://bicf-references"
referenceBase = "/project/BICF/BICF_Core/shared/gudmap/references"
referenceInfer = Channel.fromList(["ERCC","GRCh","GRCm"])
multiqcConfig = Channel.fromPath("${baseDir}/conf/multiqc_config.yaml")

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbagFetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parseMeta.py")
script_inferMeta = Channel.fromPath("${baseDir}/scripts/inferMeta.sh")
script_calculateTPM = Channel.fromPath("${baseDir}/scripts/calculateTPM.R")
script_tinHist = Channel.fromPath("${baseDir}/scripts/tinHist.py")

/*
 * trackStart: track start of pipeline
 */
params.ci = false
params.dev = false
process trackStart {
  script:
  """
  hostname
  ulimit -a
  export https_proxy=\${http_proxy}
  
  curl -H 'Content-Type: application/json' -X PUT -d \
    '{ \
      "sessionId": "${workflow.sessionId}", \
      "pipeline": "gudmap.rbk_rnaseq", \
      "start": "${workflow.start}", \
      "repRID": "${repRID}", \
      "astrocyte": false, \
      "status": "started", \
      "nextflowVersion": "${workflow.nextflow.version}", \
      "ci": ${params.ci}, \
      "dev": ${params.dev} \
    }' \
    "https://xku43pcwnf.execute-api.us-east-1.amazonaws.com/ProdDeploy/pipeline-tracking"
  """
 }

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
process getBag {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.getBag.{out,err}"

  input:
    path credential, stageAs: "credential.json" from deriva
    path derivaConfig

  output:
    path ("Replicate_*.zip") into bagit
    path ("${repRID}.getBag.err")

  script:
    """
    hostname > ${repRID}.getBag.err
    ulimit -a >> ${repRID}.getBag.err
    export https_proxy=\${http_proxy}

    # link credential file for authentication
    ln -sf `readlink -e credential.json` ~/.deriva/credential.json 1>> ${repRID}.getBag.out 2>> ${repRID}.getBag.err
    echo "LOG: deriva credentials linked" >> ${repRID}.getBag.err

    # deriva-download replicate RID
    echo "LOG: fetching deriva catalog for selected RID in GUDMAP." >> ${repRID}.getBag.err
    deriva-download-cli dev.gudmap.org --catalog 2 ${derivaConfig} . rid=${repRID} 1>> ${repRID}.getBag.out 2>> ${repRID}.getBag.err
    """
}

/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 */
process getData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.getData.{out,err}"

  input:
    path script_bdbagFetch
    path cookies, stageAs: "deriva-cookies.txt" from bdbag
    path bagit

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    path ("**/File.csv") into fileMeta
    path ("**/Experiment Settings.csv") into experimentSettingsMeta
    path ("**/Experiment.csv") into experimentMeta
    path ("${repRID}.getData.{out,err}")

  script:
    """
    hostname > ${repRID}.getData.err
    ulimit -a >> ${repRID}.getData.err
    export https_proxy=\${http_proxy}
    
    # link deriva cookie for authentication
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt 1>> ${repRID}.getData.out 2>> ${repRID}.getData.err
    echo "LOG: deriva cookie linked" >> ${repRID}.getData.err
    
    # get bagit basename
    replicate=\$(basename "${bagit}" | cut -d "." -f1) 1>> ${repRID}.getData.out 2>> ${repRID}.getData.err
    echo "LOG: \${replicate}" >> ${repRID}.getData.err
    
    # unzip bagit
    unzip ${bagit} 1>> ${repRID}.getData.out 2>> ${repRID}.getData.err
    echo "LOG: replicate bdbag unzipped" >> ${repRID}.getData.err
    
    # bagit fetch fastq"s only and rename by repRID
    sh ${script_bdbagFetch} \${replicate} ${repRID} 1>> ${repRID}.getData.out 2>> ${repRID}.getData.err
    echo "LOG: replicate bdbag fetched" >> ${repRID}.getData.err
    """
}

// Replicate raw fastq's for multiple process inputs
fastqs.into {
  fastqs_trimData
  fastqs_fastqc
}

/*
 * parseMetadata: parses metadata to extract experiment parameters
*/
process parseMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.parseMetadata.{out,err}"

  input:
    path script_parseMeta
    path fileMeta
    path experimentSettingsMeta
    path experimentMeta

  output:
    path "design.csv" into metadata
    path "${repRID}.parseMetadata.{out,err}"

  script:
    """
    hostname > ${repRID}.parseMetadata.err
    ulimit -a >> ${repRID}.parseMetadata.err

    # Check replicate RID metadata
    rep=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p repRID)
    echo "LOG: replicate RID metadata parsed: \${rep}" >> ${repRID}.parseMetadata.err

    # Get experiment RID metadata
    exp=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p expRID)
    echo "LOG: experiment RID metadata parsed: \${exp}" >> ${repRID}.parseMetadata.err
    
    # Get study RID metadata
    study=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p studyRID)
    echo "LOG: study RID metadata parsed: \${study}" >> ${repRID}.parseMetadata.err

    # Get endedness metadata
    endsMeta=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p endsMeta)
    echo "LOG: endedness metadata parsed: \${endsMeta}" >> ${repRID}.parseMetadata.err
    
    # Manually get endness
    endsManual=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p endsManual)
    echo "LOG: endedness manually detected: \${endsManual}" >> ${repRID}.parseMetadata.err

    # Get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p stranded)
    echo "LOG: strandedness metadata parsed: \${stranded}" >> ${repRID}.parseMetadata.err
    
    # Get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p spike)
    echo "LOG: spike-in metadata parsed: \${spike}" >> ${repRID}.parseMetadata.err
    
    # Get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentMeta}" -p species)
    echo "LOG: species metadata parsed: \${species}" >> ${repRID}.parseMetadata.err

    # Save design file
    echo "\${endsMeta},\${endsManual},\${stranded},\${spike},\${species},\${exp},\${study}" > design.csv
    """
}

// Split metadata into separate channels
endsMeta = Channel.create()
endsManual = Channel.create()
strandedMeta = Channel.create()
spikeMeta = Channel.create()
speciesMeta = Channel.create()
expRID = Channel.create()
studyRID = Channel.create()
metadata.splitCsv(sep: ",", header: false).separate(
  endsMeta,
  endsManual,
  strandedMeta,
  spikeMeta,
  speciesMeta,
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
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.trimData.{out,err}"

  input:
    val ends from endsManual_trimData
    path (fastq) from fastqs_trimData

  output:
    path ("*.fq.gz") into fastqsTrim
    path ("*_trimming_report.txt") into trimQC
    path ("${repRID}.trimData.{out,err}")

  script:
    """
    hostname > ${repRID}.trimData.err
    ulimit -a >> ${repRID}.trimData.err

    #Trim fastq's using trim_galore
    if [ "${ends}" == "se" ]
    then
      echo "LOG: running trim_galore using single-end settings" >> ${repRID}.trimData.err
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>> ${repRID}.trimData.out 2>> ${repRID}.trimData.err
    elif [ "${ends}" == "pe" ]
    then
      echo "LOG: running trim_galore using paired-end settings" >> ${repRID}.trimData.err
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>> ${repRID}.trimData.out 2>> ${repRID}.trimData.err
    fi
    """
}

// Replicate trimmed fastq's
fastqsTrim.into {
  fastqsTrim_alignData
  fastqsTrim_downsampleData
}

/*
  * getRefInfer: Dowloads appropriate reference for metadata inference
*/
process getRefInfer {
  tag "${refName}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.getRefInfer.{out,err}"

  input:
    val refName from referenceInfer

  output:
    tuple val (refName), path ("hisat2", type: 'dir'), path ("*.fna"), path ("*.gtf")  into refInfer
    path ("${refName}", type: 'dir') into bedInfer
    path ("${repRID}.getRefInfer.{out,err}")
 
  script:
    """
    hostname > ${repRID}.getRefInfer.err
    ulimit -a >> ${repRID}.getRefInfer.err
    export https_proxy=\${http_proxy}

    #Set the reference name
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
      echo -e "LOG: ERROR - References could not be set!\nReference found: ${referenceBase}" >> ${repRID}.getRefInfer.err
      exit 1
    fi
    mkdir ${refName}
    #Retreive appropriate reference appropriate location
    if [ ${referenceBase} == "s3://bicf-references" ]
    then
      echo "LOG: grabbing reference files from S3" >> ${repRID}.getRefInfer.err
      aws s3 cp "\${references}" /hisat2 ./ --recursive 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      aws s3 cp "\${references}" /bed ./${refName}/ --recursive 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      aws s3 cp "\${references}" /*.fna --recursive 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      aws s3 cp "\${references}" /*.gtf --recursive 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
    elif [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references" ]
    then
      echo "LOG: using pre-defined locations for reference files" >> ${repRID}.getRefInfer.err
      ln -s "\${references}"/hisat2 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      ln -s "\${references}"/bed ${refName}/bed 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      ln -s "\${references}"/genome.fna 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
      ln -s "\${references}"/genome.gtf 1>> ${repRID}.getRefInfer.out 2>> ${repRID}.getRefInfer.err
    fi

    #Make blank bed folder for ERCC
    if [ "${refName}" == "ERCC" ]
    then
      rm ${refName}/bed
      mkdir ${refName}/bed
    fi
    """
}

/*
 * downsampleData: downsample fastq's for metadata inference
 */
process downsampleData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.downsampleData.{out,err}"

  input:
    val ends from endsManual_downsampleData
    path fastq from fastqsTrim_downsampleData

  output:
    path ("sampled.1.fq") into fastqs1Sample
    path ("sampled.2.fq") into fastqs2Sample
    path ("${repRID}.downsampleData.{out,err}")

  script:
    """
    hostname > ${repRID}.downsampleData.err
    ulimit -a >> ${repRID}.downsampleData.err
    export https_proxy=\${http_proxy}

    if [ "${ends}" == "se" ]
    then
      echo "LOG: downsampling single-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *trimmed.fq.gz 100000 1> sampled.1.fq 2>> ${repRID}.downsampleData.err
      touch sampled.2.fq
    elif [ "${ends}" == "pe" ]
    then
      echo "LOG: downsampling read 1 of paired-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *1.fq.gz 1000000 1> sampled.1.fq 2>> ${repRID}.downsampleData.err
      echo "LOG: downsampling read 2 of paired-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *2.fq.gz 1000000 1> sampled.2.fq 2>> ${repRID}.downsampleData.err
    fi
    """
}

// Replicate the dowsampled fastq's and attatched to the references
inferInput = endsManual_alignSampleData.combine(refInfer.combine(fastqs1Sample.collect().combine(fastqs2Sample.collect())))

/*
 * alignSampleData: aligns the downsampled reads to a reference database
*/
process alignSampleData {
  tag "${ref}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.alignSampleData.{out,err}"

  input:
    tuple val (ends), val (ref), path (hisat2), path (fna), path (gtf), path (fastq1), path (fastq2) from inferInput

  output:
    path ("${ref}.sampled.sorted.bam") into sampleBam
    path ("${ref}.sampled.sorted.bam.bai") into sampleBai
    path ("${ref}.alignSampleSummary.txt") into alignSampleQC
    path ("${repRID}.${ref}.alignSampleData.{out,err}")

  script:
    """
    hostname > ${repRID}.${ref}.alignSampleData.err
    ulimit -a >> ${repRID}.${ref}.alignSampleData.err

    #Align the reads with Hisat 2
    if [ "${ends}" == "se" ]
    then
      echo "LOG: running Hisat2 with single-end settings" >> ${repRID}.${ref}.alignSampleData.err
      hisat2 -p `nproc` --add-chrname -S ${ref}.sampled.sam -x hisat2/genome -U ${fastq1} --summary-file ${ref}.alignSampleSummary.txt --new-summary 1>> ${repRID}.${ref}.alignSampleData.out 2>> ${repRID}.${ref}.alignSampleData.err
    elif [ "${ends}" == "pe" ]
    then
      echo "LOG: running Hisat2 with paired-end settings" >> ${repRID}.${ref}.alignSampleData.err
      hisat2 -p `nproc` --add-chrname -S ${ref}.sampled.sam -x hisat2/genome --no-mixed --no-discordant -1 ${fastq1} -2 ${fastq2} --summary-file ${ref}.alignSampleSummary.txt --new-summary 1>> ${repRID}.${ref}.alignSampleData.out 2>> ${repRID}.${ref}.alignSampleData.err
    fi
    
    #Convert the output sam file to a sorted bam file using Samtools
    echo "LOG: converting from sam to bam" >> ${repRID}.${ref}.alignSampleData.err
    samtools view -1 -@ `nproc` -F 4 -F 8 -F 256 -o ${ref}.sampled.bam ${ref}.sampled.sam 1>> ${repRID}.${ref}.alignSampleData.out 2>> ${repRID}.${ref}.alignSampleData.err;

    #Sort the bam file using Samtools
    echo "LOG: sorting the bam file" >> ${repRID}.${ref}.alignSampleData.err
    samtools sort -@ `nproc` -O BAM -o ${ref}.sampled.sorted.bam ${ref}.sampled.bam 1>> ${repRID}.${ref}.alignSampleData.out 2>> ${repRID}.${ref}.alignSampleData.err;

    #Index the sorted bam using Samtools
    echo "LOG: indexing sorted bam file" >> ${repRID}.${ref}.alignSampleData.err
    samtools index -@ `nproc` -b ${ref}.sampled.sorted.bam ${ref}.sampled.sorted.bam.bai 1>> ${repRID}.${ref}.alignSampleData.out 2>> ${repRID}.${ref}.alignSampleData.err;
    """
}

alignSampleQC.into {
  alignSampleQC_inferMetadata
  alignSampleQC_aggrQC
}

process inferMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.inferMetadata.{out,err}"

  input:
    path script_inferMeta
    path beds from bedInfer.collect()
    path bam from sampleBam.collect()
    path bai from sampleBai.collect()
    path alignSummary from alignSampleQC_inferMetadata.collect()

  output:
    path "infer.csv" into inferMetadata
    path "${repRID}.infer_experiment.txt" into inferExperiment
    path "${repRID}.inferMetadata.{out,err}" optional true

  script:
    """
    hostname > ${repRID}.inferMetadata.err
    ulimit -a >> ${repRID}.inferMetadata.err

    # collect alignment rates (round down to integers)
    align_ercc=\$(echo \$(grep "Overall alignment rate" ERCC.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_ercc=\$(echo \${align_ercc%.*})
    align_hu=\$(echo \$(grep "Overall alignment rate" GRCh.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_hu=\$(echo \${align_hu%.*})
    align_mo=\$(echo \$(grep "Overall alignment rate" GRCm.alignSampleSummary.txt | cut -f2 -d ':' | cut -f2 -d ' ' | tr -d '%'))
    align_mo=\$(echo \${align_mo%.*})

    # determine spike-in
    if [ 1 -eq \$(echo \$(expr \${align_ercc} ">=" 10)) ]
    then
      spike="yes"
    else
      spike="no"
    fi
    echo -e "LOG: Inference of strandedness results is: \${spike}" >> ${repRID}.inferMetadata.err

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
      echo -e "LOG: ERROR - Inference of species returns an ambiguous result: hu=\${align_hu} mo=\${align_mo}" >> ${repRID}.inferMetadata.err
      exit 1
    fi
    echo -e "LOG: Inference of species results in: \${species}" >> ${repRID}.inferMetadata.err

    # infer experimental setting from dedup bam
    echo "LOG: infer experimental setting from dedup bam" >> ${repRID}.inferMetadata.err
    infer_experiment.py -r "\${bed}" -i "\${bam}" 1>> ${repRID}.infer_experiment.txt 2>> ${repRID}.inferMetadata.err

    echo "LOG: determining endedness and strandedness from file" >> ${repRID}.inferMetadata.err
    ended=`bash inferMeta.sh endness ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
    fail=`bash inferMeta.sh fail ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
    if [ \${ended} == "PairEnd" ] 
    then
      ends="pe"
      percentF=`bash inferMeta.sh pef ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
      percentR=`bash inferMeta.sh per ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
    elif [ \${ended} == "SingleEnd" ]
    then
      ends="se"
      percentF=`bash inferMeta.sh sef ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
      percentR=`bash inferMeta.sh ser ${repRID}.infer_experiment.txt` 1>> ${repRID}.inferMetadata.out 2>> ${repRID}.inferMetadata.err
    fi
    if [ 1 -eq \$(echo \$(expr \${percentF#*.} ">" 2500)) ] && [ 1 -eq \$(echo \$(expr \${percentR#*.} "<" 2500)) ]
    then
      stranded="forward"
    elif [ 1 -eq \$(echo \$(expr \${percentR#*.} ">" 2500)) ] && [ 1 -eq \$(echo \$(expr \${percentF#*.} "<" 2500)) ]
    then
      stranded="reverse"

    else
      stranded="unstranded"
    fi
    echo -e "LOG: stradedness set to \${stranded}" >> ${repRID}.inferMetadata.err

    # write infered metadata to file
    echo "\${ends},\${stranded},\${spike},\${species},\${align_ercc},\${align_hu},\${align_mo},\${percentF},\${percentR},\${fail}" 1>> infer.csv 2>> ${repRID}.inferMetadata.err
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
  * getRef: Dowloads appropriate reference
*/
process getRef {
  tag "${species}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.getRef.{out,err}"

  input:
    val spike from spikeInfer_getRef
    val species from speciesInfer_getRef

  output:
    tuple path ("hisat2", type: 'dir'), path ("bed", type: 'dir'), path ("*.fna"), path ("*.gtf")  into reference
    path ("${repRID}.getRef.{out,err}")
 
  script:
    """
    hostname > ${repRID}.getRef.err
    ulimit -a >> ${repRID}.getRef.err
    export https_proxy=\${http_proxy}

    #Set the reference name
    if [ "${species}" == "Mus musculus" ]
    then
      references=\$(echo ${referenceBase}/GRCm${refMoVersion})
    elif [ '${species}' == "Homo sapiens" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
    else
      echo -e "LOG: ERROR - References could not be set!\nSpecies reference found: ${species}" >> ${repRID}.getRef.err
      exit 1
    fi
    if [ "${spike}" == "yes" ]
    then
      references=\$(echo \${reference}-S/)
    elif [ "${spike}" == "no" ]
    then
      reference=\$(echo \${references}/)
    fi
    echo "LOG: species set to \${references}" >> ${repRID}.getRef.err

    #Retreive appropriate reference appropriate location
    if [ ${referenceBase} == "s3://bicf-references" ]
    then
      echo "LOG: grabbing reference files from S3" >> ${repRID}.getRef.err
      aws s3 cp "\${references}" /hisat2 ./ --recursive 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      aws s3 cp "\${references}" /bed ./ --recursive 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      aws s3 cp "\${references}" /*.fna --recursive 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      aws s3 cp "\${references}" /*.gtf --recursive 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
    elif [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references" ]
    then
      echo "LOG: using pre-defined locations for reference files" >> ${repRID}.getRef.err
      ln -s "\${references}"/hisat2 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      ln -s "\${references}"/bed 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      ln -s "\${references}"/genome.fna 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
      ln -s "\${references}"/genome.gtf 1>> ${repRID}.getRef.out 2>> ${repRID}.getRef.err
    fi
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
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.align.{out,err}"

  input:
    val ends from endsInfer_alignData
    val stranded from strandedInfer_alignData
    path fastq from fastqsTrim_alignData
    path reference_alignData

  output:
    tuple path ("${repRID}.sorted.bam"), path ("${repRID}.sorted.bam.bai") into rawBam
    path ("*.alignSummary.txt") into alignQC
    path ("${repRID}.align.{out,err}")

  script:
    """
    hostname > ${repRID}.align.err
    ulimit -a >> ${repRID}.align.err

    #Set stranded param for hisat2
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

    #Align the reads with Hisat 2
    if [ "${ends}" == "se" ]
    then
      echo "LOG: running Hisat2 with single-end settings" >> ${repRID}.align.err
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome \${strandedParam} -U ${fastq[0]} --summary-file ${repRID}.alignSummary.txt --new-summary 1>> ${repRID}.align.out 2>> ${repRID}.align.err
    elif [ "${ends}" == "pe" ]
    then
      echo "LOG: running Hisat2 with paired-end settings" >> ${repRID}.align.err
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome \${strandedParam} --no-mixed --no-discordant -1 ${fastq[0]} -2 ${fastq[1]} --summary-file ${repRID}.alignSummary.txt --new-summary 1>> ${repRID}.align.out 2>> ${repRID}.align.err
    fi
    
    #Convert the output sam file to a sorted bam file using Samtools
    echo "LOG: converting from sam to bam" >> ${repRID}.align.err
    samtools view -1 -@ `nproc` -F 4 -F 8 -F 256 -o ${repRID}.bam ${repRID}.sam 1>> ${repRID}.align.out 2>> ${repRID}.align.err;

    #Sort the bam file using Samtools
    echo "LOG: sorting the bam file" >> ${repRID}.align.err
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.bam ${repRID}.bam 1>> ${repRID}.align.out 2>> ${repRID}.align.err;

    #Index the sorted bam using Samtools
    echo "LOG: indexing sorted bam file" >> ${repRID}.align.err
    samtools index -@ `nproc` -b ${repRID}.sorted.bam ${repRID}.sorted.bam.bai 1>> ${repRID}.align.out 2>> ${repRID}.align.err;
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
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.dedup.{out,err}"

  input:
    tuple path (bam), path (bai) from rawBam_dedupData

  output:
    tuple path ("${repRID}.sorted.deduped.bam"), path ("${repRID}.sorted.deduped.bam.bai") into dedupBam
    tuple path ("${repRID}.sorted.deduped.*.bam"), path ("${repRID}.sorted.deduped.*.bam.bai") into dedupChrBam 
    path ("*.deduped.Metrics.txt") into dedupQC
    path ("${repRID}.dedup.{out,err}")

  script:
    """
    hostname > ${repRID}.dedup.err
    ulimit -a >> ${repRID}.dedup.err

    # remove duplicated reads using Picard's MarkDuplicates
    echo "LOG: running picard MarkDuplicates to remove duplicate reads" >> ${repRID}.dedup.err
    java -jar /picard/build/libs/picard.jar MarkDuplicates I=${bam} O=${repRID}.deduped.bam M=${repRID}.deduped.Metrics.txt REMOVE_DUPLICATES=true 1>> ${repRID}.dedup.out 2>> ${repRID}.dedup.err

    # Sort the bam file using Samtools
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.deduped.bam ${repRID}.deduped.bam 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    
    # Index the sorted bam using Samtools
    samtools index -@ `nproc` -b ${repRID}.sorted.deduped.bam ${repRID}.sorted.deduped.bam.bai 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    
    # Split the deduped BAM file for multi-threaded tin calculation
    for i in `samtools view ${repRID}.sorted.deduped.bam | cut -f3 | sort | uniq`;
      do
      echo "echo \"LOG: splitting each chromosome into its own BAM and BAI files with Samtools\" >> ${repRID}.dedup.err; samtools view -b ${repRID}.sorted.deduped.bam \${i} > ${repRID}.sorted.deduped.\${i}.bam; samtools index -@ `nproc` -b ${repRID}.sorted.deduped.\${i}.bam ${repRID}.sorted.deduped.\${i}.bam.bai"
    done | parallel -j `nproc` -k 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    """
}

// Replicate dedup bam/bai for multiple process inputs
dedupBam.into {
  dedupBam_countData
  dedupBam_makeBigWig
  dedupBam_dataQC
}

/*
 *Make BigWig files for output
*/
process makeBigWig {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.makeBigWig.{out,err}"
  publishDir "${outDir}/bigwig", mode: 'copy', pattern: "${repRID}.bw"

  input:
    tuple path (bam), path (bai) from dedupBam_makeBigWig

  output:
    path ("${repRID}.bw")
    path ("${repRID}.makeBigWig.{out,err}")

  script:
    """
    hostname > ${repRID}.makeBigWig.err
    ulimit -a >> ${repRID}.makeBigWig.err

    #Run bamCoverage
    echo "LOG: Running bigWig bamCoverage" >> ${repRID}.makeBigWig.err
    bamCoverage -p `nproc` -b ${bam} -o ${repRID}.bw 1>> ${repRID}.makeBigWig.out 2>> ${repRID}.makeBigWig.err
    """
}

/*
 *Run countData and get the counts, tpm
*/
process countData {
  tag "${repRID}"
  publishDir "${outDir}/countData", mode: 'copy', pattern: "${repRID}*.countTable.csv"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.countData.{out,err}"

  input:
    path script_calculateTPM
    tuple path (bam), path (bai) from dedupBam_countData
    path ref from reference_countData
    val ends from endsInfer_countData
    val stranded from strandedInfer_countData

  output:
    path ("*.countTable.csv") into counts
    path ("*.countData.summary") into countsQC
    path ("${repRID}.countData.{out,err}")

  script:
    """
    hostname > ${repRID}.countData.err
    ulimit -a >> ${repRID}.countData.err

    #Determine strandedness and setup strandig for countData
    stranding=0;
    if [ "${stranded}" == "unstranded" ]
    then
      stranding=0
      echo "LOG: strandedness set to unstranded [0]" >> ${repRID}.countData.err
    elif [ "${stranded}" == "forward" ]
    then
      stranding=1
      echo "LOG: strandedness set to forward stranded [1]" >> ${repRID}.countData.err
    elif [ "${stranded}" == "reverse" ]
    then
      stranding=2
      echo "LOG: strandedness set to forward stranded [2]" >> ${repRID}.countData.err
    fi
    #Run countData
    echo "LOG: running countData on the data" >> ${repRID}.countData.err
    if [ "${ends}" == "se" ]
    then
      featureCounts -R SAM -p -G ./genome.fna -T `nproc` -s \${stranding} -a ./genome.gtf -o ${repRID}.countData -g 'gene_name' --primary --ignoreDup ${repRID}.sorted.deduped.bam 1>> ${repRID}.countData.out 2>> ${repRID}.countData.err
    elif [ "${ends}" == "pe" ]
    then
      featureCounts -R SAM -p -G ./genmome.fna -T `nproc` -s \${stranding} -a ./genome.gtf -o ${repRID}.countData -g 'gene_name' --primary --ignoreDup -B ${repRID}.sorted.deduped.bam 1>> ${repRID}.countData.out 2>> ${repRID}.countData.err
    fi

    #Calculate TPM from the resulting countData table
    echo "LOG: calculating TPM with R" >> ${repRID}.countData.err
    Rscript calculateTPM.R --count "${repRID}.countData" 1>> ${repRID}.countData.out 2>> ${repRID}.countData.err
    """
}

/*
 *fastqc: run fastqc on untrimmed fastq's
*/
process fastqc {
  tag "${repRID}"
  publishDir "${outDir}/fastqc", mode: 'copy', pattern: "*_fastqc.zip"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.fastqc.{out,err}"

  input:
    path (fastq) from fastqs_fastqc

  output:
    path ("*_fastqc.zip") into fastqc
    path ("${repRID}.fastqc.{out,err}")

  script:
    """
    hostname > ${repRID}.fastqc.err
    ulimit -a >> ${repRID}.fastqc.err

    # run fastqc
    echo "LOG: beginning FastQC analysis of the data" >> ${repRID}.fastqc.err
    fastqc *.fastq.gz -o . 1>> ${repRID}.fastqc.out 2>> ${repRID}.fastqc.err
    """
}

/*
 *dataQC: run RSeQC to calculate transcript integrity numbers (TIN)
*/
process dataQC {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.dataQC.{out,err}"

  input:
    path script_tinHist
    path ref from reference_dataQC
    tuple path (bam), path (bai) from dedupBam_dataQC
    tuple path (chrBam), path (chrBai) from dedupChrBam
    val ends from endsInfer_dataQC
    
  output:
    path "${repRID}.tin.hist.tsv" into tin
    path "${repRID}.insertSize.inner_distance_freq.txt" into innerDistance
    path "${repRID}.dataQC.{out,err}"

  script:
    """
    hostname > ${repRID}.dataQC.err
    ulimit -a >> ${repRID}.dataQC.err

    # calcualte TIN values per feature on each chromosome
    echo -e  "geneID\tchrom\ttx_start\ttx_end\tTIN" > ${repRID}.sorted.deduped.tin.xls
    for i in `cat ./bed/genome.bed | cut -f1 | sort | uniq`; do
      echo "echo \"LOG: running tin.py on \${i}\" >> ${repRID}.dataQC.err; tin.py -i ${repRID}.sorted.deduped.\${i}.bam  -r ./bed/genome.bed 1>>${repRID}.dataQC.log 2>>${repRID}.dataQC.err; cat ${repRID}.sorted.deduped.\${i}.tin.xls | tr -s \"\\w\" \"\\t\" | grep -P \\\"\\\\t\${i}\\\\t\\\";";
    done | parallel -j `nproc` -k 1>> ${repRID}.sorted.deduped.tin.xls 2>>${repRID}.dataQC.err

    # bin TIN values
    python3 ${script_tinHist} -r ${repRID}

    # calculate inner-distances for PE data
    if [ "${ends}" == "pe" ]
    then
      inner_distance.py -i "${bam}" -o ${repRID}.insertSize -r ./bed/genome.bed 1>>${repRID}.dataQC.out 2>>${repRID}.dataQC.err
    else
      touch ${repRID}.insertSize.inner_distance_freq.txt
    fi
    """
}

/*
 *aggrQC: aggregate QC from processes as wel as metadata and run MultiQC
*/
process aggrQC {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.aggrQC.{out,err}"

  input:
    path multiqcConfig
    path fastqc
    path trimQC
    path alignQC
    path dedupQC
    path countsQC
    path innerDistance
    path tin
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
    val expRID
    val studyRID

  output:
    path "${repRID}.aggrQC.{out,err}" optional true

  script:
    """
    hostname > ${repRID}.aggrQC.err
    ulimit -a >> ${repRID}.aggrQC.err

    echo -e "Replicate RID\tExperiment RID\tStudy RID" > rid.tsv
    echo -e "${repRID}\t${expRID}\t${studyRID}" >> rid.tsv

    echo -e "Source\tSpecies\tEnds\tStranded\tSpike-in" > metadata.tsv
    echo -e "Infered\t${speciesI}\t${endsI}\t${strandedI}\t${spikeI}" >> metadata.tsv
    echo -e "Submitter\t${speciesM}\t${endsM}\t${strandedM}\t${spikeM}" >> metadata.tsv
    echo -e "Manual\t-\t${endsManual}\t-\t-" >> metadata.tsv

    # remove inner distance report if it is empty (SE repRID)
    if [ wc -l ${innerDistance} | awk '{print\${1}}' -eq 0 ]
    then
      rm -f ${innerDistance}
    fi

    #run MultiQC
    multiqc -c ${multiqcConfig} .
    """
}
