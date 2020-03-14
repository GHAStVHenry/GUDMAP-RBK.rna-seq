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
outDir = params.outDir
logsDir = "${outDir}/Logs"

// Define fixed files
derivaConfig = Channel.fromPath("${baseDir}/conf/replicate_export_config.json")
//referenceBase = "s3://bicf-references"
referenceBase = "/project/BICF/BICF_Core/shared/gudmap/references"

// Define script files
script_bdbagFetch = Channel.fromPath("${baseDir}/scripts/bdbagFetch.sh")
script_parseMeta = Channel.fromPath("${baseDir}/scripts/parseMeta.py")
script_inferMeta = Channel.fromPath("${baseDir}/scripts/inferMeta.sh")
script_aggregateInference = Channel.fromPath("${baseDir}/scripts/aggregateInference.R")

/*
 * splitData: split bdbag files by replicate so fetch can occure in parallel, and rename files to replicate rid
 */
process getBag {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "*.getBag.err"

  input:
    path credential, stageAs: "credential.json" from deriva
    path derivaConfig

  output:
    path ("Replicate_*.zip") into bagit
    path ("${repRID}.getBag.err")

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
  publishDir "${logsDir}", mode: "copy", pattern: "*.getData.err"

  input:
    path script_bdbagFetch
    path cookies, stageAs: "deriva-cookies.txt" from bdbag
    path bagit

  output:
    path ("*.R{1,2}.fastq.gz") into fastqs
    path ("**/File.csv") into fileMeta
    path ("**/Experiment Settings.csv") into experimentSettingsMeta
    path ("**/Experiment.csv") into experimentMeta
    path ("${repRID}.getData.err")

  script:
    """
    hostname >>${repRID}.getData.err
    ulimit -a >>${repRID}.getData.err
    export https_proxy=\${http_proxy}
    
    # link deriva cookie for authentication
    ln -sf `readlink -e deriva-cookies.txt` ~/.bdbag/deriva-cookies.txt >>${repRID}.getData.err
    echo "LOG: deriva cookie linked" >>${repRID}.getData.err
    
    # get bagit basename
    replicate=\$(basename "${bagit}" | cut -d "." -f1)
    echo "LOG: \${replicate}" >>${repRID}.getData.err
    
    # unzip bagit
    unzip ${bagit} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag unzipped" >>${repRID}.getData.err
    
    # bagit fetch fastq"s only and rename by repRID
    sh ${script_bdbagFetch} \${replicate} ${repRID} 2>>${repRID}.getData.err
    echo "LOG: replicate bdbag fetched" >>${repRID}.getData.err
    """
}

// Replicate raw fastqs for multiple process inputs
fastqs.into {
  fastqs_trimData
  fastqs_fastqc
}

/*
 * parseMetadata: parses metadata to extract experiment parameters
*/
process parseMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "*.parseMetadata.err"

  input:
    path script_parseMeta
    val repRID
    path fileMeta
    path experimentSettingsMeta
    path experimentMeta

  output:
    path "design.csv" into metadata

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
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p stranded -e \${endsManual})
    echo "LOG: strandedness metadata parsed: \${stranded}" >>${repRID}.parseMetadata.err
    
    # Get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p spike)
    echo "LOG: spike-in metadata parsed: \${spike}" >>${repRID}.parseMetadata.err
    
    # Get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentMeta}" -p species)
    echo "LOG: species metadata parsed: \${species}" >>${repRID}.parseMetadata.err

    # Save design file
    echo "\${rep},\${endsMeta},\${endsManual},\${stranded},\${spike},\${species}" > design.csv
    """
}

// Split metadata into separate channels
rep = Channel.create()
endsMeta = Channel.create()
endsManual = Channel.create()
stranded = Channel.create()
spike = Channel.create()
species = Channel.create()
metadata.splitCsv(sep: ",", header: false).separate(
  rep,
  endsMeta,
  endsManual,
  stranded,
  spike,
  species
)
// Replicate metadata for multiple process inputs
endsManual.into {
  endsManual_trimData
  endsManual_alignData
}
stranded.into {
  stranded_alignData
}
spike.into{
  spike_getRef
}
species.into {
  species_getRef
}

/*
  * getRef: Dowloads appropriate reference
*/
process getRef {
  tag "${species_getRef}"
  publishDir "${logsDir}", mode: "copy", pattern: "*.getRef.err"

  input:
    val spike_getRef
    val species_getRef

  output:
    path ("*")  into reference
  
  script:
    """
    hostname >>${repRID}.getRef.err
    ulimit -a >>${repRID}.getRef.err
    export https_proxy=\${http_proxy}

    # run set the reference name
    if [ "${species_getRef}" == "Mus musculus" ]
    then
      references=\$(echo ${referenceBase}/GRCm${refMoVersion})
    elif [ '${species_getRef}' == "Homo sapiens" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
    else
      exit 1
    fi
    if [ "${spike_getRef}" == "yes" ]
    then
      references=\$(echo \${reference}-S/)
    elif [ "${spike_getRef}" == "no" ]
    then
      reference=\$(echo \${references}/)
    fi

    # retreive appropriate reference appropriate location
    if [ ${referenceBase} == "s3://bicf-references" ]
    then
      aws s3 cp "\${references}" /hisat2 ./ --recursive >>${repRID}.getRef.err
      aws s3 cp "\${references}" /bed ./ --recursive >>${repRID}.getRef.err
    elif [ ${referenceBase} == "/project/BICF/BICF_Core/shared/gudmap/references" ]
    then
      cp -R "\${references}"/hisat2 ./ >>${repRID}.getRef.err
      cp -R "\${references}"/bed ./ >>${repRID}.getRef.err
    fi
    """
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "*.trimData.*"

  input:
    val endsManual_trimData
    path (fastq) from fastqs_trimData

  output:
    path ("*.fq.gz") into fastqs_trimmed
    path ("${repRID}.trimData.log")
    path ("${repRID}.trimData.err")

  script:
    """
    hostname >>${repRID}.trimData.err
    ulimit -a >>${repRID}.trimData.err

    # trim fastqs
    if [ "${endsManual_trimData}" == "se" ]
    then
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err
    elif [ "${endsManual_trimData}" == "pe" ]
    then
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>>${repRID}.trimData.log 2>>${repRID}.trimData.err
    fi
    """
}

// Replicate reference for multiple process inputs
reference.into {
  reference_alignData
  reference_rseqc
}

/*
 * alignData: aligns the reads to a reference database
*/
process alignData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "*.align.{out,err}"

  input:
    val endsManual_alignData
    val stranded_alignData
    path fastq from fastqs_trimmed
    path reference_alignData

  output:
    tuple val ("${repRID}"), path ("${repRID}.sorted.bam"), path ("${repRID}.sorted.bam.bai") into rawBam
    path ("${repRID}.align.out")
    path ("${repRID}.align.err")

  script:
    """
    hostname >${repRID}.align.err
    ulimit -a >>${repRID}.align.err

    # align reads
    if [ "${endsManual_alignData}" == "se" ]
    then
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome ${stranded_alignData} -U ${fastq[0]} 1>${repRID}.align.out 2>${repRID}.align.err
    elif [ "${endsManual_alignData}" == "pe" ]
    then
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome ${stranded_alignData} --no-mixed --no-discordant -1 ${fastq[0]} -2 ${fastq[1]} 1>${repRID}.align.out 2>${repRID}.align.err
    fi
    
    # convert sam to bam and sort and index
    samtools view -1 -@ `nproc` -F 4 -F 8 -F 256 -o ${repRID}.bam ${repRID}.sam 1>>${repRID}.align.out 2>>${repRID}.align.err;
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.bam ${repRID}.bam 1>>${repRID}.align.out 2>>${repRID}.align.err;
    samtools index -@ `nproc` -b ${repRID}.sorted.bam ${repRID}.sorted.bam.bai 1>>${repRID}.align.out 2>>${repRID}.align.err;
    """
}

// Replicate rawBam for multiple process inputs
rawBam.into {
  rawBam_dedupData
  rawBam_inferMetadata
}

/*
 *dedupData: mark the duplicate reads, specifically focused on PCR or optical duplicates
*/
process dedupData {
  tag "${repRID}"
  publishDir "${outDir}/bam", mode: 'copy', pattern: "*.deduped.bam"
  publishDir "${logsDir}", mode: 'copy', pattern: "*.dedup.{out,err}"

  input:
    set val (repRID), path (inBam), path (inBai) from rawBam_dedupData

  output:
    tuple val ("${repRID}"), path ("${repRID}.sorted.deduped.bam"), path ("${repRID}.sorted.deduped.bam.bai") into dedupBam
    path ("${repRID}.dedup.out")
    path ("${repRID}.dedup.err")

  script:
    """
    hostname >${repRID}.dedup.err
    ulimit -a >>${repRID}.dedup.err

    # remove duplicated reads
    java -jar /picard/build/libs/picard.jar MarkDuplicates I=${inBam} O=${repRID}.deduped.bam M=${repRID}.deduped.Metrics.txt REMOVE_DUPLICATES=true 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.deduped.bam ${repRID}.deduped.bam 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    samtools index -@ `nproc` -b ${repRID}.sorted.deduped.bam ${repRID}.sorted.deduped.bam.bai 1>>${repRID}.dedup.out 2>> ${repRID}.dedup.err
    """
}

// Replicate dedup bam/bai for multiple process inputs
dedupBam.into {
  dedupBam_makeBigWig
}

/*
 *Make BigWig files for output
*/
process makeBigWig {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "*.makeBigWig.err"

  input:
    set val (repRID), path (inBam), path (inBai) from dedupBam_makeBigWig

  output:
    path ("${repRID}.bw")

  script:
    """
    bamCoverage -p `nproc` -b ${inBam} -o ${repRID}.bw
    """
}

/*
 *fastqc: run fastqc on untrimmed fastq's
*/
process fastqc {
  tag "${repRID}"
  publishDir "${outDir}/fastqc", mode: 'copy', pattern: "*_fastqc.zip"
  publishDir "${logsDir}", mode: 'copy', pattern: "*.fastq.err"

  input:
    path (fastq) from fastqs_fastqc

  output:
    path ("*_fastqc.zip") into fastqc

  script:
    """
    hostname >${repRID}.fastqc.err
    ulimit -a >>${repRID}.fastqc.err

    # run fastqc
    fastqc *.fastq.gz -o . >>${repRID}.fastqc.err
    """
}

/*
 *rseqc: run RSeQC to collect stats and infer experimental metadata
*/
process inferMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "*.rseqc.err"

  input:
    path script_inferMeta
    path script_aggregateInference
    path reference_rseqc
    set val (repRID), path (inBam), path (inBai) from rawBam_inferMetadata

  output:
    path "infer.csv" into inferedMetadata
    path "${inBam.baseName}.tin.xls" into tin
    path "${repRID}.insertSize.inner_distance_freq.txt" optional true into innerDistance


  script:
    """
    hostname >${repRID}.rseqc.err
    ulimit -a >>${repRID}.rseqc.err

    # infer experimental setting from dedup bam
    infer_experiment.py -r ./bed/genome.bed -i "${inBam}" >${repRID}.rseqc.log

    endness=`bash inferMeta.sh endness ${repRID}.rseqc.log`
    fail=`bash inferMeta.sh fail ${repRID}.rseqc.log`
    if [ \${endness} == "PairEnd" ] 
    then
      percentF=`bash inferMeta.sh pef ${repRID}.rseqc.log`
      percentR=`bash inferMeta.sh per ${repRID}.rseqc.log`
      inner_distance.py -i "${inBam}" -o ${repRID}.insertSize -r ./bed/genome.bed
    elif [ \${endness} == "SingleEnd" ]
    then
      percentF=`bash inferMeta.sh sef ${repRID}.rseqc.log`
      percentR=`bash inferMeta.sh ser ${repRID}.rseqc.log`
    fi
    if [ \$percentF -gt 0.25 ] && [ \$percentR -lt 0.25 ]
    then
      stranded="forward"
      if [ \$endness == "PairEnd" ]
      then
        strategy="1++,1--,2+-,2-+"
      else
        strategy="++,--"
      fi
    elif [ \$percentR -gt 0.25 ] && [ \$percentF -lt 0.25 ]
    then
      stranded="reverse"
      if [ \$endness == "PairEnd" ]
      then
        strategy="1+-,1-+,2++,2--"
      else
        strategy="+-,-+"
      fi
    else
      stranded="unstranded"
      strategy="us"
    fi

    # calcualte TIN values per feature
    tin.py -i "${inBam}" -r ./bed/genome.bed

    # write infered metadata to file
    echo \${endness},\${stranded},\${strategy},\${percentF},\${percentR},\${fail} > infer.csv
    """
}
