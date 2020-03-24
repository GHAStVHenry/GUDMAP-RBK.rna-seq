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
script_calculateTPM = Channel.fromPath("${baseDir}/scripts/calculateTPM.R")
script_inferMeta = Channel.fromPath("${baseDir}/scripts/inferMeta.sh")

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

// Replicate raw fastqs for multiple process inputs
fastqs.into {
  fastqs_downsampleData
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
    val repRID
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
    
    # Get endedness metadata
    endsMeta=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p endsMeta)
    echo "LOG: endedness metadata parsed: \${endsMeta}" >> ${repRID}.parseMetadata.err
    
    # Manually get endness
    endsManual=\$(python3 ${script_parseMeta} -r ${repRID} -m "${fileMeta}" -p endsManual)
    echo "LOG: endedness manually detected: \${endsManual}" >> ${repRID}.parseMetadata.err

    # Get strandedness metadata
    stranded=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p stranded -e \${endsManual})
    echo "LOG: strandedness metadata parsed: \${stranded}" >> ${repRID}.parseMetadata.err
    
    # Get spike-in metadata
    spike=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentSettingsMeta}" -p spike)
    echo "LOG: spike-in metadata parsed: \${spike}" >> ${repRID}.parseMetadata.err
    
    # Get species metadata
    species=\$(python3 ${script_parseMeta} -r ${repRID} -m "${experimentMeta}" -p species)
    echo "LOG: species metadata parsed: \${species}" >> ${repRID}.parseMetadata.err

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
  endsManual_downsampleData
  endsManual_trimData
  endsManual_alignData
  endsManual_featureCounts
}
stranded.into {
  stranded_alignData
  stranded_featureCounts
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
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.getRef.{out,err}"

  input:
    val spike_getRef
    val species_getRef

  output:
    tuple path ("hisat2", type: 'dir'), path ("bed", type: 'dir'), path ("*.fna"), path ("*.gtf")  into reference
    path ("${repRID}.getRef.{out,err}")
 
  script:
    """
    hostname > ${repRID}.getRef.err
    ulimit -a >> ${repRID}.getRef.err
    export https_proxy=\${http_proxy}

    # run set the reference name
    if [ "${species_getRef}" == "Mus musculus" ]
    then
      references=\$(echo ${referenceBase}/GRCm${refMoVersion})
    elif [ '${species_getRef}' == "Homo sapiens" ]
    then
      references=\$(echo ${referenceBase}/GRCh${refHuVersion})
    else
      echo -e "LOG: ERROR - References could not be set!\nSpecies reference found: ${species_getRef}" >> ${repRID}.getRef.err
      exit 1
    fi
    if [ "${spike_getRef}" == "yes" ]
    then
      references=\$(echo \${reference}-S/)
    elif [ "${spike_getRef}" == "no" ]
    then
      reference=\$(echo \${references}/)
    fi
    echo "LOG: species set to \${references}" >> ${repRID}.getRef.err

    # retreive appropriate reference appropriate location
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
  reference_makeFeatureCounts
  reference_inferMeta
}

/*
 * trimData: trims any adapter or non-host sequences from the data
*/
process trimData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.trimData.{out,err}"

  input:
    val endsManual_trimData
    path (fastq) from fastqs_trimData

  output:
    path ("*.fq.gz") into fastqsTrim
    path ("*_trimming_report.txt") into trimQC
    path ("${repRID}.trimData.{out,err}")

  script:
    """
    hostname > ${repRID}.trimData.err
    ulimit -a >> ${repRID}.trimData.err

    #Trim fastqs using trim_galore
    if [ "${endsManual_trimData}" == "se" ]
    then
      echo "LOG: running trim_galore using single-end settings" >> ${repRID}.trimData.err
      trim_galore --gzip -q 25 --illumina --length 35 --basename ${repRID} -j `nproc` ${fastq[0]} 1>> ${repRID}.trimData.out 2>> ${repRID}.trimData.err
    elif [ "${endsManual_trimData}" == "pe" ]
    then
      echo "LOG: running trim_galore using paired-end settings" >> ${repRID}.trimData.err
      trim_galore --gzip -q 25 --illumina --length 35 --paired --basename ${repRID} -j `nproc` ${fastq[0]} ${fastq[1]} 1>> ${repRID}.trimData.out 2>> ${repRID}.trimData.err
    fi
    """
}

// Replicate trimmed fastqs
fastqsTrim.into {
  fastqsTrim_downsampleData
  fastqsTrim_alignData
}

/*
 * downsampleData: downsample fastq's for metadata inference
 */
process downsampleData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.downsampleData.{out,err}"

  input:
    val endsManual_downsampleData
    path fastq from fastqsTrim_downsampleData

  output:
    path ("sampled.{1,2}.fq") into fastqsSample
    path ("${repRID}.downsampleData.{out,err}")

  script:
    """
    hostname > ${repRID}.downsampleData.err
    ulimit -a >> ${repRID}.downsampleData.err
    export https_proxy=\${http_proxy}

    if [ "${endsManual_downsampleData}" == "se" ]
    then
      echo "LOG: downsampling single-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *trimmed.fq.gz 10000 1> sampled.1.fq 2>> ${repRID}.downsampleData.err
    elif [ "${endsManual_downsampleData}" == "pe" ]
    then
      echo "LOG: downsampling read 1 of paired-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *1.fq.gz 1000000 1> sampled.1.fq 2>> ${repRID}.downsampleData.err
      echo "LOG: downsampling read 2 of paired-end trimmed fastq" >> ${repRID}.downsampleData.err
      seqtk sample -s100 *2.fq.gz 1000000 1> sampled.2.fq 2>> ${repRID}.downsampleData.err
    fi
    """
}

/*
 * alignData: aligns the reads to a reference database
*/
process alignData {
  tag "${repRID}"
  publishDir "${logsDir}", mode: "copy", pattern: "${repRID}.align.{out,err}"

  input:
    val endsManual_alignData
    val stranded_alignData
    path fastq from fastqsTrim_alignData
    path reference_alignData

  output:
    tuple val ("${repRID}"), path ("${repRID}.sorted.bam"), path ("${repRID}.sorted.bam.bai") into rawBam
    path ("*.alignSummary.txt") into alignQC
    path ("${repRID}.align.{out,err}")

  script:
    """
    hostname > ${repRID}.align.err
    ulimit -a >> ${repRID}.align.err

    #Align the reads with Hisat 2
    if [ "${endsManual_alignData}" == "se" ]
    then
      echo "LOG: running Hisat2 with single-end settings" >> ${repRID}.align.err
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome ${stranded_alignData} -U ${fastq[0]} --summary-file ${repRID}.alignSummary.txt --new-summary 1>> ${repRID}.align.out 2>> ${repRID}.align.err
    elif [ "${endsManual_alignData}" == "pe" ]
    then
      echo "LOG: running Hisat2 with paired-end settings" >> ${repRID}.align.err
      hisat2 -p `nproc` --add-chrname --un-gz ${repRID}.unal.gz -S ${repRID}.sam -x hisat2/genome ${stranded_alignData} --no-mixed --no-discordant -1 ${fastq[0]} -2 ${fastq[1]} --summary-file ${repRID}.alignSummary.txt --new-summary 1>> ${repRID}.align.out 2>> ${repRID}.align.err
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
    set val (repRID), path (inBam), path (inBai) from rawBam_dedupData

  output:
    tuple path ("${repRID}.sorted.deduped.bam"), path ("${repRID}.sorted.deduped.bam.bai") into dedupBam
    path ("*.deduped.Metrics.txt") into dedupQC
    path ("${repRID}.dedup.{out,err}")

  script:
    """
    hostname > ${repRID}.dedup.err
    ulimit -a >> ${repRID}.dedup.err

    # remove duplicated reads using Picard's MarkDuplicates
    echo "LOG: running picard MarkDuplicates to remove duplicate reads" >> ${repRID}.dedup.err
    java -jar /picard/build/libs/picard.jar MarkDuplicates I=${inBam} O=${repRID}.deduped.bam M=${repRID}.deduped.Metrics.txt REMOVE_DUPLICATES=true 1>> ${repRID}.dedup.out 2>> ${repRID}.dedup.err

    #Use SamTools to sort the now deduped bam file
    echo "LOG: sorting the deduped bam file" >> ${repRID}.dedup.err
    samtools sort -@ `nproc` -O BAM -o ${repRID}.sorted.deduped.bam ${repRID}.deduped.bam 1>> ${repRID}.dedup.out 2>> ${repRID}.dedup.err

    #Use SamTools to index the now sorted deduped bam file
    echo "LOG: indexing the sorted deduped bam file" >> ${repRID}.dedup.err
    samtools index -@ `nproc` -b ${repRID}.sorted.deduped.bam ${repRID}.sorted.deduped.bam.bai 1>> ${repRID}.dedup.out 2>> ${repRID}.dedup.err
    """
}

// Replicate dedup bam/bai for multiple process inputs
dedupBam.into {
  dedupBam_makeFeatureCounts
  dedupBam_makeBigWig
  dedupBam_inferMeta
}

/*
 *Make BigWig files for output
*/
process makeBigWig {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.makeBigWig.{out,err}"
  publishDir "${outDir}/bigwig", mode: 'copy', pattern: "${repRID}.bw"

  input:
    set path (inBam), path (inBai) from dedupBam_makeBigWig

  output:
    path ("${repRID}.bw")
    path ("${repRID}.makeBigWig.{out,err}")

  script:
    """
    hostname > ${repRID}.makeBigWig.err
    ulimit -a >> ${repRID}.makeBigWig.err

    #Run bamCoverage
    echo "LOG: Running bigWig bamCoverage" >> ${repRID}.makeBigWig.err
    bamCoverage -p `nproc` -b ${inBam} -o ${repRID}.bw 1>> ${repRID}.makeBigWig.out 2>> ${repRID}.makeBigWig.err
    """
}

/*
 *Run featureCounts and get the counts, tpm
*/
process makeFeatureCounts {
  tag "${repRID}"
  publishDir "${outDir}/featureCounts", mode: 'copy', pattern: "${repRID}*.countTable.csv"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.makeFetureCounts.{out,err}"

  input:
    path script_calculateTPM
    tuple path (bam), path (bai) from dedupBam_makeFeatureCounts
    path reference_makeFeatureCounts
    val endsManual_featureCounts

  output:
    path ("*.countTable.csv") into counts
    path ("*.featureCounts.summary") into countsQC
    path ("${repRID}.makeFeatureCounts.{out,err}")

  script:
    """
    hostname > ${repRID}.makeFeatureCounts.err
    ulimit -a >> ${repRID}.makeFeatureCounts.err

    #Determine strandedness and setup strandig for featureCounts
    stranding=0;
    if [ "${stranded_featureCounts}" == "--rna-strandness F" ] || [ "${stranded_featureCounts}" == "--rna-strandness FR" ]
      then
      stranding=1
      echo "LOG: strandedness set to stranded [1]" >> ${repRID}.makeFeatureCounts.err
    else
      stranding=0
      echo "LOG: strandedness set to unstranded [0]" >> ${repRID}.makeFeatureCounts.err
    fi;
    #Run featureCounts
    echo "LOG: running featureCounts on the data" >> ${repRID}.makeFeatureCounts.err
    if [ "${endsManual_featureCounts }" == "se" ]
    then
      featureCounts -R SAM -p -G ./genome.fna -T `nproc` -s \${stranding} -a ./genome.gtf -o ${repRID}.featureCounts -g 'gene_name' --primary --ignoreDup ${repRID}.sorted.deduped.bam 1>> ${repRID}.makeFeatureCounts.out 2>> ${repRID}.makeFeatureCounts.err
    elif [ "${endsManual_featureCounts }" == "pe" ]
    then
      featureCounts -R SAM -p -G ./genmome.fna -T `nproc` -s \${stranding} -a ./genome.gtf -o ${repRID}.featureCounts -g 'gene_name' --primary --ignoreDup -B ${repRID}.sorted.deduped.bam 1>> ${repRID}.makeFeatureCounts.out 2>> ${repRID}.makeFeatureCounts.err
    fi

    #Calculate TMP from the resulting featureCounts table
    echo "LOG: calculating TMP with R" >> ${repRID}.makeFeatureCounts.err
    Rscript calculateTPM.R --count "${repRID}.featureCounts" 1>> ${repRID}.makeFeatureCounts.out 2>> ${repRID}.makeFeatureCounts.err
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
 *inferMetadata: run RSeQC to collect stats and infer experimental metadata
*/

process inferMetadata {
  tag "${repRID}"
  publishDir "${logsDir}", mode: 'copy', pattern: "${repRID}.rseqc.{out,err}"

  input:
    path script_inferMeta
    path reference_inferMeta
    set path (inBam), path (inBai) from dedupBam_inferMeta

  output:
    path "infer.csv" into inferedMetadata
    path "${inBam.baseName}.tin.xls" into tin
    path "${repRID}.insertSize.inner_distance_freq.txt" optional true into innerDistance
    path "${repRID}.rseqc.{out,err}" optional true

  script:
    """
    hostname > ${repRID}.rseqc.err
    ulimit -a >> ${repRID}.rseqc.err

    # infer experimental setting from dedup bam
    echo "LOG: running inference from bed file" >> ${repRID}.rseqc.err
    infer_experiment.py -r ./bed/genome.bed -i "${inBam}" > ${repRID}.rseqc.log 2>> ${repRID}.rseqc.err

    echo "LOG: determining endedness and strandedness from file" >> ${repRID}.rseqc.err
    endness=`bash inferMeta.sh endness ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
    fail=`bash inferMeta.sh fail ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
    if [ \${endness} == "PairEnd" ] 
    then
      percentF=`bash inferMeta.sh pef ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
      percentR=`bash inferMeta.sh per ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
      inner_distance.py -i "${inBam}" -o ${repRID}.insertSize -r ./bed/genome.bed 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
    elif [ \${endness} == "SingleEnd" ]
    then
      percentF=`bash inferMeta.sh sef ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
      percentR=`bash inferMeta.sh ser ${repRID}.rseqc.log` 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err
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
    echo -e "LOG: strategy set to \${strategy}\nStranded set to ${stranded}" >> ${repRID}.rseqc.err

    # calcualte TIN values per feature
    tin.py -i "${inBam}" -r ./bed/genome.bed 1>> ${repRID}.rseqc.out 2>> ${repRID}.rseqc.err

    # write infered metadata to file
    echo \${endness},\${stranded},\${strategy},\${percentF},\${percentR},\${fail} > infer.csv 2>> ${repRID}.rseqc.err
    """
}
