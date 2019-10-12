#!/usr/bin/env nextflow

// Define input variables
params.bdbag = "${baseDir}/../test_data/Study_Q-Y4H0.zip"

params.outDir = "${baseDir}/../output"

// Parse input variables
bdbag = Channel
  .fromPath(params.bdbag)
  .ifEmpty { exit 1, "bdbag zip file not found: ${params.bdbag}" }

outDir = params.outDir


/*
 * getData: fetch study files from consortium with downloaded bdbag.zip
 * python must be loaded prior to nextflow run, because conda env create from .yml doesn't work with nextflow loaded module (either process in-line, or config file)
 */
 process getData {
     publishDir "${outDir}/temp/getData", mode: "symlink"
     conda "${baseDir}/conf/conda.env.bdbag.yml"

     input:
        file bdbag

    output:
        file("**/*.R*.fastq.gz") into fastqPaths
        file("**/File.csv") into filePaths
        file("**/Experiment Settings.csv") into experimentSettingsPaths
        file("**/Experiment.csv") into experimentPaths

    script:
        """
        hostname
        ulimit -a
        study=\$(echo "${bdbag}" | cut -d'.' -f1)
        echo LOG: \${study}
        unzip ${bdbag}
        echo LOG: bdgag unzipped
        python3 ${baseDir}/scripts/modifyFetch.py --fetchFile \${study}
        echo LOG: fetch file filtered for only .fastq.gz
        #bdbag --materialize "\$(echo "${bdbag}" | cut -d'.' -f1)"
        sh ${baseDir}/scripts/bdbagFetch.sh \${study}
        echo LOG: bdbag fetched
        sh ${baseDir}/scripts/renameFastq.sh \${study}
        echo LOG: fastq.gz files renamed to replicate RID
        """
 }