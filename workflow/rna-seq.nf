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
        file("*") into dataPaths

    script:
        """
        hostname
        ulimit -a
        unzip ${bdbag}
        python3 ${baseDir}/scripts/modifyFetch.py -f \$(echo "${bdbag}" | cut -d'.' -f1)
        bdbag --materialize "\$(echo "${bdbag}" | cut -d'.' -f1)"
        """
 }