# v1.0.3
**User Facing**
* Endness metadata "Single Read" changed to "Single End" in data-hub, pipeline updated to handle

**Background**
* Add memory limit (75%) per thread for samtools sort (#108)
* Remove parsing restrictions for submitted stranded/spike/species (#105, #106)
* Pass unidentified ends instead of overwriting it as unknown
* Move fastqc process before trim to catch fastq errors (#107)
* Only use fastq's that match *[_.]R[1-2].fastq.gz naming convention (#107)
* Add error output for no fastq's
* Update input bag export config to only fetch fastq's that match *[_.]R[1-2].fastq.gz naming convention
* Remove check for multiple fastq check in parse metadata (redundant and no longer valid)
* Handle blank submitted endness better
* Don't use file.csv from inputBag to parse manual endness, use counted from getData
* Detect malformed fastq's (#107)
* Restrict sampled alignment process to use >32GB nodes on BioHPC (#108)
* Use nproc**-1** for alignment processes (#108)

*Known Bugs*
* Override params (inputBag, fastq, species) aren't checked for integrity
* Authentication files and tokens must be active (active auth client) for the duration of the pipeline run (until long-lived token utilization included)

<hr>

# v1.0.2
**User Facing**

**Background**
* Fix spelling in config file for process of failed fastq to upload error message (#104)

*Known Bugs*
* Override params (inputBag, fastq, species) aren't checked for integrity
* Authentication files and tokens must be active (active auth client) for the duration of the pipeline run (until long-lived token utilization included)

<hr>

# v1.0.1
**User Facing**

**Background**
* Split non-metadata mismatch error handling proces into 2, 1 to handle fastq errors and one for species errors (BUG FIX #101)
* Add known errors to integration CI tests (ambiguous species, trunkated fastq, R1/R2 mismatch (#103)
* Fix pre exeuction run fails uploading of execution run RID to tracking site (#96, #97)
* Change CI replicate count badge CI to count all execution runs that match major version

*Known Bugs*
* Override params (inputBag, fastq, species) aren't checked for integrity
* Authentication files and tokens must be active (active auth client) for the duration of the pipeline run (until long-lived token utilization included)

<hr>

# v1.0.0
**User Facing**
* Add link to reference builder script
* Output median TIN to mRNA_QC table

**Background**
* Change consistency test to check if +/- 5% of standard
* Change tool version checker for badges to use latest tag
* Utilize pipeline tracking and qc AWS tables

*Known Bugs*
* Override params (inputBag, fastq, species) aren't checked for integrity
* Authentication files and tokens must be active (active auth client) for the duration of the pipeline run (until long-lived token utilization included)

<hr>

# v0.1.0
**User Facing**
* Add option to pull references from datahub
* Add option to send email on workflow error, with pipeline error message
* Add versions and paper references of software used to report
* Upload input bag
* Upload execution run
* Upload mRNA QC
* Create and upload output bag
* Add optional to not upload
* Update references to use bags
* Update to newer references (GRCh38.p13.v36 and GRCm38.p6.vM25)
* Use production server for data-hub reference call
* Error pipeline if submitted does not match infered
* Update execution run with "Success" or "Error"
* Error if fastq error (>2, if pe != 2, if se !=1)
* Error if pe and line count of R1 != R2
* Error if ambiguous species inference
* Remove non fastq from inputBag from the export bag config level

**Background**
* Remove (comment out) option to pull references from S3
* Make pull references from BioHPC default (including in biohpc.config)
* Start using new gudmaprbk dockerhub (images autobuilt)
* Moved consistency checks to be fully python
* Changed order of steps so that fastqc is done after the trim step
* Change docker images to production
* Add automated version badges
* Only calculate/report tin values on regular chromosomes (from gtf)
* Change inputBag fetch to manifest then validate (if fail fetch missing and revalidate up to 3 times)
* Retry getData and trimData processes up to once
* Make inputBag export config to create inputBag with only small txt file for CI unit test of getData (and update test)

*Known Bugs*
* Override params (inputBag, fastq, species) aren't checked for integrity

<hr>

# v0.0.3
**User Facing**
* TPM table:
  * Add Ensembl Gene ID
  * Rename columns: *GENCODE_Gene_Symbol*, *Ensembl_GeneID*, *NCBI_GeneID*
* MultiQC output custom tables (html+JSON):
  * Run table: *Session ID* and *Pipeline Version*
  * Reference Table: *Species*, *Genome Reference Consortium Build*, *Genome Reference Consortium Patch*, *GENCODE Annotation Release* (outputs both human and mouse versions)
* Add inputBag override param (`inputBagForce`) [`*.zip`]
  * Uses provided inputBag instead of downloading from data-hub
  * Still requires matching repRID input param
* Add fastq override param (`fastqsForce`) [`R1`,`R2`]
  * Uses provided fastq instead of downloading from data-hub
  * Still requires matching repRID input param and will pull inputBag from data-hub to access submitted metadata for reporting
* Add species override param (`speciesForce`) [`Mus musculus` or `Homo sapiens`]
  * forces the use of the provided species
  * ignores inferred ambiguous species

**Background**
* Add GeneSymbol/EnsemblID/EntrezID translation files to references

*Known Bugs*
* outputBag does not contain fetch for processed data
* Does not include automatic data upload
* Override params (inputBag, fastq, species) aren't checked for integrity

<hr>

# v0.0.2
**User Facing**
* Output:
  * inputBag
  * outputBag
* Remove gene details from tpm table
* Add EntrezID translation to tpm table (from version specific reference)

**Background**
* Add GeneSymbol/EnsemblID/EntrezID translation files to references

*Known Bugs*
* outputBag does not contain fetch for processed data
* Does not include automatic data upload

<hr>

# v0.0.1
**INITIAL BETA VERSION**\
Does not include automatic data upload\
This version is for initial upload of test data to GUDMAP/RBK data-hub for internal integration
<hr>
