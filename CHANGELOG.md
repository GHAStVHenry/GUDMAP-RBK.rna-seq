# v1.0.0 (in development)
**User Facing**
* Add option to pull references from datahub
* Add option to send email on workflow error, with pipeline error message
* Add versions and paper references of software used to report
* Upload input bag
* Upload execution run
* Upload mRNA QC
* Create and upload output bag
* Add optional to not upload

**Background**
* Remove (comment out) option to pull references from S3
* Make pull references from BioHPC default (including in biohpc.config)
* Start using new gudmaprbk dockerhub (images autobuilt)
* Moved consistency checks to be fully python
* Changed order of steps so that fastqc is done after the trim step

*Known Bugs*
* Datahub reference pull uses dev.gudmap.org as source until referencencs are placed on production
* outputBag does not contain fetch for processed data
* Does not include automatic data upload
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
