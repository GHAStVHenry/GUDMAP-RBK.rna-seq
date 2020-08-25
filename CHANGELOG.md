# v0.0.3 (in development)
**User Facing**
* TPM table:
  * Add Ensembl Gene ID
  * Rename columns: *GENCODE_Gene_Symbol*, *Ensembl_GeneID*, *NCBI_GeneID*
* MultiQC output custom talbes (html+JSON):
  * Run table: *Session ID* and *Pipeline Version*
  * Reference Table: *Species*, *Genome Reference Consortium Build*, *Genome Reference Consortium Patch*, *GENCODE Annotation Release* (ouputs both human and mouse versions)

**Background**
* Add GeneSymbol/EnsemblID/EntrezID translation files to references

*Known Bugs*
* outputBag does not contain fetch for processed data
* Does not include automatic data upload

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
