|master|develop|
|:-:|:-:|
|[![pipeline status](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/badges/master/pipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/commits/master)|[![pipeline status](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/badges/develop/pipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/commits/develop)|
|[![pipeline](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/masterPipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/master)|[![pipeline](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/developPipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/develop)|
|[![nextflow](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/masterNextflow.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/master)|[![nextflow](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/developNextflow.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/develop)|


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4429316.svg)](https://doi.org/10.5281/zenodo.4429316)


RNA-Seq Analytic Pipeline for GUDMAP/RBK
========================================

Introduction
------------
This pipeline was created to be a standard mRNA-sequencing analysis pipeline which integrates with the GUDMAP and RBK consortium data-hub. It is designed to run on the HPC cluster ([BioHPC](https://portal.biohpc.swmed.edu)) at UT Southwestern Medical Center (in conjunction with the standard nextflow profile: config `biohpc.config`)

Authentication:
----------------
The consortium server used must be authentificated with the [deriva authentication client](https://github.com/informatics-isi-edu/gudmap-rbk/wiki/), and remain authentificated till the end of the pipeline run. Prematurely closing the client will result in invalidation of the tokens, and may result in the pipeline failure. The use of long-lived "globus" tokens is on the roadmap for use in the future.

To Run:
-------
* Available parameters:
  * `--deriva` active **credential.json** file from [deriva-auth](https://github.com/informatics-isi-edu/gudmap-rbk/wiki/Uploading-files-via-Deriva-client-tools#from-a-remote-server)
  * `--bdbag` active **cookies.txt** file from [deriva-auth](https://github.com/informatics-isi-edu/gudmap-rbk/wiki/Uploading-files-via-Deriva-client-tools#from-a-remote-server)
  * `--repRID` mRNA-seq replicate RID
  * `--source` consortium server source
    * **dev** = [dev.gudmap.org](dev.gudmap.org) (default, does not contain all data)
    * **staging** = [staging.gudmap.org](staging.gudmap.org) (does not contain all data)
    * **production** = [www.gudmap.org](www.gudmap.org) (***does contain  all data***)
  * `--refMoVersion` mouse reference version ***(optional, default = 38.p6.vM25)***
  * `--refHuVersion` human reference version ***(optional, default = 38.p13.v36)***
  * `--refERCCVersion` human reference version ***(optional, default = 92)***
  * `--upload` option to not upload output back to the data-hub ***(optional, default = false)***
    * **true** = upload outputs to the data-hub
    * **false** = do *NOT* upload outputs to the data-hub
  * `-profile` config profile to use ***(optional)***:
    * defaut = processes on BioHPC cluster
    * **biohpc** = process on BioHPC cluster
    * **biohpc_max** = process on high power BioHPC cluster nodes (=> 128GB nodes), for resource testing
    * **aws_ondemand** = AWS Batch on-demand instant requests
    * **aws_spot** = AWS Batch spot instance requests
  * `--email` email address(es) to send failure notification (comma separated) ***(optional)***:
    * e.g: `--email 'Venkat.Malladi@utsouthwestern.edu,Gervaise.Henry@UTSouthwestern.edu'`

* NOTES:
  * once deriva-auth is run and authenticated, the two files above are saved in ```~/.deriva/``` (see official documents from [deriva](https://github.com/informatics-isi-edu/deriva-client#installer-packages-for-windows-and-macosx) on the lifetime of the credentials)
  * reference version consists of Genome Reference Consortium version, patch release and GENCODE annotation release # (leaving the params blank will use the default version tied to the pipeline version)
    * *current mouse* **38.p6.vM25** = GRCm38.p6 with GENCODE annotation release M25
    * *current human* **38.p13.v36** = GRCh38.p13 with GENCODE annotation release 36
* ***Optional*** input overrides
  * `--refSource` source for pulling references
    * **biohpc** = source references from BICF_Core gudmap reference local location (workflow must be run on BioHPC system)
    * **datahub** = source references from GUDMAP/RBK reference_table location (currently uses dev.gudmap.org)
  * `--inputBagForce` utilizes a local replicate inputBag instead of downloading from the data-hub (still requires accurate repRID input)
    * eg: `--inputBagForce test_data/bag/Q-Y5F6_inputBag_xxxxxxxx.zip` (must be the expected bag structure, this example will not work because it is a test bag)
  * `--fastqsForce` utilizes local fastq's instead of downloading from the data-hub (still requires accurate repRID input)
    * eg: `--fastqsForce 'test_data/fastq/small/Q-Y5F6_1M.R{1,2}.fastq.gz'` (note the quotes around fastq's which must me named in the correct standard [*\*.R1.fastq.gz and/or \*.R2.fastq.gz*] and in the correct order, also consider using `endsForce` if the endness doesn't match submitted value)
  * `--speciesForce` forces the species to be "Mus musculus" or "Homo sapiens", it bypasses a metadata mismatch or an ambiguous species error
    * eg: `--speciesForce 'Mus musculus'`
  * `--endsForce` forces the endness to be "se", or "pe", it bypasses a metadata mismatch error
    * eg: `--endsForce 'pe'`
  * `--strandedForce` forces the strandedness to be "forward", "reverse" or "unstranded", it bypasses a metadata mismatch error
    * eg: `--strandedForce 'unstranded'`
  * `--spikeForce` forces the spike-in to be "false", or "true", it bypasses a metadata mismatch error
    * eg: `--spikeForce 'true'`
* Tracking parameters ([Tracking Site](http://bicf.pipeline.tracker.s3-website-us-east-1.amazonaws.com/)):
  * `--ci` boolean (default = false)
  * `--dev` boolean (default = true)

FULL EXAMPLE:
-------------
```
nextflow run workflow/rna-seq.nf --repRID Q-Y5JA --source production --deriva ./data/credential.json --bdbag ./data/cookies.txt --dev false --upload true -profile biohpc
```
<hr>
Cloud Compatibility:
--------------------
This pipeline is also capable of being run on AWS and DNAnexus. To do so:
* The Nextflow binary needs to contain a custom scm config to allow nextflow to pull the pipeline from the UTSW self-hosted GitLab server (git.biohpc.swmed.edu)
  ```
  providers {
    bicf {
        server = 'https://git.biohpc.swmed.edu'
        platform = 'gitlab'
    }
  }
  ```
  This is required for the use of `nextflow run` or `nextflow pull` pointed directly to the git repo, but also the use in AWS or DNAnexus environments as those both use `nextflow run` directly to that repo. To get around this requirement, there is a clone of the repo hosted on [GitHub](https://github.com/utsw-bicf/gudmap_rbk.rna-seq) which can be used... but the currency of that clone cannot be guarnteed!
### [AWS](https://aws.amazon.com/)
* Build a AWS batch queue and environment either manually or with a template, such as: [Genomics Workflows on AWS](https://docs.opendata.aws/genomics-workflows/)
* The user must have awscli configured with an appropriate authentication (with `aws configure` and access keys) in the environment which nextflow
* Follow the instructions from [AWS](https://docs.aws.amazon.com/cli/latest/reference/batch/submit-job.html) about launching runs, using AWS cli. A template *json* file has been included ([awsExample.json](docs/awsExample.json))
  * `[version]` should be replaced with the pipeline version required (eg: `v2.0.0`)
  * `[credential.json]` should be replaced with the location of the credential file outpted by authentification with Deriva
  * `[cookies.txt]` should be replaced with the location of the cookies file outpted by authentification with Deriva for BDBag
  * `[repRID]` should be replaced with the replicate RID to be analized (eg: `Q-Y5F6`)
  * `[outDir]` should be replaced with the location to save local outputs of the pipeline

  example `aws batch submit-job` command (replaceing the parameters in `[]` with the appropriate values)
  ```
  aws batch submit-job\
    --job-name [Job Name]\
    --job-queue [Queue]\
    --job-definition [Job Definition]\
    --container-overrides command=`cat ./docs/nxf_aws-ci-test.json`
  ```
### [DNAnexus](https://dnanexus.com/) (utilizes the [DNAnexus extension package for Nextflow (XPACK-DNANEXUS)](https://github.com/seqeralabs/xpack-dnanexus))
* Follow the istructions from [XPACK-DNANEXUS](https://github.com/seqeralabs/xpack-dnanexus) about installing and authenticating (a valid license must be available for the extension package from Seqera Labs, as well as a subsription with DNAnexus)
* Follow the instructions from [XPACK-DNANEXUS](https://github.com/seqeralabs/xpack-dnanexus) about launching runs. A template *json* file has been included ([dnanexusExample.json](docs/dnanexusExample.json))
  * `[version]` should be replaced with the pipeline version required (eg: `v2.0.0`)
  * `[credential.json]` should be replaced with the location of the credential file outpted by authentification with Deriva
  * `[cookies.txt]` should be replaced with the location of the cookies file outpted by authentification with Deriva for BDBag
  * `[repRID]` should be replaced with the replicate RID to be analized (eg: `Q-Y5F6`)
  * `[outDir]` should be replaced with the location to save local outputs of the pipeline

  example `dx-run` command
  ```
  dx run nf-dxapp-bicf \
    --delay-workspace-destruction \
    --instance-type mem1_ssd1_v2_x16 \
    --input-json "$(envsubst < ./docs/nxf_dnanexus-ci-test.json)"
  ```
### NOTE:
* File locations used in cloud deployments (auth files and output folder) need to be accessible in that environment (eg s3 location, or DNAnexus location). Local paths cannot be read local locations.
<hr>
To generate you own references or new references:
------------------------------------------
Download the [reference creation script](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/snippets/31).
This script will auto create human and mouse references from GENCODE. It can also create ERCC92 spike-in references as well as concatenate them to GENCODE references automatically. In addition, it can create references from manually downloaded FASTA and GTF files.
<hr>
Errors:
-------
Error reported back to the data-hub are (they aren't thrown on the command line by the pipeline, but rather are submitted (if `--upload true`) to the data-hub for that replicate in the execution run submission):

|Error|Descripton|
|:-|:-:|
|**Too many fastqs detected (>2)**|Data-hub standards and that of this pipeline is for one read-1 fastq and if paired-end, one read\-2 fastq. As a result, the maximum number of fastq's per replicate cannot be more than 2.|
|**No valid fastqs detected (may not match {_.}R{12}.fastq.gz convention)**|Files may be missing, or named in a way that doesn't match the data-hub convention.|
|**Number of fastqs detected does not match submitted endness**|Single-end sequenced replicates can only have one fastq, while paried\-end can only have two (see above).|
|**Number of reads do not match for R1 and R2**|For paired\-end sequenced studies the number of reads in read\-1 fastq must match that of read\-2. This error is usually indicative of uploading of currupted, trunkated, or wrong fastq files.|
|**There is an error with the structure of the fastq**|The fastq's fail a test of their structure. This error is usually indicative of uploading of currupted, trunkated, or wrong fastq files.|
|**Infered species does not match for R1 and R2**|The species inferred from each read does not match. This error is usually indicative of uploading of wrong fastq files.|
|**Infered species confidence is low**|The confidence of the species inferrence call is low. This is usually indicative of very low quality samples.|
|**Infered sequencing type is not mRNA-seq**|The sequence type inferred is not mRNA-seq. This is usually indicative of uploading wrong fastq files.|
|**Infered sequencing type does not match for R1 and R2**|The sequencing type inferred from each read does not match. This error is usually indicative of uploading of wrong fastq files.|
|**Infered species confidence is low**|The confidence of the species inferrence call is low AND 3 sets of a random sampling of the fastq's do not match. This is usually indicative of very low quality samples.|
|**Submitted metadata does not match inferred**|All required metadata for analysis of the data is internally inferred by the pipeline, if any of those do not match the submitted metadata, this error is detected to notify of a potential error. The mismatched metadata will be listed.|

<hr>
[**CHANGELOG**](CHANGELOG.md)
<hr>

Credits
=======
This workflow is was developed by [Bioinformatic Core Facility (BICF), Department of Bioinformatics](http://www.utsouthwestern.edu/labs/bioinformatics/)

PI
--
Venkat S. Malladi\
*Faculty Associate & Director*\
Bioinformatics Core Facility\
UT Southwestern Medical Center\
<a href="https://orcid.org/0000-0002-0144-0564" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">orcid.org/0000-0002-0144-0564</a>\
[venkat.malladi@utsouthwestern.edu](mailto:venkat.malladi@utsouthwestern.edu)


Developers
----------
Gervaise H. Henry\
*Computational Biologist*\
Department of Urology\
UT Southwestern Medical Center\
<a href="https://orcid.org/0000-0001-7772-9578" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">orcid.org/0000-0001-7772-9578</a>\
[gervaise.henry@utsouthwestern.edu](mailto:gervaise.henry@utsouthwestern.edu)

Jonathan Gesell\
*Computational Biologist*\
Bioinformatics Core Facility\
UT Southwestern Medical Center\
<a href="https://orcid.org/0000-0001-5902-3299" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">orcid.org/0000-0001-5902-3299</a>\
[johnathan.gesell@utsouthwestern.edu](mailto:jonathn.gesell@utsouthwestern.edu)

Jeremy A. Mathews\
*Computational Biologist*\
Bioinformatics Core Facility\
UT Southwestern Medical Center\
<a href="https://orcid.org/0000-0002-2931-1430" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">orcid.org/0000-0002-2931-1430</a>\
[jeremy.mathews@utsouthwestern.edu](mailto:jeremy.mathews@utsouthwestern.edu)

Please cite in publications: Pipeline was developed by BICF from funding provided by **Cancer Prevention and Research Institute of Texas (RP150596)**.

<hr>
<hr>

Pipeline Directed Acyclic Graph
-------------------------------

![dag](docs/dag.png "DAG")
