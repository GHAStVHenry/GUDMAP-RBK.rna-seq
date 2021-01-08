|master|develop|
|:-:|:-:|
|[![pipeline status](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/badges/master/pipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/commits/master)|[![pipeline status](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/badges/develop/pipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/commits/develop)|
|[![pipeline](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/masterPipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/master)|[![pipeline](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/developPipeline.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/develop)|
|[![nextflow](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/masterNextflow.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/master)|[![nextflow](https://gudmap_rbk.pages.biohpc.swmed.edu/rna-seq/badges/developNextflow.svg)](https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/tree/develop)|

<!--
[![DOI]()]()
-->
RNA-Seq Analytic Pipeline for GUDMAP/RBK
========================================

Introduction
------------
This pipeline was created to be a standard mRNA-sequencing analysis pipeline which integrates with the GUDMAP and RBK consortium data-hub. It is designed to run on the HPC cluster ([BioHPC](https://portal.biohpc.swmed.edu)) at UT Southwestern Medical Center (in conjunction with the standard nextflow profile: config `biohpc.config`)

Cloud Compatibility:
--------------------
This pipeline is also capable of being run on AWS. To do so:
* Build a AWS batch queue and environment either manually or with [aws-cloudformantion](https://console.aws.amazon.com/cloudformation/home?#/stacks/new?stackName=Nextflow&templateURL=https://s3.amazonaws.com/aws-genomics-workflows/templates/nextflow/nextflow-aio.template.yaml)
* Edit one of the aws configs in workflow/config/
  * Replace workDir with the S3 bucket generated
  * Change region if different
  * Change queue to the aws batch queue generated
* The user must have awscli configured with an appropriate authentication (with `aws configure` and access keys) in the environment which nextflow will be run
* Add `-profile` with the name aws config which was customized

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
    * eg: `--fastqsForce 'test_data/fastq/small/Q-Y5F6_1M.R{1,2}.fastq.gz'` (note the quotes around fastq's which must me named in the correct standard [*\*.R1.fastq.gz and/or \*.R2.fastq.gz*] and in the correct order)
  * `--speciesForce` forces the species to be "Mus musculus" or "Homo sapiens", it bypasses ambiguous species error
    * eg: `--speciesForce 'Mus musculus'`
* Tracking parameters ([Tracking Site](http://bicf.pipeline.tracker.s3-website-us-east-1.amazonaws.com/)):
  * `--ci` boolean (default = false)
  * `--dev` boolean (default = true)

FULL EXAMPLE:
-------------
```
nextflow run workflow/rna-seq.nf --repRID Q-Y5JA --source production --deriva ./data/credential.json --bdbag ./data/cookies.txt --dev false --upload true -profile biohpc
```

To run a set of replicates from study RID:
------------------------------------------
Run in repo root dir:
* `sh workflow/scripts/splitStudy.sh [studyRID]`
It will run in parallel in batches of 5 replicatesRID with 30 second delays between launches.\
NOTE: Nextflow "local" processes for all replicates will run on the node/machine the bash script is launched from... consider running the study script on the BioHPC's SLURM cluster (use `sbatch`).

Errors:
-------
Error reported back to the data-hub are (they aren't thrown on the command line by the pipeline, but rather are submitted (if `--upload true`) to the data-hub for that replicate in the execution run submission):

|Error|Descripton|
|:-|:-:|
|**Too many fastqs detected (>2)**|Data-hub standards and that of this pipeline is for one read-1 fastq and if paired-end, one read\-2 fastq. As a result, the maximum number of fastq's per replicate cannot be more than 2.|
|**Number of fastqs detected does not match submitted endness**|Single-end sequenced replicates can only have one fastq, while paried\-end can only have two (see above).|
|**Number of reads do not match for R1 and R2**|For paired\-end sequenced studies the number of reads in read\-1 fastq must match that of read\-2. This error is usually indicative of uploading of currupted, trunkated, or wrong fastq files|
|**Inference of species returns an ambiguous result**|Species of the replicate is done by aligning a random subset of 1 million reads from the data to both the human and mouse reference genomes. If there isn't a clear difference between the alignment rates (`>=40%` of one species, but `<40%` of the other), then this error is detected.|
|**Submitted metadata does not match inferred**|All required metadata for analysis of the data is internally inferred by the pipeline, if any of those do not match the submitted metadata, this error is detected to notify of a potential error.|

<hr>
[**CHANGELOG**](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/blob/develop/CHANGELOG.md)
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
