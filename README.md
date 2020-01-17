|*master*|*develop*|
|:-:|:-:|
|[![Build Status](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/badges/master/build.svg)](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/commits/master)|[![Build Status](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/badges/develop/build.svg)](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/commits/develop)|
<!--
[![DOI]()]()
-->
GUDMAP/RBK RNA-Seq Pipeline
===========================

Introduction
------------
This pipeline was created to be a standard mRNA-sequencing analysis pipeline which integrates with the GUDMAP and RBK consortium data-hub.

To Run:
-------
* Available parameters:
  * *--deriva* active ```credential.json``` file from [deriva-auth](https://github.com/informatics-isi-edu/gudmap-rbk/wiki/Uploading-files-via-Deriva-client-tools#from-a-remote-server)
  * *--bdbag* active ```cookies.txt``` file from [deriva-auth](https://github.com/informatics-isi-edu/gudmap-rbk/wiki/Uploading-files-via-Deriva-client-tools#from-a-remote-server)
  * *--repRID* mRNA-seq replicate RID\
  note: once deriva-auth is run and authenticated, the two files above are saved in ```~/.deriva/``` (see official documents from [deriva](https://github.com/informatics-isi-edu/deriva-client#installer-packages-for-windows-and-macosx) on the lifetime of the credentials)
* FULL EXAMPLE:
  ```
  nextflow run workflow/rna-seq.nf --deriva ./data/credential.json --bdbag ./data/cookies.txt --repRID Q-Y5JA
  ```



[**CHANGELOG**](https://git.biohpc.swmed.edu/BICF/gudmap_rbk/rna-seq/blob/develop/CHANGELOG.md)

---

Credits
=======
This worklow is was developed by [Bioinformatic Core Facility (BICF), Department of Bioinformatics](http://www.utsouthwestern.edu/labs/bioinformatics/)

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
*Computational Intern*\
Bioinformatics Core Facility\
UT Southwestern Medical Center\
<a href="https://orcid.org/0000-0002-2931-1430" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">orcid.org/0000-0002-2931-1430</a>\
[jeremy.mathews@utsouthwestern.edu](mailto:jeremy.mathews@utsouthwestern.edu)

Please cite in publications: Pipeline was developed by BICF from funding provided by **Cancer Prevention and Research Institute of Texas (RP150596)**.
