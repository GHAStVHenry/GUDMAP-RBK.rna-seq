#!/usr/bin/env R
#convertGeneSymbols.R
#*
#* --------------------------------------------------------------------------
#* Licensed under MIT (https://git.biohpc.swmed.edu/gudmap_rbk/rna-seq/-/blob/14a1c222e53f59391d96a2a2e1fd4995474c0d15/LICENSE)
#* --------------------------------------------------------------------------
#*

gc()
library(optparse)

option_list=list(
  make_option("--repRID",action="store",type='character',help="Replicate RID")
)
opt=parse_args(OptionParser(option_list=option_list))
rm(option_list)

countTable <- read.csv(paste0(opt$repRID,".countTable.csv"), stringsAsFactors=FALSE)
geneID <- read.delim("geneID.tsv", header=FALSE, stringsAsFactors=FALSE)
Entrez <- read.delim("Entrez.tsv", header=FALSE, stringsAsFactors=FALSE)

convert <- data.frame(gene_name=countTable$gene_name)
convert <- merge(x=convert,y=geneID[,1:2],by.x="gene_name",by.y="V2",all.x=TRUE)
convert <- merge(x=convert,y=Entrez,by.x="V1",by.y="V1",all.x=TRUE)
convert[is.na(convert$V2),3] <- ""
convert <- convert[,-1]
colnames(convert) <- c("GeneID","EntrezID")
convert <- unique(convert)

output <- merge(x=convert,y=countTable[,c("gene_name","gene_id","count","tpm")],by.x="GeneID",by.y="gene_name",all.x=TRUE)
colnames(output) <- c("GENCODE_Gene_Symbol","NCBI_GeneID","Ensembl_GeneID","count","tpm")
output <- output[,c(1,3,2,4:5)]

write.table(output,file=paste0(opt$repRID,"_tpmTable.csv"),sep=",",row.names=FALSE,quote=FALSE)
