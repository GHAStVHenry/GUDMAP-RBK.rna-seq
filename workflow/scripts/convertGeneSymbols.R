gc()
library(optparse)

option_list=list(
  make_option("--repRID",action="store",type='character',help="Replicate RID")
)
opt=parse_args(OptionParser(option_list=option_list))
rm(option_list)

countTable <- read.csv(paste0(opt$repRID,".countData.countTable.csv"), stringsAsFactors=FALSE)
geneID <- read.delim("geneID.tsv", header=FALSE, stringsAsFactors=FALSE)
Entrez <- read.delim("Entrez.tsv", header=FALSE, stringsAsFactors=FALSE)
convert <- merge(x=Entrex,y=geneID,by.x="V1",by.y="V1",all.x=TRUE)



convert <- data.frame(geneID=countTable$Geneid)
convert <- merge(x=convert,y=geneID[,1:2],by.x="geneID",by.y="V2",all.x=TRUE)
convert <- merge(x=convert,y=Entrez,by.x="V1",by.y="V1",all.x=TRUE)
convert[is.na(convert$V2),3] <- ""
convert <- convert[,-1]
colnames(convert) <- c("GeneID","EntrezID")
convert <- unique(convert)

output <- merge(x=convert,y=countTable,by.x="GeneID",by.y="Geneid",all.x=TRUE)

write.table(output,file=paste0(opt$repRID,".tpmTable.csv"),sep=",",row.names=FALSE,quote=FALSE)
