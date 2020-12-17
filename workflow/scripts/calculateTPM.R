gc()
library(optparse)

option_list=list(
  make_option("--count",action="store",type='character',help="Count File")
)
opt=parse_args(OptionParser(option_list=option_list))
rm(option_list)

if (!("count" %in% names(opt))){
  stop("No count file passed, exiting.")
} else if (!file.exists(opt$count)) {
  stop("Count file doesn't exist, exiting.")
}

repRID <- basename(gsub(".countData","",opt$count))

count <- read.delim(opt$count, comment.char="#") # if featureCounts file changes structure, be sure to update count and Length columns below
colnames(count)[1] <- "gene_name"
colnames(count)[7] <- "gene_id"
colnames(count)[8] <- "count"
count <- count[,c(1,7,2:6,8)]

rpk <- count$count/count$Length/1000

scale <- sum(rpk)/1000000

tpm <- rpk/scale

output <- cbind(count,tpm)
colnames(output)[7] <- "count"

write.table(output,file=paste0(repRID,"_countTable.csv"),sep=",",row.names=FALSE,quote=FALSE)
