gc()
library(optparse)

option_list=list(
  make_option("--endness",action="store",type='character',help="Infered Endness"),
  make_option("--stranded",action="store",type='character',help="Infered Strandedness"),
  make_option("--strategy",action="store",type='character',help="Infered Sequencing Strategy"),
  make_option("--percentF",action="store",type='double',help="Percent Forward Aligned Reads"),
  make_option("--percentR",action="store",type='double',help="Percent Reverse Aligned Reads"),
  make_option("--percentFail",action="store",type='double',help="Percent Failed Aligned Reads"),
  make_option("--tin",action="store",type='character',help="TIN File")
)
opt=parse_args(OptionParser(option_list=option_list))
rm(option_list)

if (length(setdiff(c("endness","stranded","strategy","percentF","percentR","percentFail","tin","help"),names(opt))) != 0){
  stop(paste0("No input missing ",setdiff(c("endness","stranded","strategy","percentF","percentR","percentFail","tin","help"),names(opt)),", exiting."))
} else if (!file.exists(opt$tin)){
  stop("No tin file passed, exiting.")
}

if (opt$endness == "PairEnd"){
  endness <- "pe"
} else if (opt$endness == "SingleEnd"){
  endness <- "se"
}

percentF <- round(opt$percentF,digits=2)
percentR <- round(opt$percentR,digits=2)
percentFail <- round(opt$percentFail,digits=2)

repRID <- basename(gsub(".sorted.deduped.tin.xls","",opt$tin))

tin <- read.delim(opt$tin,sep="\t")

tin.min <- round(min(tin$TIN),digits=2)
tin.med <- round(median(tin$TIN),digits=2)
tin.max <- round(max(tin$TIN),digits=2)
tin.sd <- round(sd(tin$TIN),digits=2)

output <- paste(endness,opt$stranded,opt$strategy,percentF,percentR,percentFail,tin.min,tin.med,tin.max,tin.sd,sep=",")

write(output,file="infer.csv")
