library(RISmed)
library("optparse")


option_list <- list(
  make_option(c("-k", "--keyword"), type = "character"),
  make_option(c("-n", "--mindate"), type = "integer"),
  make_option(c("-x", "--maxdate"), type = "integer"),
  make_option(c("-o", "--output_file_path"), type = "character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



# key='"Animals"[MeSH Terms] AND "Single-Cell Analysis"'
data=EUtilsSummary(opt$keyword,db="pubmed",retmax=122451,mindate=opt$mindate,maxdate=opt$maxdate)

pid=data@PMID

output_data <-data.frame("PMID"=pid)
#  summary(pid)

write.csv(output_data,file=opt$output_file_path,row.names = F ,quote = F )