## /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/R

args <- commandArgs(T)

temp_pdf_path  <- args[1]
db_path <-args[2]


i=0
all_file = length(list.files(temp_pdf_path))
for(each_file in list.files(temp_pdf_path)){
  
  pmid = gsub(".pdf","",each_file)
  file.copy(  paste0(temp_pdf_path, "/",each_file) ,  paste0(db_path,"/",pmid,"/",each_file))
  i=i+1
  print(paste0("done ", i ," from total  " , all_file))
}