## /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/R

library(easyPubMed)
library(tidypmc)
library(tidyverse)
library(europepmc)
library(rentrez)
args <- commandArgs(T)

#paper_id_path <-"d:/download_collection/manual_pubmedid.csv"
#output_path <-"d:/download_collection/"

paper_id_path <- args[1]
output_path <-args[2]
fail_path <-args[3]

paper_id <-read.csv(paper_id_path,header=T)

if(nrow(paper_id)>0){
  
  #########  extract all paper abstract info
  for(each_paper_id in paper_id$PMID){
  
    file_path =paste0(output_path,"/",each_paper_id,"/")
    output_file= paste0(file_path,"abstract.tsv")
  
    if(file.exists(output_file)){
      next
    }
    if(!file.exists(file_path)){
      dir.create(file_path)
    }
  
    ############ 提取摘要和 标题
    loop_count=0
    while(loop_count<5){
  
      tryCatch({
  
        my_pubmed_ids <- get_pubmed_ids(each_paper_id)
        my_data <- fetch_pubmed_data(my_pubmed_ids, encoding = "ASCII")
        df <- table_articles_byAuth(my_data,
                                    included_authors = "first",
                                    max_chars = 2000,
                                    encoding = "ASCII",getKeywords=T)
  
        output_data =paste0("Paper_title: " ,df$title   ,"\nPaper_abstract: " ,df$abstract,"\nPaper_keyword: " ,df$keywords,"\nPaper_date: ", df$year,"\nPaper_journal: ",df$journal)
        write.table(output_data,file=output_file,row.names = F,col.names = F,quote = F,fileEncoding = "UTF-8")
        loop_count <<- 6
  
      }, error = function(e){
        loop_count <<-loop_count+1
      })
  
    }
  
  }

  ##########  get PMC full info 
  
  success_paper <-c()
  for(each_paper_id in paper_id$PMID){
    
    file_path =paste0(output_path,"/",each_paper_id,"/")
    output_file= paste0(file_path,"PMC.tsv")
    if(file.exists(output_file)){
      success_paper <-c(success_paper, each_paper_id)
      next
    }
    print(paste0(each_paper_id," : get pmc start"))
    output_datafrme<- data.frame()
    
    boolFalse<-T
    while(boolFalse){
      tryCatch({
        results <- entrez_link(dbfrom = 'pubmed', id = each_paper_id, db = 'pmc')
        boolFalse <-F
      },error=function(e){
        boolFalse <-T
      })
    }
    ############ 提取 PMC 开源的文献 ，生成PMC.tsv
    temp  <- results$links['pubmed_pmc']
    if(is.null(names(temp))){
      next
    } 
    
    if(is.na(names(temp))){
      next
    }else{
      
      loop_count<-0
      while(loop_count<4){
        
        tryCatch({ 

          doc <- pmc_xml(paste0("PMC",temp))
          txt <- pmc_text(doc)
          fig = pmc_caption(doc)
          
          all_ref_segment = str_extract_all(as.character(doc) ,paste0(" [a-zA-Z]+?<xref.*?>\\d+?</xref>"))
          
          for(each_all_ref_segment in unlist(all_ref_segment)){
            
            tryCatch({
              keyword <-gsub("<xref","",str_extract(each_all_ref_segment, "\\S+?<xref"))
              
              first_ref <-gsub("</xref|>","",str_extract(each_all_ref_segment , '>\\d+?</xref'))
              
              grep_index <- grep(paste0(keyword,first_ref,"\\D*?"),txt$text)
              if(length(grep_index)>0){
                for(each_grep_index in grep_index){
                  txt$text[each_grep_index] <-gsub(paste0(keyword,first_ref,"\\D+?"),keyword    ,txt$text[each_grep_index])
                }
              }
            },error=function(e){
              print(1)
            })
          }
          keyword_line_list <-list()
          all_line= seq(1,nrow(txt))
          for(each_keyword in c("title","abstract","back|intro","result|case","discussion","method")){
          
            each_keyword_row =sort( unique(c(grep(each_keyword,substr(txt$section,1,30),ignore.case = T),grep(toupper(each_keyword),substr(txt$section,1,100),ignore.case = F))))
            if(length(each_keyword_row)>0){
              keyword_line_list[[each_keyword]] <- intersect(each_keyword_row , all_line)
              all_line <-setdiff(all_line,each_keyword_row )
            }
          }
          
          
          keyword = c("title","abstract","method","discussion")
          for(each_keyword in keyword){
            
            if(!is.null(keyword_line_list[[each_keyword]])){
              output_datafrme <-rbind(output_datafrme , data.frame("col1"= paste0("Paper_",each_keyword,": " ,paste0(txt$text[keyword_line_list[[each_keyword]]],collapse = ", ") ))) 
              output_datafrme <-rbind(output_datafrme , data.frame("col1"= txt$text[keyword_line_list[[each_keyword]]]))
            }
          }
          each_keyword= "back|intro"
          if(!is.null(keyword_line_list[[each_keyword]])){
            
            output_datafrme <-rbind(output_datafrme , data.frame("col1"= paste0("Paper_intro: @@" ,paste0(txt$text[keyword_line_list[[each_keyword]]],collapse = " @@ ") ))) 
            output_datafrme <-rbind(output_datafrme , data.frame("col1"= txt$text[keyword_line_list[[each_keyword]]]))
          }
          each_keyword=  "result|case"
          if(!is.null(keyword_line_list[[each_keyword]])){
            subset_txt = txt[keyword_line_list[[each_keyword]],]
            output_datafrme <-rbind(output_datafrme , data.frame("col1"= paste0("Paper_section1: " ,paste0(subset(subset_txt[keyword_line_list[[each_keyword]],],paragraph==1 )$text,collapse = ", ") ))) 
            if( 2 %in% unique(subset_txt$paragraph)){
              output_datafrme <-rbind(output_datafrme , data.frame("col1"= paste0("Paper_section2: " ,paste0(subset(subset_txt[keyword_line_list[[each_keyword]],],paragraph==2 )$text,collapse = ", ") ))) 
            }
            output_datafrme <-rbind(output_datafrme , data.frame("col1"= txt$text[keyword_line_list[[each_keyword]]]))
          }
            if(nrow(fig)>0){
              output_datafrme <-rbind(output_datafrme , data.frame("col1"= fig$text))
            }
  
          
          
          write.table(output_datafrme,file=output_file,col.names = F,row.names = F,quote = F,fileEncoding = "UTF-8")
          success_paper <-c(success_paper, each_paper_id)
          loop_count  <<- 6
        }, error = function(e){
          loop_count <<-loop_count+1
        })
        
      }
    }
  }
  error_dataframe <-data.frame("PMID"=setdiff(paper_id$PMID,success_paper ) )
  write.table(error_dataframe,file=fail_path,col.names = F,row.names = F,quote = F,sep=",",fileEncoding = "UTF-8")
  
  
  
}

















