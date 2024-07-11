#使用R脚本从tidyPMC下载XML格式的文献。下载成功则更新至数据库目录下，下载失败则将PMID保存至update_PMID_checked_PMC_download_fail.csv
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/extract_paper_info.R /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//update_PMID_checked.csv /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/DB_source /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//update_PMID_checked_PMC_download_fail.csv &&\


#使用scihub，将update_PMID_checked_PMC_download_fail.csv中的文献利用scihub下载PDF格式的文献
##下载量小时用 ： /mnt/icfs/personal/xiaolingzhang/.local/bin/scihub -ns -ow N -s /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//update_PMID_checked_PMC_download_fail.csv  -O /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/temp_pdf/  -u https://sci-hub.ren
#cd  /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_0.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_1.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_2.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_3.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_4.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_5.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_6.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_7.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_8.sh  
# qsub /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/run_pdf_9.sh  

#将PDF格式的文献转存至数据库目录下
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/mv_pdf_download_fromTemp.R /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//sub_step3_pdf_download/temp_pdf/ /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/DB_source
