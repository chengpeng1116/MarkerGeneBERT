#根据关键词从PUBMED得到PMID列表，保存于update_pmid.csv 
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/batch_download_PMID.R  --keyword '"Animals"[MeSH Terms] AND ("Single-Cell Analysis"[MeSH Terms] OR ("single-cell" AND "expression"))' --mindate 2020  --maxdate 2020  --output_file_path "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/update_pmid_cellmarker2.csv" &&\


#将update_pmid.csv文件中的PMID与DB_PMID_status.csv进行比较，输出本数据库未收录的PMID列表，保存于文件update_PMID_checked.csv
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/update_DB.py  --update_pmid "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/update_pmid_cellmarker2.csv" --DB_PMID_status /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//DB_PMID_status.csv --output_path /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//update_PMID_checked.csv
