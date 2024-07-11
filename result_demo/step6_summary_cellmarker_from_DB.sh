#将分析路径下的所有文献进行cell-marker的提取
head -n 1 /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//paper_process_result/28576768/cellmarker_result.csv > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//cellmarker_result_all.csv
cat /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//paper_process_result/*/cellmarker_result.csv | grep -v "PMID"  >> /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test//cellmarker_result_all.csv
