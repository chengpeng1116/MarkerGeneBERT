#Download XML formatted literature from tidyPMC using R script. If the download is successful, it will be updated to the database directory. If the download fails, the PMID will be saved to update-PMID_checked_PMC_download_fail.csv
Rscript ~./script/extract_paper_info.R ~/result_demo/update_PMID_checked.csv ~./result_demo/DB_source ~./result_demo/update_PMID_checked_PMC_download_fail.csv &&\


#Using scihub, download the literature in update-PMID_checked_PMC_download_fail.csv in PDF format using scihub
~./bin/scihub -ns -ow N -s ~./result_demo/update_PMID_checked_PMC_download_fail.csv  -O ~./result_demo/sub_step3_pdf_download/temp_pdf/  -u https://sci-hub.ren


#Transfer PDF format literature to database directory
Rscript ~./script//pipeline_sh/mv_pdf_download_fromTemp.R ~./result_demo/sub_step3_pdf_download/temp_pdf/ ~./result_demo/DB_source
