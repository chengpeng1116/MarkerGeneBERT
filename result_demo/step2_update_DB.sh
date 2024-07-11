#Obtain the PMID list from PUBMED based on keywords and save it in update_pmid.csv
Rscript ~./script/batch_download_PMID.R  --keyword '"Animals"[MeSH Terms] AND ("Single-Cell Analysis"[MeSH Terms] OR ("single-cell" AND "expression"))' --mindate 2017  --maxdate 2022  --output_file_path "~./result_demo/demo_update_pmid_cellmarker2.csv" &&\


#Compare the PMID in the update_pmid.csv file with DB-PMID_status.csv, output a list of PMIDs not included in this database, and save it in the file update-PMID_check.csv
python ~./script/update_DB.py  --update_pmid "~./result_demo/demo_update_pmid_cellmarker2.csv" --DB_PMID_status ~./result_demo/DB_PMID_status.csv --output_path ~./result_demo/update_PMID_checked.csv
