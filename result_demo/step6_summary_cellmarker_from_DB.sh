#Parse the PDF of the literature under the analysis path into txt format # Extract all the literature under the analysis path using cell markers
head -n 1 ~./result_demo//paper_process_result/26940531/cellmarker_result.csv > ~./result_demo/cellmarker_result_all.csv
cat ~./result_demo/paper_process_result/*/cellmarker_result.csv | grep -v "PMID"  >> ~./result_demo/cellmarker_result_all.csv
