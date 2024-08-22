# MarkerGeneBERT
a natural language processing (NLP) system designed to extract critical information from the literature regarding species, tissues, cell types, and cell marker genes in the context of single-cell sequencing studies

## docker image
docker pull chengpeng1116/markergenebert:V1<br> 
docker run  run --rm -i -t  --ulimit core=0   chengpeng1116/markergenebert:V1<br> 

We provide Docker images (markergenebert.image) and Dockerfile for building, and we recommend using Docker containers that have already been set up, as some models that build MarkerGeneBert have a large amount of data and are affected by download speed<br> 
When building a container using Dockerfile, we suggest downloading the wget models in advance and building the Docker container by docker run -v your_local_download_path:/opt/spacy_model 
## Conda environment setup
### Requisites
python(v3.9.0)<br>
R(v4.3.2)<br>
sci-hub : https://github.com/suqingdong/scihub<br>
grobid : https://github.com/kermitt2/grobid<br>

### Download the NER_Model from scispacy
en_core_sci_scibert :https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_scibert-0.5.4.tar.gz<br>
en_ner_jnlpba_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_jnlpba_md-0.5.4.tar.gz<br>
en_ner_craft_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_craft_md-0.5.4.tar.gz<br>
en_ner_bionlp13cg_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz<br>
en_ner_bc5cdr_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bc5cdr_md-0.5.4.tar.gz<br>

### Instalation NER Model
pip install pandas<br>
pip install spacy<br>
pip install scispacy<br>
pip install NER_Model_local_URL<br>



##  run MarkerGeneBERT
cd ~./result_demo/<br>
python /opt/docker_file/MarkerGeneBERT/script/CellMarker_sh.py --cfg CellMarker.cfg<br> 
- CellMarker_sh.py will return step1-6.sh , the user needs to run the scripts for step1-step6 in sequence
- The Marker-related sentence classification model of MarkerGeneBERT is stored in https://github.com/chengpeng1116/MarkerGeneBERT/Database/MarkerGeneBERT/mode-best and embedded in the step5_cellmarker_extract.sh script
```
(1)  bash step1_check_PMID_status.sh    ## Used to detect the analysis status of various literature under the result_demo, and save the analysis status results to DB-PMID_status.csv 
(2)  bash step2_update_DB.sh    ## Obtain the PMID list from PUBMED based on keywords and save it in update_pmid.csv
(3)  bash step3_download_PMID.sh    ## Download XML formatted literature from tidyPMC using R script. If the download is successful, it will be updated to the database directory. If the download fails, the PMID will be saved to update_PMID_checked_PMC_download_fail.csv
(4)  bash step4_pdf2txt_a04.sh    ## Parse the PDF of the literature under the result_demo into txt format
(5)  bash step5_cellmarker_extract.sh    ## Extract cell markers from all literature under the result_demo
(6)  bash step6_summary_cellmarker_from_DB.sh    ## Summarize all results under the result_demo
```





