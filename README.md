# MarkerGeneBERT
a natural language processing (NLP) system designed to extract critical information from the literature regarding species, tissues, cell types, and cell marker genes in the context of single-cell sequencing studies


## Requisites
python 3.9<br>
R 4.4< br >
sci-hub : https://github.com/suqingdong/scihub
grobid : https://github.com/kermitt2/grobid


## Downloading the Model
en_core_sci_scibert :https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_scibert-0.5.4.tar.gz
en_ner_jnlpba_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_jnlpba_md-0.5.4.tar.gz
en_ner_craft_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_craft_md-0.5.4.tar.gz
en_ner_bionlp13cg_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz
en_ner_bc5cdr_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bc5cdr_md-0.5.4.tar.gz

## Instalation
pip install scispacy
pip install <Model local URL>


##  run MarkerGeneBERT
python CellMarker_sh.py --cfg CellMarker.cfg 


