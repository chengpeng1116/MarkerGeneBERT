# MarkerGeneBERT
a natural language processing (NLP) system designed to extract critical information from the literature regarding species, tissues, cell types, and cell marker genes in the context of single-cell sequencing studies


## Requisites
python 3.9<br>
R 4.4<br>
sci-hub : https://github.com/suqingdong/scihub<br>
grobid : https://github.com/kermitt2/grobid<br>


## Downloading the Model
en_core_sci_scibert :https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_scibert-0.5.4.tar.gz<br>
en_ner_jnlpba_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_jnlpba_md-0.5.4.tar.gz<br>
en_ner_craft_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_craft_md-0.5.4.tar.gz<br>
en_ner_bionlp13cg_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bionlp13cg_md-0.5.4.tar.gz<br>
en_ner_bc5cdr_md : https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_ner_bc5cdr_md-0.5.4.tar.gz<br>

## Instalation
pip install scispacy<br>
pip install <Model local URL><br>


##  run MarkerGeneBERT
python CellMarker_sh.py --cfg CellMarker.cfg<br> 


