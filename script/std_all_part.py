import argparse
from dataclasses import replace
from operator import index
import sys
parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--work_path', type = str, default = None)
parser.add_argument('--paper_id_path', type = str, default = None)
parser.add_argument('--paper_id_start', type = str, default = None)
parser.add_argument('--paper_id_end', type = str, default = None)

root_sh_dir = sys.path[0] 
# root_sh_dir = "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh"
args = parser.parse_args()
### /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python

global nlp_model_jnlpba_extra_tissue,nlp_model_jnlpba_extra_gene_cell,nlp_model_scibert,nlp_model_scibert_abbr,nlp_model_jnlpba,nlp_model_craft,nlp_model_bionlp13cg,nlp_model_bc5cdr
global  punctuation_word,gene_collect_path,cell_collect_path,tissue_collect_path,nlp_model_scibert_NER


work_path=args.work_path
paper_id_path=args.paper_id_path
paper_id_start=args.paper_id_start
paper_id_end=args.paper_id_end
#work_path = "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cellmarker_pdf_part/"
#paper_id_path = "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/test_folder/paper_id.csv"

#do_part= ["pdf2txt","txt2segment","extract_potentiall_valuable_segment","add_segment_text_clean","add_predict_result","cell_marker_extract_main"]
do_part= ["txt2segment","extract_potentiall_valuable_segment","add_segment_text_clean","add_predict_result","cell_marker_extract_main"]
##### common import 
gene_collect_path = root_sh_dir + "/../Database/Extra_named_entity_database/gene_collect.csv"
cell_collect_path = root_sh_dir +"/../Database/Extra_named_entity_database/cell_ontology_cellname.txt"
tissue_collect_path =root_sh_dir + "/../Database/Extra_named_entity_database/cell_ontology_tissue.txt"
tissue_cancer_collect_path = root_sh_dir +"/../Database/Extra_named_entity_database/tissue_cancer.txt"
nlp_model_path=root_sh_dir+"/../Database/MarkerGeneBert/model-best/model-best" 
from tkinter.tix import COLUMN
import pandas as pd 
import os 
import os 
from tqdm import tqdm
import spacy
import codecs
import scispacy
import pandas as pd 
import re
import copy
from copy import deepcopy
from collections import Counter
import pickle
from spacy.lang.en.stop_words import STOP_WORDS
import string
from spacy.tokens import DocBin
from scispacy.abbreviation import AbbreviationDetector
import warnings
warnings.filterwarnings("ignore")
import multiprocessing
import csv
from spacy.matcher import Matcher
import inflect

nlp_model_jnlpba_extra_tissue = spacy.load('en_ner_jnlpba_md')
nlp_model_jnlpba_extra_tissue.remove_pipe("ner")
nlp_model_jnlpba_extra_gene_cell = spacy.load('en_ner_jnlpba_md')
nlp_model_jnlpba_extra_gene_cell.remove_pipe("ner")
nlp_model_scibert = spacy.load('en_core_sci_scibert')
nlp_model_scibert_NER = spacy.load('en_core_sci_scibert', disable = ["tagger", "parser", "attribute_ruler", "lemmatizer"])
nlp_model_jnlpba = spacy.load('en_ner_jnlpba_md', disable = ["tagger", "parser", "attribute_ruler", "lemmatizer"])
nlp_model_craft = spacy.load('en_ner_craft_md', disable = ["tagger", "parser", "attribute_ruler", "lemmatizer"])
nlp_model_bionlp13cg = spacy.load('en_ner_bionlp13cg_md', disable = ["tagger", "parser", "attribute_ruler", "lemmatizer"])
nlp_model_bc5cdr = spacy.load("en_ner_bc5cdr_md", disable = ["tagger", "parser", "attribute_ruler", "lemmatizer"])
nlp_model_scibert_abbr = spacy.load('en_core_sci_scibert')
nlp_model_scibert_abbr.add_pipe("abbreviation_detector")


gene_collect = pd.read_csv(gene_collect_path)
gene_patterns = [{"label": "manual_marker", "pattern": x}  for x in gene_collect['gene_collect']]
#### import cell
cell_contact = pd.read_table(cell_collect_path,header=None,encoding='utf-8')
cell_contact.columns=["index_name","cell_collect"]
cell_collect = []
for x in cell_contact['cell_collect']:
    x=str(x)
    cell_collect.append(x)
    cell_collect.append(x.lower())
    cell_collect.append(x.title())
    if x.endswith('s'):
        temp = x[:-1]
        cell_collect.append(temp)
        cell_collect.append(temp.lower())
        cell_collect.append(temp.title())
    else:
        temp = x+"s"
        cell_collect.append(temp)
        cell_collect.append(temp.lower())
        cell_collect.append(temp.title())
#### make gene & cell pattern        
cell_pattern = [{"label": "manual_cell", "pattern": [{"LOWER": x}]} for x in set(cell_collect) if x not in ['cells','cell',"Cell",'Cells']]
ruler = nlp_model_jnlpba_extra_gene_cell.add_pipe("entity_ruler")
ruler.add_patterns(gene_patterns)
ruler.add_patterns(cell_pattern)
#### tissue pattern 
p = inflect.engine()
tissue_contact = pd.read_table(tissue_collect_path,header=None)
delete_keyword = list(set([ x     for x  in  tissue_contact.iloc[:,1]  if  x.endswith("system")  ]))
delete_keyword +=["cell cluster","tissue"]
tissue_contact= tissue_contact.drop(tissue_contact[tissue_contact.iloc[:,1].isin(delete_keyword)].index)
tissue_contact.drop_duplicates( keep='first',inplace=True)
tissue_collect = [ x.lower()  for x in tissue_contact.iloc[:,1] ]
tissue_collect =tissue_collect + [ p.plural(x)  for x in tissue_collect]
tissue_collect =tissue_collect + [ p.singular_noun(x)  for x in tissue_collect]
tissue_collect = tissue_collect+ ["peripheral blood"]
tissue_pattern =  [{"label": "tissue_from_cellontology", "pattern": x}  for x in list(set(tissue_collect))]
tissue_cancer_contact = pd.read_table(tissue_cancer_collect_path,header=None)
tissue_cancer_contact = [ x.lower()  for x in tissue_cancer_contact.iloc[:,0] ]
tissue_cancer_contact =tissue_cancer_contact + [ p.plural(x)  for x in tissue_collect]
tissue_cancer_contact =tissue_cancer_contact + [ p.singular_noun(x)  for x in tissue_collect]
tissue_pattern = tissue_pattern + [{"label": "tissue_cancer", "pattern": x}  for x in list(set(tissue_cancer_contact))]
ruler = nlp_model_jnlpba_extra_tissue.add_pipe("entity_ruler")
ruler.add_patterns(tissue_pattern)


################ function defination 
def cheak_mkdir(inpath):
    if os.path.exists(inpath):
        pass
    else:
        os.makedirs(inpath)


def abbreviation2fullname(each_line,abbreviations_dict):
    each_line_abbreviation_temp=[]
    abbreviation_word = list()
    for doc in nlp_model_craft(each_line) :
        if doc.text  in abbreviations_dict.keys():
            each_line_abbreviation_temp.append(abbreviations_dict[doc.text])
        else:
            each_line_abbreviation_temp.append(doc.text)
    return(" ".join(each_line_abbreviation_temp))


def fullname2abreviation(each_line,abbreviations_dict):
    each_line_abbreviation_temp=[]
    front_word_content=""
    for doc in nlp_model_craft(each_line) :
        if doc.text not in abbreviations_dict.keys():
            if doc.text not in ["(",")"]:
                each_line_abbreviation_temp.append(doc.text)
                front_word_content = doc.text
        else:
            if  re.search(pattern=re.escape(front_word_content),string=abbreviations_dict[doc.text] )!=None:
                each_line_abbreviation_temp=each_line_abbreviation_temp[:-len(abbreviations_dict[doc.text].split(" "))]
                each_line_abbreviation_temp.append(abbreviations_dict[doc.text])
            else:
                each_line_abbreviation_temp.append(abbreviations_dict[doc.text])
    each_line_change_abbreviation = " ".join(each_line_abbreviation_temp)
    return(each_line_change_abbreviation)





def pdf_fig2txt(work_path,paper_id_path):
    temp_order = open(work_path+"/paper_process_result/pdf_fig2txt.sh","w")
    print("bash "+ root_sh_dir+"/GetTextFromPDF_20230228.sh "+ work_path+"/paper_process_result/ &",file=temp_order ) 
    print("python3 "+root_sh_dir+"/extract_txt2cellname.py -w "+ work_path+"/paper_process_result/" + ' -p '+ paper_id_path,file=temp_order ) 
    temp_order.close()


def txt2segment(file_path,output_file_path,abbreviations_dict_path):
    if os.path.exists(file_path):
        if os.path.exists(output_file_path) and os.stat(output_file_path).st_size>1024:
            print(output_file_path + "  has exist")
        else :
            print(file_path + " :start txt2segment ")
            output_file = codecs.open(output_file_path,"w","utf-8")
            tt =  codecs.open(file_path,'r','utf-8')
            mylist = tt.read().splitlines() 
            tt.close()
            doc_collect=list()
            for each_mylist in mylist:
                try:
                    if each_mylist.startswith("Paper_"):
                        print(each_mylist.strip(),file=output_file)
                        doc_1 = nlp_model_scibert_abbr(re.sub("\., ",". ",re.sub("Paper_.*?: ","",each_mylist)).strip())  
                        doc_collect.append(doc_1)
                        for sent in doc_1.sents:
                            print(sent.text,file=output_file)
                    else:
                        doc_1 = nlp_model_scibert_abbr(each_mylist)  
                        doc_collect.append(doc_1)
                        for sent in doc_1.sents:
                            print(sent.text,file=output_file)
                except :
                    pass
            output_file.close()
            abbreviations_dict =dict()
            for each_doc_collect in doc_collect:
                temp_abrv= each_doc_collect._.abbreviations
                for abrv in temp_abrv:
                    if not str(abrv) in abbreviations_dict.keys() and  not re.search("\(|\)|\[|\]\^" , str(abrv._.long_form)) and not re.search("\(|\)|\[|\]\^" , str(abrv)) and str(abrv)!="cell" and str(abrv)!="cells" and len(str(abrv))>1:
                        try: 
                            if not re.search( str(abrv).lower() ,str(abrv._.long_form).split(" ")[0].lower()):
                                abbreviations_dict[str(abrv)]=str(abrv._.long_form)
                                abbreviations_dict[str(abrv)+'s']=str(abrv._.long_form)
                                if str(abrv).endswith("s"):
                                    if str(abrv._.long_form).endswith("s"):
                                        abbreviations_dict[str(abrv)[:-1]]=str(abrv._.long_form)[:-1]
                                    else:
                                        abbreviations_dict[str(abrv)[:-1]]=str(abrv._.long_form)
                        except :
                            pass
            dict_file_output = open(abbreviations_dict_path , 'wb')
            pickle.dump(abbreviations_dict,dict_file_output)
            dict_file_output.close()
            print(file_path + ":done ")
    else:
        print("file do not exist : "+file_path)


def info_extract( Paper_info_dict ,extract_keyword, model,candidata_list=[],back_list=[],output_file=None):
    p = inflect.engine()
    candidata_list = [p.singular_noun(x).lower()  if p.singular_noun(x)  else x.lower()  for x  in  list(set(candidata_list)-set([""]))] 
    candidata_list =list( set(candidata_list)-set(["tissue","organ","pe","tumor","tumour","tissues"]) )
    if extract_keyword=="tissue":
        error_extract_tissue = ["tissue","organ","pe","tumor","tissues","organs","specimens","specimen","tissue specimens","tissue specimen","organoid","tissue explant","cancer","left","organoid","neoplasm","singlecell","cancer tissue","biopsy","column","body","droplet""biopsy","column","body","seed","digit"]
        #### get need cell name 
        label_collect =[]
        entity_collect=[]
        cell_entity_output_paper_tissue= []
        for ent in nlp_model_craft(Paper_info_dict['Paper_tissue']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_tissue']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)         
        for ent in nlp_model_jnlpba(Paper_info_dict['Paper_tissue']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)      
        for ent in nlp_model_jnlpba_extra_gene_cell(Paper_info_dict['Paper_tissue']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)      
        if "CL" in label_collect :
            cell_entity_output_paper_tissue = cell_entity_output_paper_tissue + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CL")     ]))
        if "CELL" in label_collect:
            cell_entity_output_paper_tissue = cell_entity_output_paper_tissue + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL")     ]))
        if "CELL_TYPE" in label_collect or "CELL_LINE" in label_collect :
            if "CELL_TYPE" in label_collect:
                cell_entity_output_paper_tissue = cell_entity_output_paper_tissue + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_TYPE")     ]))
            else:
                cell_entity_output_paper_tissue = cell_entity_output_paper_tissue + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
        if "manual_cell" in label_collect:
            cell_entity_output_paper_tissue = cell_entity_output_paper_tissue +  list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_cell")     ]))
        cell_entity_output_paper_tissue = [x for x in cell_entity_output_paper_tissue if not re.search("\)|\(",x) ]
        cell_entity_output_paper_tissue = [x for x in cell_entity_output_paper_tissue if not x.startswith("cell") ]
        cell_entity_output_paper_tissue = [x for x in cell_entity_output_paper_tissue if not re.search("single cell",x) and not re.search("single-cell",x) ]
        cell_entity_output_paper_tissue = list(set(cell_entity_output_paper_tissue)-set(['cells','cell',"Cell",'Cells']) )
        cell_entity_output_paper_tissue = list(set([ x  if not p.singular_noun(x)  else p.singular_noun(x) for x in cell_entity_output_paper_tissue ]))
        #print(cell_entity_output_paper_tissue)
#
        label_collect =[]
        entity_collect=[]
        cell_entity_output_paper_method= []
        for ent in nlp_model_craft(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)         
        for ent in nlp_model_jnlpba(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)      
        for ent in nlp_model_jnlpba_extra_gene_cell(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)      
        if "CL" in label_collect :
            cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CL")     ]))
        if "CELL" in label_collect:
            cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL")     ]))
        if "CELL_TYPE" in label_collect or "CELL_LINE" in label_collect :
            if "CELL_TYPE" in label_collect:
                cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_TYPE")     ]))
            else:
                cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
        if "manual_cell" in label_collect:
            cell_entity_output_paper_method = cell_entity_output_paper_method +  list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_cell")     ]))
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not re.search("\)|\(",x) ]
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not x.startswith("cell") ]
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not re.search("single cell",x) and not re.search("single-cell",x) ]
        cell_entity_output_paper_method = list(set(cell_entity_output_paper_method)-set(['cells','cell',"Cell",'Cells']) )
        cell_entity_output_paper_method = list(set([ x  if not p.singular_noun(x)  else p.singular_noun(x) for x in cell_entity_output_paper_method ]))
        #print(cell_entity_output_paper_method)
        cell_entity_output = list(set(cell_entity_output_paper_tissue).intersection(cell_entity_output_paper_method))
        if len(cell_entity_output)==0:
            cell_entity_output = cell_entity_output_paper_tissue
        #####  candidate tissue
        candidate_tissue_collect =  []
        for x in candidata_list:
            for ent in nlp_model_bionlp13cg(x).ents:
                    if ent.label_ not in ['AMINO_ACID',"CELL","CELLULAR_COMPONENT","GENE_OR_GENE_PRODUCT","ORGANISM","ORGANISM_SUBDIVISION","ORGANISM_SUBSTANCE","SIMPLE_CHEMICAL"]:
                        if ent.text ==x:
                            candidate_tissue_collect.append(x.lower())
        candidate_tissue_collect = [ x  for x in candidate_tissue_collect if not re.search("tumor|tumour",x,re.IGNORECASE)  and len(x)>2]
        candidate_tissue_collect = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(candidate_tissue_collect))] 
        candidate_tissue_collect=list(set(candidate_tissue_collect)-set(error_extract_tissue))
        ####
        Paper_tissue_from_line = tissue_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_tissue']),model=nlp_model_jnlpba_extra_tissue)
        Paper_tissue_from_line = [ x  for x in Paper_tissue_from_line if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_tissue_from_line = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_tissue_from_line))] 
        #
        Paper_tissue_collect_model = []
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_tissue']).ents:
            if ent.label_  in ["ANATOMICAL_SYSTEM", "DEVELOPING_ANATOMICAL_STRUCTURE","IMMATERIAL_ANATOMICAL_ENTITY","MULTI-TISSUE_STRUCTURE","ORGAN","TISSUE"]:
                Paper_tissue_collect_model.append(ent.text )   
        Paper_tissue_collect_model = [ x  for x in Paper_tissue_collect_model if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_tissue_collect_model = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_tissue_collect_model))] 
        #
        Paper_tissue_from_line = [ x  for x in Paper_tissue_from_line if x not in   Paper_tissue_collect_model ]
        Paper_tissue_collect = [  each_Paper_tissue_from_line        for each_Paper_tissue_from_line in Paper_tissue_from_line if not  re.search(each_Paper_tissue_from_line ," ".join(Paper_tissue_collect_model) ) ] +Paper_tissue_collect_model
        Paper_tissue_collect = set(Paper_tissue_collect)-set(error_extract_tissue)
        ####
        Paper_title_collect = []
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_title']).ents:
            if ent.label_  in ["ANATOMICAL_SYSTEM", "DEVELOPING_ANATOMICAL_STRUCTURE","IMMATERIAL_ANATOMICAL_ENTITY","MULTI-TISSUE_STRUCTURE","ORGAN","TISSUE"]:
                Paper_title_collect.append(ent.text )   
        Paper_title_collect = [ x  for x in Paper_title_collect if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_title_collect = [p.singular_noun(x).lower()  if len(x)>3 and   p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_title_collect))] 
        Paper_title_collect = set(Paper_title_collect)-set(error_extract_tissue)
        ####        
        Paper_section1_from_line = tissue_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_section1']),model=nlp_model_jnlpba_extra_tissue)
        Paper_section1_from_line = [ x  for x in Paper_section1_from_line if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_section1_from_line = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_section1_from_line))] 
        #
        Paper_section1_collect_model = []
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_section1']).ents:
            if ent.label_  in ["ANATOMICAL_SYSTEM", "DEVELOPING_ANATOMICAL_STRUCTURE","IMMATERIAL_ANATOMICAL_ENTITY","MULTI-TISSUE_STRUCTURE","ORGAN","TISSUE"]:
                Paper_section1_collect_model.append(ent.text )   
        Paper_section1_collect_model = [ x  for x in Paper_section1_collect_model if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_section1_collect_model = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_section1_collect_model))] 
        #
        Paper_section1_from_line = [ x  for x in Paper_section1_from_line if x not in   Paper_section1_collect_model ]
        Paper_section1_collect = [  each_Paper_section1_from_line      for each_Paper_section1_from_line in Paper_section1_from_line if not  re.search(each_Paper_section1_from_line ," ".join(Paper_section1_collect_model) ) ] +Paper_section1_collect_model
        Paper_section1_collect = set(Paper_section1_collect)-set(error_extract_tissue)
        ####
        Paper_method_from_line = tissue_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_method']),model=nlp_model_jnlpba_extra_tissue)
        Paper_method_from_line = [ x  for x in Paper_method_from_line if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_method_from_line = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_method_from_line))] 
        #
        Paper_method_collect_model = []
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_method']).ents:
            if ent.label_  in ["ANATOMICAL_SYSTEM", "DEVELOPING_ANATOMICAL_STRUCTURE","IMMATERIAL_ANATOMICAL_ENTITY","MULTI-TISSUE_STRUCTURE","ORGAN","TISSUE"]:
                Paper_method_collect_model.append(ent.text )   
        Paper_method_collect_model = [ x  for x in Paper_method_collect_model if not re.search("tumor|tumour",x,re.IGNORECASE) and len(x)>2]
        Paper_method_collect_model = [p.singular_noun(x).lower()  if len(x)>3 and  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_method_collect_model))]                 
        #
        Paper_method_from_line = [ x  for x in Paper_method_from_line if x not in   Paper_method_collect_model ]
        Paper_method_collect = [  each_Paper_method_from_line        for each_Paper_method_from_line in Paper_method_from_line if not  re.search(each_Paper_method_from_line ," ".join(Paper_method_collect_model) ) ] +Paper_method_collect_model
        Paper_method_collect = set(Paper_method_collect)-set(error_extract_tissue)
        #    
        temp = list(Paper_title_collect) + list(candidate_tissue_collect) +list(Paper_tissue_collect)+list(Paper_section1_collect)+list(Paper_method_collect)
        all_tissue = Counter([x for x  in temp  if x not in cell_entity_output_paper_tissue+cell_entity_output_paper_method])
        #print(all_tissue)
        extra_info_table = list(Paper_tissue_collect)+list(Paper_title_collect)+list(candidate_tissue_collect) +["temp"]
        if Counter(extra_info_table).most_common(1)[0][1]==3:
            auto_select_tissue =   [  x     for x in Counter(extra_info_table) if Counter(extra_info_table)[x] ==3]
        else:
            if len(all_tissue)==0:
                return("not detect")
            else:
                sort_tissue_dataframe =  pd.DataFrame({"tissue_name": [x for x in all_tissue.keys()] ,"tissue_freq":[all_tissue[x] for x in all_tissue] })  
                sort_tissue_dataframe["keyword_freq"] = 0
                must_use_tissue = [x for x in list(Paper_title_collect) +list(candidate_tissue_collect) if x in all_tissue.keys()]
                Paper_tissue_line = Paper_info_dict['Paper_tissue'].split("&&")
                for i_row in range(0,sort_tissue_dataframe.shape[0]):
                    try:
                        temp_count=0
                        grep1 = [i  for i in range(0,len(Paper_tissue_line)) if re.search(sort_tissue_dataframe["tissue_name"][i_row], Paper_tissue_line[i],re.I)]
                        grep1_2  = [i  for i in range(0,len(Paper_tissue_line)) if re.search(p.plural_noun(sort_tissue_dataframe["tissue_name"][i_row]), Paper_tissue_line[i])]
                        grep1 =list(set(grep1+grep1_2))
                        for each_keyword in ["10x","single cell","single-cell","single-nucleus","single nucleus","umap","tsne","suspen","specimen","dissociation","dissociate","scrna","snrnA","dissected","isolat","surgery","model","digest","deposit","were collected","culture","harvest"]:
                            grep2 = [i  for i in range(0,len(Paper_tissue_line)) if re.search(each_keyword, Paper_tissue_line[i],re.I)]
                            if len(set(grep1).intersection(grep2))>0:
                                temp_count = temp_count+1
                                sort_tissue_dataframe.loc[i_row, 'keyword_freq']= temp_count
                    except :
                        temp_count = 0
                        sort_tissue_dataframe.loc[i_row, 'keyword_freq']= temp_count
                ##
                sort_tissue_dataframe["cell_freq"] = 0
                for i_row in range(0,sort_tissue_dataframe.shape[0]):
                    try:
                        temp_count=0
                        grep1 = [i  for i in range(0,len(Paper_tissue_line)) if re.search(sort_tissue_dataframe["tissue_name"][i_row], Paper_tissue_line[i],re.I)]
                        grep1_2  = [i  for i in range(0,len(Paper_tissue_line)) if re.search(p.plural_noun(sort_tissue_dataframe["tissue_name"][i_row]), Paper_tissue_line[i])]
                        grep1 =list(set(grep1+grep1_2))
                        for each_keyword in cell_entity_output:
                            grep2 = [i  for i in range(0,len(Paper_tissue_line)) if re.search(each_keyword, Paper_tissue_line[i],re.I)]
                            if len(set(grep1).intersection(grep2))>0:
                                temp_count = temp_count+1
                                sort_tissue_dataframe.loc[i_row, 'cell_freq']= temp_count
                    except :
                        temp_count = 0
                        sort_tissue_dataframe.loc[i_row, 'cell_freq']= temp_count
                ##
                if len(must_use_tissue)>0:
                    if  max(sort_tissue_dataframe[sort_tissue_dataframe["tissue_name"].isin(must_use_tissue)]["keyword_freq"].tolist())==0:
                        sort_tissue_dataframe.loc[sort_tissue_dataframe["tissue_name"].isin(must_use_tissue), 'keyword_freq']= max(max(sort_tissue_dataframe["keyword_freq"].tolist()),1)
                        sort_tissue_dataframe.loc[sort_tissue_dataframe["tissue_name"].isin(must_use_tissue), 'tissue_freq']=max(max(sort_tissue_dataframe["tissue_freq"].tolist()),1)
                if  max(sort_tissue_dataframe["keyword_freq"].tolist())>0 and  Counter(sort_tissue_dataframe["keyword_freq"].tolist())[max(sort_tissue_dataframe["keyword_freq"].tolist())]!=sort_tissue_dataframe[sort_tissue_dataframe["keyword_freq"]>0].shape[0]:
                    change_index = sort_tissue_dataframe["keyword_freq"]==max(sort_tissue_dataframe["keyword_freq"].tolist())
                    sort_tissue_dataframe.loc[change_index, 'tissue_freq']= max(sort_tissue_dataframe["tissue_freq"].tolist())
                #if  max(sort_tissue_dataframe["tissue_freq"].tolist())>0 and  Counter(sort_tissue_dataframe["tissue_freq"].tolist())[max(sort_tissue_dataframe["tissue_freq"].tolist())]!=sort_tissue_dataframe[sort_tissue_dataframe["tissue_freq"]>0].shape[0]:
                #    change_index = sort_tissue_dataframe["tissue_freq"]==max(sort_tissue_dataframe["tissue_freq"].tolist())
                #    sort_tissue_dataframe.loc[change_index, 'keyword_freq']= max(sort_tissue_dataframe["keyword_freq"].tolist())
                sort_tissue_dataframe["tissue_rank"] = sort_tissue_dataframe["tissue_freq"].rank(method="dense").astype("int")
                sort_tissue_dataframe["keyword_rank"] = sort_tissue_dataframe["keyword_freq"].rank(method="dense").astype("int")
                sort_tissue_dataframe["cell_rank"] = sort_tissue_dataframe["cell_freq"].rank(method="dense").astype("int")
                sort_tissue_dataframe["rank_stat"] = sort_tissue_dataframe["tissue_rank"]+sort_tissue_dataframe["keyword_rank"] +sort_tissue_dataframe["cell_rank"] 
                if max(sort_tissue_dataframe["keyword_freq"].tolist())>0:
                    sort_tissue_dataframe = sort_tissue_dataframe[sort_tissue_dataframe["keyword_freq"]>0]
                else:
                    return("not detect")
                sort_tissue_dataframe.sort_values(by="rank_stat",inplace=True,ascending=False)
                sort_tissue_dataframe["tissue_name_merge"] = sort_tissue_dataframe['tissue_name'] +":" +sort_tissue_dataframe['rank_stat'].astype("str")
                auto_select_tissue = sort_tissue_dataframe.loc[sort_tissue_dataframe["rank_stat"].nlargest(2,keep="all").index , "tissue_name"].tolist()
                sort_tissue_dataframe.to_csv(output_file,index=False,encoding='utf_8_sig')
                #print(sort_tissue_dataframe)
        return("|".join(auto_select_tissue))
    elif extract_keyword=="disease":
        disease_type_output=[]
        temp_line=Paper_info_dict["Paper_title"]    
        try:
            for temp in model(temp_line).ents:
                if temp.label_ in ["DISEASE"]:
                    disease_type_output.append(temp.text.upper()) 
        except Exception as e :
            print("extract disease from model error ,pass extract from model"  )    
        disease_type_output= list(set(disease_type_output))
        if len(disease_type_output)==1:
            return(disease_type_output[0])
        elif   len(disease_type_output)==0:
            return("normal")
        else :
            return_name=[]
            if candidata_list!=[]:
                if len(candidata_list)>0:
                    for each_candidata_list in candidata_list:
                        if any(each_candidata_list in s for s in disease_type_output):    
                            return_name.append(each_candidata_list)  
                    for each_disease_type_output in disease_type_output:
                        if any(each_disease_type_output in s for s in candidata_list):    
                            return_name.append(each_disease_type_output)  
            else:
                if back_list!=[]:
                    for each_back_list in back_list:
                        if any(each_back_list in s for s in disease_type_output):    
                            return_name.append(each_back_list)  
                    for each_disease_type_output in disease_type_output:
                        if any(each_disease_type_output in s for s in back_list):    
                            return_name.append(each_disease_type_output)      
            if len( return_name)>0:
                return(Counter(return_name).most_common(1)[0][0])
            else:
                return(" or ".join(disease_type_output))
    elif extract_keyword=="speices":
        if re.search("patient|individual|human" ,Paper_info_dict["Paper_title"] ,re.IGNORECASE):
            return("human")   
        if re.search("mouse|murine|rat" ,Paper_info_dict["Paper_title"] ,re.IGNORECASE):
            return("mouse")  
        label_collect =[]
        entity_collect=[]
        cell_entity_output_paper_species= []
        for ent in nlp_model_craft(Paper_info_dict['Paper_species']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)
    #
    #
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_species']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)   
    #
    #       
        for ent in nlp_model_jnlpba(Paper_info_dict['Paper_species']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)  
    #
    #     
        for ent in nlp_model_jnlpba_extra_gene_cell(Paper_info_dict['Paper_species']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)   
    #
    #    
        if "CL" in label_collect :
            cell_entity_output_paper_species = cell_entity_output_paper_species + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CL")     ]))
    #
        if "CELL" in label_collect:
            cell_entity_output_paper_species = cell_entity_output_paper_species + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL")     ]))
    #
        if "CELL_TYPE" in label_collect or "CELL_LINE" in label_collect :
            if "CELL_TYPE" in label_collect:
                cell_entity_output_paper_species = cell_entity_output_paper_species + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_TYPE")     ]))
            else:
                cell_entity_output_paper_species = cell_entity_output_paper_species + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
    #
        if "manual_cell" in label_collect:
            cell_entity_output_paper_species = cell_entity_output_paper_species +  list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_cell")     ]))
    #
        cell_entity_output_paper_species = [x for x in cell_entity_output_paper_species if not re.search("\)|\(",x) ]
        cell_entity_output_paper_species = [x for x in cell_entity_output_paper_species if not x.startswith("cell") ]
        cell_entity_output_paper_species = [x for x in cell_entity_output_paper_species if not re.search("single cell",x) and not re.search("single-cell",x) ]
        cell_entity_output_paper_species = list(set(cell_entity_output_paper_species)-set(['cells','cell',"Cell",'Cells']) )
        cell_entity_output_paper_species = list(set([ x  if not p.singular_noun(x)  else p.singular_noun(x) for x in cell_entity_output_paper_species ]))
    #
        label_collect =[]
        entity_collect=[]
        cell_entity_output_paper_method= []
        for ent in nlp_model_craft(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)
    #
        for ent in nlp_model_bionlp13cg(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)         
    #
        for ent in nlp_model_jnlpba(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)     
    # 
        for ent in nlp_model_jnlpba_extra_gene_cell(Paper_info_dict['Paper_method']).ents:
            label_collect.append(ent.label_)
            entity_collect.append(ent.text)   
    #   
        if "CL" in label_collect :
            cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CL")     ]))
    #
        if "CELL" in label_collect:
            cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL")     ]))
    #
        if "CELL_TYPE" in label_collect or "CELL_LINE" in label_collect :
            if "CELL_TYPE" in label_collect:
                cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_TYPE")     ]))
            else:
                cell_entity_output_paper_method = cell_entity_output_paper_method + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
    #
        if "manual_cell" in label_collect:
            cell_entity_output_paper_method = cell_entity_output_paper_method +  list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_cell")     ]))
    #
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not re.search("\)|\(",x) ]
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not x.startswith("cell") ]
        cell_entity_output_paper_method = [x for x in cell_entity_output_paper_method if not re.search("single cell",x) and not re.search("single-cell",x) ]
        cell_entity_output_paper_method = list(set(cell_entity_output_paper_method)-set(['cells','cell',"Cell",'Cells']) )
        cell_entity_output_paper_method = list(set([ x  if not p.singular_noun(x)  else p.singular_noun(x) for x in cell_entity_output_paper_method ]))
        #print(cell_entity_output_paper_method)
        cell_entity_output = list(set(cell_entity_output_paper_species).intersection(cell_entity_output_paper_method))
        if len(cell_entity_output)==0:
            cell_entity_output = cell_entity_output_paper_species
    #
    #
        Paper_title_collect = species_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_title']),model=nlp_model_craft)
        Paper_title_collect = [p.singular_noun(x).lower()  if  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_title_collect))] 
        if Paper_title_collect ==[] and re.search("patient|individual" ,Paper_info_dict["Paper_title"] ,re.IGNORECASE):
            Paper_title_collect=["human"]
    #
        Paper_abstract_collect = species_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_abstract']),model=nlp_model_craft)
        Paper_abstract_collect = [p.singular_noun(x).lower()  if  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_abstract_collect))] 
        if Paper_abstract_collect ==[] and re.search("patient|individual" ,Paper_info_dict["Paper_abstract"] ,re.IGNORECASE):
            Paper_abstract_collect=["human"]
    #
    #
        Paper_species_collect = species_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_species']),model=nlp_model_craft)
        Paper_species_collect = [p.singular_noun(x).lower()  if  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_species_collect))] 
        if Paper_species_collect ==[] and re.search("patient|individual" ,Paper_info_dict["Paper_species"] ,re.IGNORECASE):
            Paper_species_collect=["human"]
    #
        Paper_section1_collect = species_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_section1']),model=nlp_model_craft)
        Paper_section1_collect = [p.singular_noun(x).lower()  if  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_section1_collect))] 
        if Paper_section1_collect ==[] and re.search("patient|individual" ,Paper_info_dict["Paper_section1"] ,re.IGNORECASE):
            Paper_section1_collect=["human"]
    #
        Paper_method_collect = species_extract(each_line=re.sub("-"," ",Paper_info_dict['Paper_method']),model=nlp_model_craft)
        Paper_method_collect = [p.singular_noun(x).lower()  if  p.singular_noun(x)  else x.lower()  for x  in  list(set(Paper_method_collect))] 
        if Paper_method_collect ==[] and re.search("patient|individual" ,Paper_info_dict["Paper_method"] ,re.IGNORECASE):
            Paper_method_collect=["human"]
    #
    #
        temp = list(Paper_title_collect) + list(candidata_list) +list(Paper_species_collect)+list(Paper_section1_collect)+list(Paper_method_collect)+list(Paper_abstract_collect)
        all_species =Counter([x for x  in temp  if x not in cell_entity_output_paper_species+cell_entity_output_paper_method])
        #print(all_species)
        extra_info_table = list(Paper_species_collect)+list(Paper_title_collect)+list(candidata_list) +list(Paper_abstract_collect)+["temp"]
        auto_select_species =   [  x     for x in Counter(extra_info_table) if Counter(extra_info_table)[x] ==Counter(extra_info_table).most_common(1)[0][1]]
        #print(auto_select_species)
        if Counter(extra_info_table).most_common(1)[0][1]>=3 and len(auto_select_species)==1 and len(Paper_abstract_collect)<2:
            auto_select_species =   auto_select_species
            #print(1)
        else:
            if len(all_species)==0:
                return("not detect")
            else:
                sort_species_dataframe =  pd.DataFrame({"species_name": [x for x in all_species.keys()] ,"species_freq":[all_species[x] for x in all_species] })  
                sort_species_dataframe["keyword_freq"] = 0
                must_use_tissue = [x for x in list(Paper_title_collect) +list(candidata_list) if x in all_species.keys()]
                Paper_species_line = Paper_info_dict['Paper_species'].split("&&")
                for i_row in range(0,sort_species_dataframe.shape[0]):
                    try:
                        temp_count=0
                        grep1 = [i  for i in range(0,len(Paper_species_line)) if re.search(sort_species_dataframe["species_name"][i_row], Paper_species_line[i],re.I)]
                        grep1_2  = [i  for i in range(0,len(Paper_species_line)) if re.search(p.plural_noun(sort_species_dataframe["species_name"][i_row]), Paper_species_line[i])]
                        grep1 =list(set(grep1+grep1_2))
                        for each_keyword in ["10x","single cell","single-cell","single-nucleus","single nucleus","umap","tsne","suspen","specimen","dissociation","dissociate","scrna","snrnA","dissected","isolat","surgery","model","digest","deposit","were collected","culture","harvest","genome","aligned","obtain"]:
                            grep2 = [i  for i in range(0,len(Paper_species_line)) if re.search(each_keyword, Paper_species_line[i],re.I)]
                            if len(set(grep1).intersection(grep2))>0:
                                temp_count = temp_count+1
                                sort_species_dataframe.loc[i_row, 'keyword_freq']= temp_count
                    except :
                        temp_count = 0
                        sort_species_dataframe.loc[i_row, 'keyword_freq']= temp_count
                ##
                sort_species_dataframe["cell_freq"] = 0
                for i_row in range(0,sort_species_dataframe.shape[0]):
                    try:
                        temp_count=0
                        #grep1 = [i  for i in range(0,len(Paper_species_line)) if re.search(sort_species_dataframe["species_name"][i_row], Paper_species_line[i],re.I)]
                        #grep1_2  = [i  for i in range(0,len(Paper_species_line)) if re.search(p.plural_noun(sort_species_dataframe["species_name"][i_row]), Paper_species_line[i])]
                        #grep1 =list(set(grep1+grep1_2))
                        #for each_keyword in cell_entity_output:
                        # grep2 = [i  for i in range(0,len(Paper_species_line)) if re.search(each_keyword, Paper_species_line[i],re.I)]
                            #if len(set(grep1).intersection(grep2))>0:
                            # temp_count = temp_count+1
                                #sort_species_dataframe.loc[i_row, 'cell_freq']= temp_count
                        temp_count = Paper_info_dict['Paper_species'].count(sort_species_dataframe["species_name"][i_row])
                        sort_species_dataframe.loc[i_row, 'cell_freq']= temp_count
                    except :
                        temp_count = 0
                        sort_species_dataframe.loc[i_row, 'cell_freq']= temp_count
                ##
                if len(must_use_tissue)>0:
                    if  max(sort_species_dataframe[sort_species_dataframe["species_name"].isin(must_use_tissue)]["keyword_freq"].tolist())!=0:
                        sort_species_dataframe = sort_species_dataframe[sort_species_dataframe["species_name"].isin(must_use_tissue)]
                if  max(sort_species_dataframe["keyword_freq"].tolist())>0 and  Counter(sort_species_dataframe["keyword_freq"].tolist())[max(sort_species_dataframe["keyword_freq"].tolist())]!=sort_species_dataframe[sort_species_dataframe["keyword_freq"]>0].shape[0]:
                    change_index = sort_species_dataframe["keyword_freq"]==max(sort_species_dataframe["keyword_freq"].tolist())
                    sort_species_dataframe.loc[change_index, 'species_freq']= max(sort_species_dataframe["species_freq"].tolist())
                #if  max(sort_tissue_dataframe["tissue_freq"].tolist())>0 and  Counter(sort_tissue_dataframe["tissue_freq"].tolist())[max(sort_tissue_dataframe["tissue_freq"].tolist())]!=sort_tissue_dataframe[sort_tissue_dataframe["tissue_freq"]>0].shape[0]:
                #    change_index = sort_tissue_dataframe["tissue_freq"]==max(sort_tissue_dataframe["tissue_freq"].tolist())
                #    sort_tissue_dataframe.loc[change_index, 'keyword_freq']= max(sort_tissue_dataframe["keyword_freq"].tolist())
                sort_species_dataframe.loc[sort_species_dataframe.species_name.isin(["patient","individual"]),"species_name"]="human"
                sort_species_dataframe = sort_species_dataframe.groupby(['species_name'],as_index=False).agg({'species_freq': 'sum','keyword_freq': 'sum',"cell_freq":"sum"})
                sort_species_dataframe["species_rank"] = sort_species_dataframe["species_freq"].rank(method="dense").astype("int")
                sort_species_dataframe["keyword_rank"] = sort_species_dataframe["keyword_freq"].rank(method="dense").astype("int")
                sort_species_dataframe["cell_rank"] = sort_species_dataframe["cell_freq"].rank(method="dense").astype("int")
                sort_species_dataframe["rank_stat"] = sort_species_dataframe["species_rank"]+sort_species_dataframe["keyword_rank"] +sort_species_dataframe["cell_rank"] 
                if max(sort_species_dataframe["keyword_freq"].tolist())>0:
                    sort_species_dataframe = sort_species_dataframe[sort_species_dataframe["keyword_freq"]>0]
                else:
                    return("not detect")
                sort_species_dataframe.sort_values(by="rank_stat",inplace=True,ascending=False)
                sort_species_dataframe["species_name_merge"] = sort_species_dataframe['species_name'] +":" +sort_species_dataframe['rank_stat'].astype("str")
                sort_species_dataframe.to_csv(output_file,index=False,encoding='utf_8_sig')
                auto_select_species = sort_species_dataframe.loc[sort_species_dataframe["rank_stat"].nlargest(1,keep="all").index , "species_name"].tolist()
                #print(sort_species_dataframe)
                if not len(auto_select_species)==1:
                    temp_sort_species_dataframe = sort_species_dataframe[sort_species_dataframe['species_name'].isin(auto_select_species)]
                    auto_select_species = temp_sort_species_dataframe.loc[temp_sort_species_dataframe["keyword_freq"].nlargest(1,keep="all").index , "species_name"].tolist()
                #   
                #sort_species_dataframe.to_csv(output_file,index=False,encoding='utf_8_sig')
    #
        return("|".join(auto_select_species))



def species_extract(each_line, model=nlp_model_craft):
    p = inflect.engine()
    stat_word_list=[]
    word_entity_dict=dict()
    entity_word= ["TAXON"]
    for each_word in nlp_model_scibert(each_line).ents:
        for temp in model(each_word.text).ents:
            append_word =temp.text.upper()
            if p.singular_noun(append_word):
                append_word=p.singular_noun(append_word)
            stat_word_list.append(append_word)                    
            word_entity_dict.update({append_word:temp.label_})     
    counter_std=0
    outputname=[]
    counter_table = Counter(stat_word_list)  
    for each_word_entity_dict in word_entity_dict.keys():
        if word_entity_dict[each_word_entity_dict] in entity_word:
            if counter_std < counter_table[each_word_entity_dict]:
                outputname=[each_word_entity_dict]
                counter_std = counter_table[each_word_entity_dict]
            elif counter_std == counter_table[each_word_entity_dict]:
                outputname=outputname + [each_word_entity_dict]
    outputname = list(set(outputname)-set(['ANIMAL',"MAMMLIAN"]))
    return(outputname)


def species_extract(each_line,model):
    p = inflect.engine()
    species_type=[]
    for ent in model(each_line.lower()).ents:
        if ent.label_ in ["TAXON"]:
            if not p.singular_noun(ent.text):
                    species_type.append(ent.text)
            else:
                species_type.append(p.singular_noun(ent.text))
    return(list(set(species_type)-set(['ANIMAL',"MAMMLIAN"])))


def tissue_extract(each_line,model):
    p = inflect.engine()
    tissue_type=[]
    for ent in model(each_line.lower()).ents:
        if ent.label_ in ["tissue_from_cellontology","tissue_cancer"]:
            if not p.singular_noun(ent.text):
                    tissue_type.append(ent.text)
            else:
                tissue_type.append(p.singular_noun(ent.text))
    return(list(set(tissue_type)))



def extract_potentiall_valuable_segment(open_file_path,output_file_path,output_fig_segment_file_path,paper_id,abbr_dict_path, append=False):
    p = inflect.engine()
    if  os.path.exists(open_file_path):
        if  os.path.exists(output_file_path) and os.stat(output_file_path).st_size>1024:
            print("dont run : file  exist : "+output_file_path)  
        else:
            print(str(paper_id) +" : start extract_potentiall_valuable_evidence")
            output =  codecs.open(output_file_path,'w','utf-8')
            output_fig =  codecs.open(output_fig_segment_file_path,'w','utf-8')
            print("PMID"+"\t"+"Title"+"\t"+"Date"+"\t"+"Journal"+"\t"+"Species"+"\t"+"Disease"+"\t"+"Tissue_type"+"\t"+"Evidence"+"\t"+"Weight"+"\t"+"Gene_name"+"\t"+"Cell_name"+"\t"+"context_background", file=output )
            #tt =  codecs.open(open_file_path,'r','utf-8')
            #mylist = tt.read().splitlines() 
            #all_row =(" ").join(mylist)
            #fullname_entity_collect = list(set([x.text for x in nlp_model_scibert_NER(all_row).ents]  ))
            #tt.close()
            #species_name = select_disease_species(each_line=all_row, select="species")
            ##### make abbreviations dict{ abbreviations:Full name , abbreviationss:Full name,abbreviation:Full name   }
            if os.path.exists(abbr_dict_path):
                abbreviations_dict = pickle.load(open(abbr_dict_path, 'rb'))
            else:
                print(abbr_dict_path +" not exist")
                abbreviations_dict=dict()
            input_file = codecs.open(open_file_path,"r","utf-8")
            Paper_info_dict=dict(Paper_title="", Paper_abstract="",Paper_intro="",Paper_other_back="",Paper_section1="",Paper_section2="",Paper_method="",Paper_keyword="" ,Paper_species="",Paper_tissue="" ,Paper_date="",Paper_journal=""   ) 
            for each_line in input_file:
                #  if each_line.startswith("Paper_title") :
                #      Paper_info_dict['Paper_title'] =each_line.strip()      
                #  elif each_line.startswith("Paper_abstract") :
                #      Paper_info_dict['Paper_abstract'] =each_line.strip()              
                if each_line.startswith("Paper_intro") :
                    Paper_info_dict['Paper_intro'] =abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line.strip() )), abbreviations_dict=  abbreviations_dict )       
                elif each_line.startswith("Paper_other_back") :
                    Paper_info_dict['Paper_other_back'] =abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line.strip() )), abbreviations_dict=  abbreviations_dict )          
                elif  each_line.startswith("Paper_section1")  :
                    Paper_info_dict['Paper_section1'] =abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line.strip() )), abbreviations_dict=  abbreviations_dict )       
                elif  each_line.startswith("Paper_section2")  :
                    Paper_info_dict['Paper_section2'] =abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line.strip() )), abbreviations_dict=  abbreviations_dict )       
                elif  each_line.startswith("Paper_method")  :
                    Paper_info_dict['Paper_method'] =abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line.strip() )), abbreviations_dict=  abbreviations_dict )          
            input_file.close()   
            input_file = codecs.open(open_file_path,"r","utf-8")
            temp_tissue_line= []
            for all_each_line in input_file:
                for each_line in re.split(r'\.,\s*(?![^()]*\))', all_each_line):  
                    each_line= each_line.strip()
                    if not each_line.startswith("@@") and  not each_line.startswith("Paper_") :
                        if   re.search("10x|single cell|single-cell|single-nucleus|single nucleus|umap|tsne|suspen|specimen|dissociation|dissociate|scrna|snrnA|dissected|isolat|surgery|model|digest|deposit|were collected|harvest",each_line,re.IGNORECASE)   :
                            if  (Paper_info_dict['Paper_intro'] ).find(each_line) <0 :
                                each_line_full = abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line)), abbreviations_dict=  abbreviations_dict )
                                temp_tissue_line.append(each_line_full)
            Paper_info_dict['Paper_tissue'] = " && " +  " && ".join(list(set(temp_tissue_line)))  + " && "        
            input_file.close()
            input_file = codecs.open(open_file_path,"r","utf-8")
            temp_species_line= []
            for all_each_line in input_file:
                for each_line in re.split(r'\.,\s*(?![^()]*\))', all_each_line):  
                    each_line= each_line.strip()
                    if not each_line.startswith("@@") and  not each_line.startswith("Paper_") :
                        if   re.search("10x|single cell|single-cell|single-nucleus|single nucleus|umap|tsne|suspen|specimen|dissociation|dissociate|scrna|snrnA|dissected|isolat|surgery|model|digest|deposit|were collected|harvest|genome|aligned|obtain",each_line,re.IGNORECASE)   :
                            if  (Paper_info_dict['Paper_intro'] ).find(each_line) <0 :
                                each_line_full = abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_line)), abbreviations_dict=  abbreviations_dict )
                                temp_species_line.append(each_line_full)
            input_file.close()
            Paper_info_dict['Paper_species'] = " && " +  " && ".join(list(set(temp_species_line)))  + " && "    
            if  os.path.exists(re.sub("segment.txt","abstract.tsv",open_file_path)):
                input_file = codecs.open(re.sub("segment.txt","abstract.tsv",open_file_path),"r","utf-8")
                for each_line in input_file:
                    if each_line.startswith("Paper_title") :
                        Paper_info_dict['Paper_title'] =each_line.strip()      
                    elif each_line.startswith("Paper_abstract") :
                        Paper_info_dict['Paper_abstract'] =each_line.strip()     
                    elif each_line.startswith("Paper_keyword") :
                        Paper_info_dict['Paper_keyword'] =each_line.strip() 
                    elif each_line.startswith("Paper_date") :
                        Paper_info_dict['Paper_date'] =each_line.strip() 
                    elif each_line.startswith("Paper_journal") :
                        Paper_info_dict['Paper_journal'] =each_line.strip() 
                input_file.close()
            ##### if  segment was like  :   The trophoblastic clusters (P10–P12) could be annotated as  P10 (EVTBs). it will abbreviation like  the trophoblastic clusters (P10–P12)  could be annotated as extravillous trophoblasts
            # stat_word_line_change_abbreviation =abbreviation2fullname(each_line=stat_word_line, abbreviations_dict=abbreviations_dict)
            keyword_collect = list(set([x.strip() for x in re.split(";|,",re.sub("Paper_keyword:","",Paper_info_dict['Paper_keyword']))]))
            ####  组织解决方案
            tissue_type_output = info_extract(Paper_info_dict=Paper_info_dict, extract_keyword="tissue",model= nlp_model_jnlpba_extra_tissue,candidata_list= keyword_collect,output_file=os.path.dirname(output_file_path)+"/tissue_all_possible.csv")
            ####### 疾病解决方案：
            candidate_disease_list =  []
            for x in keyword_collect:
                for ent in nlp_model_bc5cdr(x).ents:
                        if ent.label_  in ['DISEASE']:
                            candidate_disease_list.append(x.lower())
            candidate_disease_list=list(set(candidate_disease_list))
            disease_name_output = info_extract(Paper_info_dict=Paper_info_dict, extract_keyword="disease",model= nlp_model_bc5cdr,candidata_list= candidate_disease_list,back_list=keyword_collect)
            ###### 物种解决方案 ：
            candidate_species_list =  []
            for x in keyword_collect:
                for ent in nlp_model_craft(x).ents:
                        if ent.label_  in ['TAXON']:
                            candidate_species_list.append(x.upper())
            candidate_species_list = [x  if not p.singular_noun(x)  else p.singular_noun(x) for x in candidate_species_list]
            candidate_species_list = list(set(candidate_species_list) - set(["ANIMAL","MAMMLIAN"]))
            species_name_output = info_extract(Paper_info_dict=Paper_info_dict, extract_keyword="speices",model= nlp_model_craft,candidata_list= candidate_species_list,output_file=os.path.dirname(output_file_path)+"/species_all_possible.csv")
            ##### match paper tissue disease type etc
            input_file = codecs.open(open_file_path,"r","utf-8")
            ##### for each row ,save three sentences wfrom upstream
            context_background_list=list()
            for each_line_org in input_file:
                if each_line_org.startswith("Paper_"):
                    continue
                each_line_org=each_line_org.strip()
                each_line=copy.deepcopy(each_line_org)
                each_line=re.sub("\s\s+"," ",re.sub("\)"," ) ",re.sub("\("," ( ",re.sub("-","-",re.sub("/"," or ",each_line)))))
                if context_background_list==[]:
                    context_background=""   
                else:
                    context_background=". ".join(context_background_list)
                if len(context_background_list)==2:
                    context_background_list.append(each_line)
                    context_background_list=context_background_list[1:4]
                else:
                    context_background_list.append(each_line)
                #end_time = time.time()   
                #print("prepare part : " +str(end_time - start_time) )
                #start_time = time.time()
                each_line_change_abbreviation =abbreviation2fullname(each_line=re.sub("-","-",re.sub("/"," or ",each_line)), abbreviations_dict=  abbreviations_dict )
                ##### match NER model : 
                ######  cell from NER, weight = 2
                ######cell from manual collect, weight =2.6
                ###### segment were not recognized, but have keyword "marker", weight = 1.2                   
                #fullname_entity_collect = list(set([x.text for x in nlp_model_scibert_NER(each_line).ents]  ))
                gene_InOrNot=0
                cell_InOrNot=0
                label_collect=[]
                entity_collect=[]
                temp_cell_collect_1=[]
                temp_cell_collect_2=[]
                temp_cell_collect_3=[]
                temp_cell_collect_4 = []
                gene_entity_output=[]
                cell_entity_output=[]
                for ent in nlp_model_craft(each_line).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)
                for ent in nlp_model_craft(each_line_change_abbreviation).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)       
                for ent in nlp_model_bionlp13cg(each_line).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)      
                for ent in nlp_model_bionlp13cg(each_line_change_abbreviation).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)      
                for ent in nlp_model_jnlpba(each_line).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)      
                for ent in nlp_model_jnlpba(each_line_change_abbreviation).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)   
                for ent in nlp_model_jnlpba_extra_gene_cell(each_line).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)      
                for ent in nlp_model_jnlpba_extra_gene_cell(each_line_change_abbreviation).ents:
                    label_collect.append(ent.label_)
                    entity_collect.append(ent.text)   
                if "CL" in label_collect :
                    cell_InOrNot=1
                    temp_cell_collect_1 = list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CL")     ]))
                    cell_entity_output = cell_entity_output + temp_cell_collect_1
                if "CELL" in label_collect:
                    cell_InOrNot=1
                    temp_cell_collect_2 = list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL")     ]))
                    cell_entity_output = cell_entity_output + temp_cell_collect_2
                if "CELL_TYPE" in label_collect or "CELL_LINE" in label_collect :
                    cell_InOrNot=1
                    temp_cell_collect_3 = []
                    if "CELL_TYPE" in label_collect:
                        temp_cell_collect_3 = list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_TYPE")     ]))
                        cell_entity_output = cell_entity_output + temp_cell_collect_3
                    if "CELL_LINE" in label_collect:
                        temp_cell_collect_3 +=list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
                        cell_entity_output = cell_entity_output + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"CELL_LINE")     ]))
                if "manual_marker" in label_collect:
                    gene_entity_output = gene_entity_output + list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_marker")     ]))
                    gene_entity_output = set(gene_entity_output) - set([x  for x in abbreviations_dict.keys() ])
                    if len(gene_entity_output)>0:
                        gene_InOrNot=1 
                if "manual_cell" in label_collect:
                    cell_InOrNot=1.6
                    temp_cell_collect_4= list(set([ entity_collect[each_index]      for each_index in find_indices(label_collect,"manual_cell")     ]))
                    cell_entity_output +=  temp_cell_collect_4
                if cell_InOrNot+gene_InOrNot>=2:
                    #fullname_entity_collect = list(set(list(set([x.text for x in nlp_model_scibert_NER(each_line_change_abbreviation).ents]  ) ) +list(set([x.text for x in nlp_model_scibert_NER(each_line).ents]  ) )))
                    #fullname_entity_collect = list(set([x.text for x in nlp_model_scibert_NER(each_line_change_abbreviation).ents]  ) )   
                    #for x in  range(0,len(fullname_entity_collect)):
                    #    for y in abbreviations_dict.keys():
                    #        if abbreviations_dict[y] in fullname_entity_collect[x] or fullname_entity_collect[x] in abbreviations_dict[y] :
                    #            fullname_entity_collect[x] =abbreviations_dict[y] 
                    cell_entity_output = [x for x in cell_entity_output if len(x)>2 ]
                    cell_entity_output = [x for x in cell_entity_output if not re.search("\)|\(",x) ]
                    cell_entity_output = [x for x in cell_entity_output if not x.startswith("cell") ]
                    cell_entity_output = [x for x in cell_entity_output if not re.search("single cell",x) and not re.search("single-cell",x) ]
                    cell_entity_output = [   abbreviations_dict[x]  if  x in abbreviations_dict.keys() else x for x in cell_entity_output]
                    cell_entity_output_temp =[]
                    for each_cell_entity_output in cell_entity_output:
                        split_part = each_cell_entity_output.split(" ")
                        cell_entity_output_temp.append(" ".join([abbreviations_dict[x]   if x in abbreviations_dict.keys() else x  for x in split_part]))
                    cell_entity_output =cell_entity_output_temp
                    cell_entity_output = list(set(cell_entity_output))
                    #
                    cell_fullname_entity_output_list = []
                    delete_fullname=[]
                    for  each_cell_entity_output in cell_entity_output:
                        fullname_index = [i for i, s in enumerate(cell_entity_output) if each_cell_entity_output in s]
                        for each_fullname_index in fullname_index:
                            if each_cell_entity_output!= cell_entity_output[each_fullname_index]:
                                if re.search("\S+?\s"+re.escape(each_cell_entity_output), cell_entity_output[each_fullname_index]) or re.search("\s\S+?"+re.escape(each_cell_entity_output), cell_entity_output[each_fullname_index]) or re.search("\S+?\s"+re.escape(each_cell_entity_output)+"\s\S+?", cell_entity_output[each_fullname_index]):
                                    temp_mark =0
                                    for each_temp_list in [temp_cell_collect_1,temp_cell_collect_2,temp_cell_collect_3,temp_cell_collect_4]:
                                        if all(item in each_temp_list  for item in [each_cell_entity_output,cell_entity_output[each_fullname_index]]):
                                            temp_mark=1
                                    if temp_mark==0:
                                        delete_fullname.append(each_cell_entity_output)    
                    cell_fullname_entity_output_list  = list(set(cell_entity_output) -set(delete_fullname))
                    #cell_fullname_entity_output_list = get_fullname_entity(part_entity_list=cell_entity_output, fullname_entity_list=fullname_entity_collect,abbreviations_dict=abbreviations_dict)
                    cell_fullname_entity_output_list = list(set(cell_fullname_entity_output_list)-set(['cells','cell',"Cell",'Cells',"Single-cell","single-cell"]) -gene_entity_output) 
                    cell_fullname_entity_output_list = list(set([ x  if not p.singular_noun(x)  else p.singular_noun(x) for x in cell_fullname_entity_output_list ]))
                    Paper_title = re.sub("Paper_title:","",Paper_info_dict['Paper_title']).strip()
                    Paper_date = re.sub("Paper_date:","",Paper_info_dict['Paper_date']).strip()
                    Paper_journal = re.sub("Paper_journal:","",Paper_info_dict['Paper_journal']).strip()
                    if len(cell_fullname_entity_output_list)>0:
                        print( paper_id+"\t"+ Paper_title+"\t"+ Paper_date+"\t"+Paper_journal +"\t"+species_name_output +"\t"+disease_name_output+"\t"+tissue_type_output +"\t" +each_line_org+"\t"+str(cell_InOrNot+gene_InOrNot)+"\t"+" & ".join(list(gene_entity_output))+"\t"+" & ".join(list(set(cell_fullname_entity_output_list)))+"\t"+context_background, file=output )
            input_file.close()                
            output.close()
            output_fig.close()
    else:
        print("file do not exist : "+open_file_path)


def is_float_digit(n: str) -> bool:
    try:
        float(n)
        return True
    except ValueError:
        return False


def is_gene(n: str,use_gene_collect) -> bool:
     try:
        temp_str = str(n)
        if temp_str in use_gene_collect:
            return True
        else :
            return False
     except ValueError:
         return False

def is_fig_anno(n: str) -> bool:
    try:
        temp_str = str(n)
        pattern1 = re.compile("^([0-9]+[a-z]+)")
        pattern2 = re.compile("^(s[0-9]+[a-z]+)")
        if temp_str.startswith("fig") or temp_str.startswith("Fig") or temp_str.startswith("Extend") or temp_str.startswith("extend") or temp_str.startswith("supp") or temp_str.startswith("[") or  pattern1.match(temp_str)!=None or  pattern2.match(temp_str)!=None:
            return True
        else:
            return False
    except ValueError:
        return False


def preprocessing(text,use_gene_collect,abbr_dict_path):
    if os.path.exists(abbr_dict_path):
        abbreviations_dict = pickle.load(open(abbr_dict_path, 'rb'))
    else:
        abbreviations_dict=dict()
    try:
        text = abbreviation2fullname(re.sub("\s\s+"," ",re.sub('\)'," ) ",re.sub('\('," ( ",re.sub('\[.*?\]',"",text)))), abbreviations_dict)
        all_combine_cell_cluster =[]
        for i in range(1,10):
            for cell_name in abbreviations_dict.keys():
                all_combine_cell_cluster.append(cell_name+str(i))
        each_file_use_gene_collect = use_gene_collect - set(all_combine_cell_cluster )
        tokens = [token.lemma_ if is_gene(token,use_gene_collect=each_file_use_gene_collect)==False else str(token) for token in nlp_model_scibert(text)]
        tokens = [t for t in tokens if 
                    t not in STOP_WORDS and 
                    t not in punctuation_word and 
                    t.strip()!="" and 
                    t.isnumeric()==False and 
                    is_float_digit(t)==False and 
                    is_fig_anno(t)==False ] 
        tokens = [t.lower() if is_gene(t,use_gene_collect=each_file_use_gene_collect)==False else t for t in tokens  ] 
        tokens_remove_meanless_parentheses =tokens
        istart = []  # stack of indices of opening parentheses
        d = {}
        for i, c in enumerate(tokens_remove_meanless_parentheses):
            if c == '(':
                istart.append(i)
            if c == ')':
                try:
                    d[istart.pop()] = i
                except IndexError:
                    print('Too many closing parentheses')
        #
        #
        delete_d = []
        for each_d in d.keys():
            for each_d_temp in d.keys():
                if int(each_d)>int(each_d_temp):
                    if int(d[each_d] <d[each_d_temp] ):
                        delete_d.append(each_d)
        #
        #
        for each_d_use in d.keys():
            if  each_d_use not in delete_d:
                if  len(set(each_file_use_gene_collect).intersection(tokens_remove_meanless_parentheses[each_d_use :  (d[each_d_use]+1) ] ))>0:
                    tokens_remove_meanless_parentheses[(each_d_use+1) :  (d[each_d_use]) ] =[""  for x in tokens_remove_meanless_parentheses[(each_d_use+1) :  (d[each_d_use]) ] ]
                else:
                    tokens_remove_meanless_parentheses[(each_d_use) :  (d[each_d_use]+1) ] =[""  for x in tokens_remove_meanless_parentheses[(each_d_use) :  (d[each_d_use]+1) ] ]
#
        tokens_remove_meanless_parentheses = []
        parentheses_index=0
        temp_token=[]
        front_word_content = ""
        for each_tokens in tokens:
            if each_tokens!="(" and each_tokens!=")" :
                if parentheses_index!=1 :
                    front_word_content=each_tokens
                    if each_tokens not in each_file_use_gene_collect:
                        tokens_remove_meanless_parentheses.append(each_tokens)
                else:
                    temp_token.append(each_tokens)
            elif each_tokens=="(":
                temp_token=[]
                parentheses_index=1
            elif each_tokens==")":
                if len( set(temp_token) & each_file_use_gene_collect )>0 :
                    tokens_remove_meanless_parentheses.append("(")
                    #for x in temp_token:
                        #tokens_remove_meanless_parentheses.append(x)
                    tokens_remove_meanless_parentheses.append(")")
                elif  front_word_content in each_file_use_gene_collect and  temp_token!=[]:
                    tokens_remove_meanless_parentheses.append("(")
                    tokens_remove_meanless_parentheses.append(")")
                each_tokens=[]
                parentheses_index=0 
                #print(temp_token)
            # print(tokens_remove_meanless_parentheses)       
        return(re.sub("\s\s+"," ",re.sub("\[*\d.*?\]"," "," ".join(tokens_remove_meanless_parentheses))))
    except:
        return("")


def add_segment_clean(open_file_path,use_gene_collect,abbr_dict_path):
    output_file_path = re.sub("evidence","evidence_clean",open_file_path)
    if os.path.exists(output_file_path) and os.stat(output_file_path).st_size>1024:
        print(output_file_path + " : has exist")
    elif not os.path.exists(open_file_path):
        print(open_file_path + " : not exist")
    else:
        input_data = pd.read_table(open_file_path, delimiter="\t", quoting=csv.QUOTE_NONE)
        temp = input_data["Evidence"].apply(preprocessing,use_gene_collect=use_gene_collect,abbr_dict_path=abbr_dict_path)
        input_data["Evidence_clean"] = temp
        input_data.to_csv(output_file_path,sep="\t",index=False,encoding='utf_8_sig')


def add_segment_text_clean(open_file_path,abbr_dict_path):
    if not os.path.exists(open_file_path) :
        print(open_file_path + " : not exist")
    else:
        gene_collect = pd.read_csv(gene_collect_path)
        use_gene_collect = set([str(x) for x in list(gene_collect['gene_collect'])])
        print(open_file_path+" : start add clean evidence")
        add_segment_clean(open_file_path=open_file_path,use_gene_collect=use_gene_collect,abbr_dict_path=abbr_dict_path)
        print(open_file_path +" : done ")


def model_result(x,keyword):
    temp=nlp_model(str(x[keyword])).cats
    if temp['True']>temp['False']:
        return("True")
    else:
        return("False")


def model_result_true(x,keyword):
    temp=nlp_model(str(x[keyword])).cats
    return(temp['True'])


def model_result_false(x,keyword):
    temp=nlp_model(str(x[keyword])).cats
    return(temp['False'])


def predict_textcat(paper_id_path,work_path,input_file_name):
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    for file in tqdm(paper_id_file):
        file= str(int(file))
        open_file_path = work_path+"/paper_process_result/"+file+"/"+input_file_name
        input_file = pd.read_table(open_file_path,encoding='utf-8',delimiter="\t",quoting=csv.QUOTE_NONE)
        Predict_result= input_file.apply( model_result,axis=1)
        input_file['Predict_result'] =Predict_result
        Predict_result_true_p= input_file.apply( model_result_true,axis=1)
        input_file['Predict_result_true_p'] =Predict_result
        Predict_result_false_p= input_file.apply( model_result_false,axis=1)
        input_file['Predict_result_false_p'] =Predict_result        
        input_file.to_table(work_path+"/paper_process_result/"+file+"/predict_textcat.txt")


def find_indices(list_to_check,item_to_find):
    return [idx for idx, value in enumerate(list_to_check) if value==item_to_find]


def get_fullname_entity(part_entity_list, fullname_entity_list,abbreviations_dict):
    output_list=list()
    for each_entity_output in part_entity_list:
        for  each_fullname_entity_collect in fullname_entity_list:
            split_each_entity_output = each_entity_output.split(" ")
            if len(split_each_entity_output) >1 :
                for each_split_each_entity_output in split_each_entity_output:
                    if each_split_each_entity_output in abbreviations_dict.keys():
                        temp_entity = abbreviations_dict[each_split_each_entity_output]
                    else:
                        temp_entity = each_split_each_entity_output
                    if len(re.findall('[A-z]', temp_entity))>1 and  re.match(".*\s"+re.escape(temp_entity) ,each_fullname_entity_collect ) or re.match(re.escape(temp_entity)+"\s.*"  ,each_fullname_entity_collect ):
                        output_list.append(each_fullname_entity_collect)
            else:
                if len(re.findall('[A-z]', each_entity_output))>1 and str(each_fullname_entity_collect).startswith(each_entity_output) or str(each_fullname_entity_collect).endswith(each_entity_output):
                    output_list.append(each_fullname_entity_collect)
    return(list(set(output_list)))


def model_result_true(x,keyword,nlp_model):
    temp=nlp_model(str(x[keyword])).cats
    return(temp['True'])


def add_predict_result(input_filepath ,output_filepath, nlp_model,keyword="Evidence_clean"):
    if os.path.exists(output_filepath) and os.stat(output_filepath).st_size>1024:
        print(output_filepath + " : has exist")
    elif not os.path.exists(input_filepath):
        print(input_filepath + " : not exist")
    else:
        input_data = pd.read_table(input_filepath,encoding='ISO-8859-1' ,delimiter="\t")
        input_data.insert(loc=input_data.shape[1], column='Predict_prob', value=input_data.apply( model_result_true,keyword=keyword,nlp_model=nlp_model,axis=1) )
        input_data = input_data.sort_values(by="Predict_prob",ascending=False)
        input_data.to_csv(output_filepath,index=False,sep="\t",encoding='ISO-8859-1')


def punct_segment_extract(doc):
    matcher = Matcher(nlp_model_scibert.vocab)
    pattern = [ {'DEP':'amod', 'OP':'{0,5}'}, 
                {'DEP':'nmod:npmod', 'OP':'{0,5}'}, 
                {'DEP':'amod', 'OP':'{0,5}'}, 
                {'DEP':'compound', 'OP':'{0,5}'}, 
                {'DEP':'amod', 'OP':'{0,5}'}, 
                {'POS':'NOUN', 'OP':'{0,2}'}, # 形容词修饰
                {'POS':'ADJ', 'OP':'{0,1}'}, # 形容词修饰
                {'DEP':'nmod', 'OP':'{0,1}'}, # 形容词修饰
                {'POS':'NOUN', 'OP':'{0,1}'}, # 形容词修饰
                {'TEXT':'('}]
    matcher.add("segment_with_punct", [pattern])
    matches = matcher(doc)
    segment_with_punct_pattern = []
    for i in range(0,len(matches)):
        segment_with_punct_pattern.append(str(doc[matches[i][1]:matches[i][2]]))
    delete_index =[]
    for index in range(1,len(segment_with_punct_pattern)):
        if re.search(re.escape(segment_with_punct_pattern[-index]),segment_with_punct_pattern[-index-1]):
            delete_index.append(-index)
    #from operator import itemgetter
    if len(delete_index)>0:
        segment_with_punct_pattern = list( set(segment_with_punct_pattern)  -set([segment_with_punct_pattern[i] for i in delete_index])   )
    output_dataframe = pd.DataFrame()
    for each_segment_with_punct_pattern in segment_with_punct_pattern:
        true_index_instr = doc.text.index(each_segment_with_punct_pattern)
        all_index = [m.start() for m in re.finditer(re.escape(re.sub("\s*\(","",each_segment_with_punct_pattern)), doc.text)]
        if len(all_index)==0:
            continue
        true_index =all_index.index(true_index_instr)
        count_index=0
        p1 = re.compile(re.escape(each_segment_with_punct_pattern)+r'(.*?)[)]', re.S)
        str_in_punct = re.findall(p1, doc.text)
        if len(str_in_punct)==0:
            continue
        split_str_in_punct = [x.strip() for x in re.split(",|and|/" , str_in_punct[0])]
        for i,token in enumerate(doc):
            if re.sub("\s*\(","",each_segment_with_punct_pattern).endswith(token.text) and token.idx >= true_index_instr  and (token.idx -  true_index_instr)<len(each_segment_with_punct_pattern) :
                if count_index !=true_index:
                    true_index+=1
                else:
                 #   if token.head.pos_ =="NOUN"   :
                #        temp_dataframe = pd.DataFrame({"left_col":token.head.text, "right_col" :str_in_punct})
                  #      output_dataframe = pd.concat([output_dataframe,temp_dataframe])
                 #   for each_child in token.head.children:
                  #      if each_child.pos_ =="NOUN" and doc.text.index(each_child.text)<doc.text.index(each_segment_with_punct_pattern)  :
                  #          temp_dataframe = pd.DataFrame({"left_col":each_child.text, "right_col" :str_in_punct})
                  #          output_dataframe = pd.concat([output_dataframe,temp_dataframe])
                    temp_dataframe = pd.DataFrame({"left_col":re.sub("@@.*","",fill_token_name(token)['full_name'][0]), "right_col" :split_str_in_punct,"distance":1,"same_ancestors":1,"subtree_number":0})
                    output_dataframe = pd.concat([output_dataframe,temp_dataframe])
      #  temp_dataframe = pd.DataFrame({"left_col":re.sub("\s*\(","",each_segment_with_punct_pattern), "right_col" :str_in_punct})
      #  output_dataframe = pd.concat([output_dataframe,temp_dataframe])   
    output_dataframe = output_dataframe.drop_duplicates()
    return(output_dataframe)


def another_subtree_judge_rootlevel(temp_tok):
    return_index=0
    if temp_tok.dep_ in "conj" and temp_tok.n_rights+temp_tok.n_lefts>0:
        for x in temp_tok.subtree:
            if x.dep_ in  ['nsubj',"nsubjpass","nmod","nmod:npmod","appos","dep"]:
                return_index=1
    return(return_index)


def another_subtree_judge_sublevel(temp_tok):
    return_index=0
    if  temp_tok.n_rights+temp_tok.n_lefts>0:
        for x in temp_tok.subtree:
            if x!=temp_tok and x.dep_ in  ['nsubj',"nsubjpass","nmod","nmod:npmod","appos","dep"]:
                return_index=1
    return(return_index)


def complete_subtree_judge(temp_tok):
    return_index=0
    if  temp_tok.n_rights+temp_tok.n_lefts>0:
        for x in temp_tok.subtree:
            if x.dep_ in  ['nsubj',"nsubjpass","nmod","nmod:npmod","appos","dep"]:
                return_index=1
        if temp_tok.dep_ =="conj":
                return_index=2
    return(return_index)


def subtree_all(doc):
    output_dataframe=pd.DataFrame()
    puzzle_datafrme=pd.DataFrame()
    puzzle_keyword = "|".join([" small "," non "," not ","non-","respectively"," neither ","Neither ","distinguish"," nor "," no "," no-",".*?:.*?;.*?:.*"])
    pass_keyword=[ 'decreased',"decrease","downregulated","downregulate","weak"," absent"," absented"," low "," lack "," lacks"," lower "," few "]
    if  len(re.findall(puzzle_keyword,doc.text ))>0:
        temp_dataframe = pd.DataFrame({"Puzzle_evidence":[doc.text]})     
        puzzle_datafrme = pd.concat([puzzle_datafrme,temp_dataframe])     
        return({"output_dataframe":output_dataframe,"puzzle_datafrme":puzzle_datafrme})
    for i,tok in enumerate(doc):
        #if tok.dep_.endswith("obj") or tok.head.text==tok.text or tok.dep_ in ['nsubj',"nsubjpass","nmod","csubj","csubjpass","conj","ccomp","appos"] : 更少更准
        if tok.text in pass_keyword:
            continue
        if tok.dep_.endswith("obj") or tok.head.text==tok.text or tok.dep_ in ['nsubj',"nsubjpass","nmod","csubj","csubjpass","conj","ccomp"] :
            dict_dict=dict()
            level_dict=dict()
            complete_subtree_dict =dict()
            ancestors_dict=dict()
            temp = fill_token_name(tok)
            token_name = temp["full_name"]
            dict_dict['root_name']=token_name
            level_dict[token_name[0]]=0
            complete_subtree_dict[token_name[0]]=0
            complete_subtree_dict[tok.text+"@@"+str(tok.idx)]=0
            ancestors_dict[token_name[0]]=temp['ancestors']
            root_child = [x    for x in tok.children ]
            root_child_text = [x.text   for x in tok.children ] 
            for each_erroneous_judgemen_word in pass_keyword:
                if each_erroneous_judgemen_word in root_child_text:
                    root_child = root_child[:root_child_text.index(each_erroneous_judgemen_word)]   
            root_child = [ x  for x in root_child if x not in temp['delete_child']    ]
            for x in temp['delete_child']:
                complete_subtree_dict[x.text+"@@"+str(x.idx)]=complete_subtree_dict[tok.text+"@@"+str(tok.idx)]
            for each_root_child in root_child:
                #pass_or_not = another_subtree_judge_rootlevel(each_root_child)
                #if pass_or_not==1:
                #    continue
                complete_subtree_dict[each_root_child.text+"@@"+str(each_root_child.idx)] = complete_subtree_judge(each_root_child)
                temp_list=[]
                return_dict_org= make_level_dict(input_token=each_root_child, current_level=0,level_dict=level_dict,ancestors_dict=ancestors_dict,complete_subtree_dict=complete_subtree_dict,temp_list=temp_list,father_tok_text=token_name[0])
                temp = fill_token_name(each_root_child)
                token_name = temp["full_name"]
                level_dict = return_dict_org['level_dict']
                ancestors_dict = return_dict_org['ancestors_dict']
                complete_subtree_dict = return_dict_org['complete_subtree_dict']
                temp_list  = return_dict_org['temp_list']
                if complete_subtree_dict[each_root_child.text+"@@"+str(each_root_child.idx)]>1:
                    continue
                if each_root_child.n_rights+each_root_child.n_lefts>0:
                    sub1_child = [x    for x in each_root_child.children ]   
                    sub1_child_text = [x.text   for x in each_root_child.children ]  
                    for each_erroneous_judgemen_word in pass_keyword:
                        if each_erroneous_judgemen_word in sub1_child_text:
                            sub1_child = sub1_child[:sub1_child_text.index(each_erroneous_judgemen_word)]                   
                    sub1_child = [ x  for x in sub1_child if x not in temp['delete_child']    ]
                    for x in temp['delete_child']:
                        complete_subtree_dict[x.text+"@@"+str(x.idx)]=complete_subtree_dict[each_root_child.text+"@@"+str(each_root_child.idx)]
                    for each_sub1_child in sub1_child:
                        #if another_subtree_judge_sublevel(each_sub1_child)==1 or each_sub1_child.dep_=="dep" or each_sub1_child.dep_=="appos":   更少更准
                        complete_subtree_dict[each_sub1_child.text+"@@"+str(each_sub1_child.idx)] = complete_subtree_judge(each_sub1_child)
                       # if another_subtree_judge_sublevel(each_sub1_child)==1 or each_sub1_child.dep_=="dep" :
                           # continue
                        return_dict= make_level_dict(input_token=each_sub1_child, current_level=1,level_dict=level_dict,ancestors_dict=ancestors_dict,complete_subtree_dict=complete_subtree_dict,temp_list=temp_list,father_tok_text=token_name[0])
                        temp = fill_token_name(each_root_child)
                        token_name = temp["full_name"]
                        level_dict = return_dict['level_dict']
                        ancestors_dict = return_dict['ancestors_dict']
                        complete_subtree_dict = return_dict['complete_subtree_dict']
                        temp_list  = return_dict['temp_list']
                        if complete_subtree_dict[each_sub1_child.text+"@@"+str(each_sub1_child.idx)]>1:
                            continue
                        if each_sub1_child.n_rights+each_sub1_child.n_lefts>0:
                            sub2_child = [x    for x in each_sub1_child.children ]   
                            sub2_child_text = [x.text   for x in each_sub1_child.children ]  
                            for each_erroneous_judgemen_word in pass_keyword:
                                if each_erroneous_judgemen_word in sub2_child_text:
                                    sub2_child = sub2_child[:sub2_child_text.index(each_erroneous_judgemen_word)]  
                            sub2_child = [ x  for x in sub2_child if x not in temp['delete_child']    ]  
                            for x in temp['delete_child']:
                                complete_subtree_dict[x.text+"@@"+str(x.idx)]=complete_subtree_dict[each_sub1_child.text+"@@"+str(each_sub1_child.idx)]           
                            for each_sub2_child in sub2_child:
                                # if another_subtree_judge_sublevel(each_sub2_child)==1 or each_sub2_child.dep_=="dep" or each_sub2_child.dep_=="appos":   ##更少更准
                                #if another_subtree_judge_sublevel(each_sub2_child)==1 or each_sub2_child.dep_=="dep" :  
                                #    continue
                                complete_subtree_dict[each_sub2_child.text+"@@"+str(each_sub2_child.idx)] = complete_subtree_judge(each_sub2_child)
                                return_dict1= make_level_dict(input_token=each_sub2_child, current_level=2,level_dict=level_dict,ancestors_dict=ancestors_dict,complete_subtree_dict=complete_subtree_dict,temp_list=temp_list,father_tok_text=token_name[0])
                                temp = fill_token_name(each_sub2_child)
                                token_name = temp["full_name"]
                                level_dict = return_dict1['level_dict']
                                ancestors_dict = return_dict1['ancestors_dict']
                                complete_subtree_dict = return_dict1['complete_subtree_dict']
                                temp_list  = return_dict1['temp_list']
                                if complete_subtree_dict[each_sub2_child.text+"@@"+str(each_sub2_child.idx)]>1:
                                    continue
                                if each_sub2_child.n_rights+each_sub2_child.n_lefts>0:
                                    sub3_child = [x    for x in each_sub2_child.children ]   
                                    sub3_child_text = [x.text   for x in each_sub2_child.children ]      
                                    for each_erroneous_judgemen_word in pass_keyword:
                                        if each_erroneous_judgemen_word in sub3_child_text:
                                            sub3_child = sub3_child[:sub3_child_text.index(each_erroneous_judgemen_word)]         
                                    sub3_child = [ x  for x in sub3_child if x not in temp['delete_child']    ] 
                                    for x in temp['delete_child']:
                                        complete_subtree_dict[x.text+"@@"+str(x.idx)]=complete_subtree_dict[each_sub2_child.text+"@@"+str(each_sub2_child.idx)]   
                                    for each_sub3_child in sub3_child:
                                        # if another_subtree_judge_sublevel(each_sub3_child)==1 or each_sub3_child.dep_=="dep" or each_sub3_child.dep_=="appos":#更少更准
                                        #if another_subtree_judge_sublevel(each_sub3_child)==1 or each_sub3_child.dep_=="dep" :
                                        #    continue
                                        complete_subtree_dict[each_sub3_child.text+"@@"+str(each_sub3_child.idx)] = complete_subtree_judge(each_sub3_child)
                                        return_dict2= make_level_dict(input_token=each_sub3_child, current_level=3,level_dict=level_dict,ancestors_dict=ancestors_dict,complete_subtree_dict=complete_subtree_dict,temp_list=temp_list,father_tok_text=token_name[0])
                                        temp = fill_token_name(each_sub3_child)
                                        token_name = temp["full_name"]
                                        level_dict = return_dict2['level_dict']
                                        ancestors_dict = return_dict2['ancestors_dict']
                                        complete_subtree_dict = return_dict2['complete_subtree_dict']
                                        temp_list  = return_dict2['temp_list']
                                        if complete_subtree_dict[each_sub3_child.text+"@@"+str(each_sub3_child.idx)]>1:
                                            continue
                                        if each_sub3_child.n_rights+each_sub3_child.n_lefts>0:
                                            sub4_child = [x    for x in each_sub3_child.children ]   
                                            sub4_child_text = [x.text   for x in each_sub3_child.children ]    
                                            for each_erroneous_judgemen_word in pass_keyword:
                                                if each_erroneous_judgemen_word in sub4_child_text:
                                                    sub4_child = sub4_child[:sub4_child_text.index(each_erroneous_judgemen_word)]    
                                            sub4_child = [ x  for x in sub4_child if x not in temp['delete_child']    ]  
                                            for x in temp['delete_child']:
                                                complete_subtree_dict[x.text+"@@"+str(x.idx)]=complete_subtree_dict[each_sub3_child.text+"@@"+str(each_sub3_child.idx)]                                                          
                                            for each_sub4_child in sub4_child:
                                                #if another_subtree_judge_sublevel(each_sub4_child)==1 or each_sub4_child.dep_=="dep" or each_sub4_child.dep_=="appos": 更少更准
                                                #if another_subtree_judge_sublevel(each_sub4_child)==1 or each_sub4_child.dep_=="dep" :
                                                #    continue
                                                complete_subtree_dict[each_sub4_child.text+"@@"+str(each_sub4_child.idx)] = complete_subtree_judge(each_sub4_child)
                                                return_dict4= make_level_dict(input_token=each_sub4_child, current_level=4,level_dict=level_dict,ancestors_dict=ancestors_dict,complete_subtree_dict=complete_subtree_dict,temp_list=temp_list,father_tok_text=token_name[0])
                                                temp = fill_token_name(each_sub4_child)
                                                token_name = temp["full_name"]
                                                level_dict = return_dict4['level_dict']
                                                ancestors_dict = return_dict4['ancestors_dict']
                                                complete_subtree_dict = return_dict4['complete_subtree_dict']
                                                temp_list  = return_dict4['temp_list']
                dict_dict[each_root_child] = temp_list
                if len(dict_dict)>50:
                    continue
                #### make output 
                dict_key=list(dict_dict.keys())
                for each_key_index in range(0,(len(dict_dict)-1)):
                    for each_key_index_use in range(each_key_index+1, min(each_key_index+3 ,len(dict_dict )  )):
                        for temp_i_index in range(0,len(dict_dict[dict_key[each_key_index]])):
                            for temp_j_index in range(0,len(dict_dict[dict_key[each_key_index_use]])):
                                temp_i = dict_dict[dict_key[each_key_index]][temp_i_index]
                                temp_j = dict_dict[dict_key[each_key_index_use]][temp_j_index]
                                if abs(level_dict[temp_j]-level_dict[temp_i])<4 :
                                    if  all(item in ancestors_dict[temp_i] for item in ancestors_dict[temp_j]) or all(item in ancestors_dict[temp_j] for item in ancestors_dict[temp_i]) :
                                        temp_dataframe = pd.DataFrame({"left_col":[re.sub("@@.*","",temp_i)], "right_col" :[re.sub("@@.*","",temp_j)],"distance":abs(level_dict[temp_j]-level_dict[temp_i]),"same_ancestors":1,"subtree_number":abs(complete_subtree_dict[temp_j]+complete_subtree_dict[temp_i])})     
                                        output_dataframe = pd.concat([output_dataframe,temp_dataframe])  
                                    else:
                                        temp_dataframe = pd.DataFrame({"left_col":[re.sub("@@.*","",temp_i)], "right_col" :[re.sub("@@.*","",temp_j)],"distance":abs(level_dict[temp_i]-level_dict[temp_i]),"same_ancestors":0,"subtree_number":abs(complete_subtree_dict[temp_j]+complete_subtree_dict[temp_i])})     
                                        output_dataframe = pd.concat([output_dataframe,temp_dataframe])  
        elif  tok.dep_=="dep" :
            if not tok.head.dep_.endswith("obj") or tok.head.head.text==tok.head.text or tok.head.dep_ in ['nsubj',"nsubjpass","nmod","csubj","csubjpass","conj"]:
                temp_dataframe = pd.DataFrame({"left_col":[re.sub("@@.*","",fill_token_name(tok.head)['full_name'][0])], "right_col" :[tok.text],"distance":1,"same_ancestors":1,"subtree_number":0})     
                output_dataframe = pd.concat([output_dataframe,temp_dataframe])      
                if tok.n_rights+tok.n_lefts>0:
                    for each_child in tok.children:
                        if each_child.dep_ in ["dep","conj","appos"]:
                            temp_dataframe = pd.DataFrame({"left_col":[re.sub("@@.*","",fill_token_name(tok.head)['full_name'][0])], "right_col" :[each_child.text],"distance":1,"same_ancestors":1,"subtree_number":0})     
                            output_dataframe = pd.concat([output_dataframe,temp_dataframe])                       
    return({"output_dataframe":output_dataframe,"puzzle_datafrme":puzzle_datafrme})


def make_level_dict(input_token, current_level,level_dict,ancestors_dict,complete_subtree_dict,temp_list,father_tok_text=None):
    temp = fill_token_name(input_token)
    if input_token.dep_ not in ["conj", "appos"] :
        temp_list.append(temp["full_name"][0])
        complete_subtree_dict[temp["full_name"][0]] = complete_subtree_dict[input_token.text+"@@"+str(input_token.idx)]
        level_dict[temp["full_name"][0]]=current_level+1
        ancestors_dict[temp["full_name"][0]]=temp['ancestors']
    elif input_token.dep_ =="conj":
        if (input_token.n_rights+input_token.n_lefts)==0 :
            level_dict[temp["full_name"][0]]=current_level
            ancestors_dict[temp["full_name"][0]]=ancestors_dict[father_tok_text]
            temp_list.append(temp["full_name"][0])
            complete_subtree_dict[temp["full_name"][0]] = complete_subtree_dict[father_tok_text]
        else:
            level_dict[temp["full_name"][0]]=current_level+1
            ancestors_dict[temp["full_name"][0]]=temp['ancestors']
            complete_subtree_dict[temp["full_name"][0]] = complete_subtree_dict[father_tok_text]
    elif input_token.dep_ =="appos":
        level_dict[temp["full_name"][0]]=level_dict[father_tok_text]
        ancestors_dict[temp["full_name"][0]]=ancestors_dict[father_tok_text]
        temp_list.append(temp["full_name"][0])
        complete_subtree_dict[temp["full_name"][0]] = complete_subtree_dict[father_tok_text]
    return_dict={"temp_list":temp_list,"level_dict":level_dict, "ancestors_dict":ancestors_dict,"complete_subtree_dict":complete_subtree_dict }
    return(return_dict)


def fill_token_name(tok):
    return_name = tok.text+"@@"+str(tok.idx)
    delete_child=[]
    ancestors=[ x.text+"@@"+str(x.idx) for x in tok.ancestors]
    output_form_marker=4
    if tok.n_rights+tok.n_lefts >0:
        sub_name=""
        for each_tok_child in tok.children:
            if each_tok_child.dep_ in ['amod', "nummod", "compound"]:
                temp_child_text = each_tok_child.text.strip()
                delete_child.append(each_tok_child)
                if each_tok_child.n_rights+each_tok_child.n_lefts >0:
                    temp_child_text_amod=""
                    for each_tok_child_child in each_tok_child.children:
                        if each_tok_child_child.dep_ in ['amod', "nummod", "compound"]:
                            temp_child_text_amod= temp_child_text_amod+" "+each_tok_child_child.text
                            delete_child.append(each_tok_child_child)
                    temp_child_text = temp_child_text_amod+" "+temp_child_text    
                sub_name =  sub_name+" " +temp_child_text
                #delete_child.append(temp_child_text)
                output_form_marker=0
         #   elif each_tok_child.dep_ in ['conj']:
          #      temp_child_text = fill_token_name(each_tok_child)
          #      sub_name = sub_name+"," +temp_child_text['full_name'][0]
          #      delete_child.append(each_tok_child.text)
          #     output_form_marker=1
            elif each_tok_child.dep_ in ['cop',"punct","cc","det",'mark',"det","case"] :
                temp_child_text = each_tok_child.text.strip()
                delete_child.append(each_tok_child)
    if output_form_marker==0:
            return_name = sub_name.strip()+" " +return_name
    elif output_form_marker==1:
            return_name = sub_name.strip()+"," +return_name
    return({"full_name":[return_name], "delete_child":delete_child,"ancestors":ancestors})


def cell_marker_extract_main(file_path,output_path,abbreviations_dict_path_root,keyword):
    if not os.path.exists(file_path):
        print(file_path + " : not exist")
    else:
        segment_data =pd.read_table(file_path, delimiter="\t", quoting=csv.QUOTE_NONE)
        use_segment_data =segment_data[segment_data['Predict_prob']>0.7]
        output_pd =pd.DataFrame()
        puzzle_output_pd =pd.DataFrame()
        for i in tqdm(range(0,use_segment_data.shape[0])):
            Segment_line = use_segment_data[[keyword]].iloc[i,0]
            if Segment_line.startswith("Paper_"):
                pass
            pmid = use_segment_data[['PMID']].iloc[i,0]
            abbreviations_dict_path = abbreviations_dict_path_root+"/"+str(pmid)+"/abbreviations_dict.pickle"
            if os.path.exists(abbreviations_dict_path):
                abbreviations_dict = pickle.load(open(abbreviations_dict_path, 'rb'))
            else:
                abbreviations_dict=dict()
            if isinstance( use_segment_data[['Gene_name']].iloc[i,0],float) or isinstance( use_segment_data[['Cell_name']].iloc[i,0],float):
                continue
            gene_name_collect = [x.strip() for x in use_segment_data[['Gene_name']].iloc[i,0].split("&")]
            cell_name_collect = [x.strip() for x in use_segment_data[['Cell_name']].iloc[i,0].split("&")]
            #cell_name_collect = cell_name_collect +   [x    for x in abbreviations_dict.keys() if abbreviations_dict[x] in  cell_name_collect ]
            cell_name_collect = [x for x in cell_name_collect if not re.search("\(|\)",x)]
            for each_Segment_line in re.split(r';\s*(?![^()]*\))', Segment_line):
                for temp_each_Segment_line in re.split("while|whereas", each_Segment_line):
                    if len(temp_each_Segment_line)<5 or len(re.findall('[A-z]', temp_each_Segment_line))<(len(re.findall('\S', temp_each_Segment_line))*0.5):
                        continue
                    for x in gene_name_collect: 
                        if re.search(x+"\s{0,1}\^{0,1}\-|\−",temp_each_Segment_line):
                            temp_each_Segment_line= re.sub(x+"\s{0,1}\^{0,1}\-|\−","",temp_each_Segment_line)
                    return_dict = extract_cell_marker(each_Segment_line=temp_each_Segment_line,gene_name_collect=gene_name_collect, cell_name_collect=cell_name_collect,abbreviations_dict=abbreviations_dict)
                    temp_oututdata=pd.DataFrame()
                    temp_oututdata=return_dict['temp_oututdata']
                    if temp_oututdata.shape[0]>0:
                        temp_oututdata[['PMID']]=pmid
                        temp_oututdata[['Title']]=use_segment_data[['Title']].iloc[i,0]
                        temp_oututdata[['Date']]=use_segment_data[['Date']].iloc[i,0]
                        temp_oututdata[['Journal']]=use_segment_data[['Journal']].iloc[i,0]
                        temp_oututdata[['Species']]=use_segment_data[['Species']].iloc[i,0]
                        temp_oututdata[['Tissue']]=use_segment_data[['Tissue_type']].iloc[i,0]
                        temp_oututdata[['Disease']]=use_segment_data[['Disease']].iloc[i,0]
                        temp_oututdata[[keyword]]=Segment_line
                        temp_oututdata[["Evidence_prob"]]=use_segment_data[['Predict_prob']].iloc[i,0]
                        temp_oututdata[["context_background"]]=use_segment_data[['context_background']].iloc[i,0]
                        temp_oututdata = temp_oututdata[temp_oututdata['same_ancestors']==1  ]
                        temp_oututdata= temp_oututdata[['PMID','Title','Date','Journal',"Species","Tissue","Disease",keyword,"Evidence_prob",'Org_Cell_type',"decorate_Cell_type",'Marker',"distance","subtree_number","context_background"]]
                        temp_oututdata =temp_oututdata.drop_duplicates()
                        output_pd = pd.concat([temp_oututdata,output_pd])
                    puzzle_dataframe=pd.DataFrame()
                    puzzle_dataframe=return_dict['puzzle_dataframe']
                    if puzzle_dataframe.shape[0]>0:
                        puzzle_dataframe[['PMID']]=pmid
                        puzzle_dataframe[['Title']]=use_segment_data[['Title']].iloc[i,0]
                        puzzle_dataframe[['Date']]=use_segment_data[['Date']].iloc[i,0]
                        puzzle_dataframe[['Journal']]=use_segment_data[['Journal']].iloc[i,0]
                        puzzle_dataframe[['Species']]=use_segment_data[['Species']].iloc[i,0]
                        puzzle_dataframe[['Tissue']]=use_segment_data[['Tissue_type']].iloc[i,0]
                        puzzle_dataframe[['Disease']]=use_segment_data[['Disease']].iloc[i,0]
                        puzzle_dataframe[['Evidence_prob']]=use_segment_data[['Predict_prob']].iloc[i,0]
                        puzzle_dataframe= puzzle_dataframe[['PMID','Title','Date','Journal',"Species","Tissue","Disease","Puzzle_evidence","Evidence_prob"]]     
                    puzzle_output_pd = pd.concat([puzzle_dataframe,puzzle_output_pd])
        summary_output = output_pd.drop_duplicates()
        summary_output.to_csv(output_path,index=False,na_rep='',encoding='utf_8_sig')
        puzzle_output_pd.drop_duplicates().to_csv(re.sub("\\.","_puzzle.",output_path),index=False,na_rep='',encoding='utf_8_sig')
        summary_output['distance'] = summary_output.groupby(['Marker','decorate_Cell_type'])['distance'].transform('min')
        summary_output['subtree_number'] = summary_output.groupby(['Marker','decorate_Cell_type'])['subtree_number'].transform('min')
        summary_output = summary_output.drop_duplicates()
        standard_cellmarker_output = summary_output[summary_output["subtree_number"]<=2]
        if standard_cellmarker_output.shape[0]>0:
            standard_cellmarker_output = standard_cellmarker_output[['PMID',keyword,"decorate_Cell_type" ,"Marker","context_background" ,"Species","Tissue","Disease",'Title','Date','Journal']]
            standard_cellmarker_output.rename(columns={"decorate_Cell_type":"Cell_type","context_background":"Context_background"},inplace=True)
            standard_cellmarker_output.drop_duplicates().to_csv(os.path.dirname(output_path)+"/cellmarker_result.csv",index=False,na_rep='',encoding='utf_8_sig')


def extract_cell_marker(each_Segment_line,gene_name_collect="[adsaddsada]",cell_name_collect="[adsaddsada]",abbreviations_dict=dict()):
    p = inflect.engine()
    #doc = nlp_model_scibert(re.sub("-"," ",each_Segment_line))
    for each_gene_name_collect in gene_name_collect:
        each_Segment_line =re.sub(each_gene_name_collect+"\s{0,1}\^{0,1}\-|\−", "",each_Segment_line  )
    doc = nlp_model_scibert(each_Segment_line)
    tt_dataframe=pd.DataFrame()
    tt_dataframe = pd.concat([tt_dataframe,punct_segment_extract(doc)])  
    punct_index= tt_dataframe.shape[0]
    return_dict=subtree_all(doc)
    tt_dataframe = pd.concat([tt_dataframe,return_dict["output_dataframe"]]) 
    tt_dataframe['extract_source'] = ""
    tt_dataframe['extract_source'][0:punct_index]="punct_src"
    tt_dataframe['extract_source'][punct_index:]="nlp_src"
    cell_count_in_segment = 0
    gene_count_in_segment = 0
    cell_name_in_segment = ""
    gene_name_in_segment = ""
    for each_cell_name_collect in cell_name_collect:
        if re.search( re.escape(each_cell_name_collect) , abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_Segment_line.strip() )), abbreviations_dict=  abbreviations_dict ) ) :
            cell_count_in_segment+=1
            cell_name_in_segment=each_cell_name_collect
    if cell_count_in_segment==1:
        for each_gene_name_collect in gene_name_collect:
            if re.search( re.escape(each_gene_name_collect) , abbreviation2fullname(each_line=re.sub("-"," ",re.sub("/"," or ",each_Segment_line.strip() )), abbreviations_dict=  abbreviations_dict ) ) :
                gene_count_in_segment +=1
                gene_name_in_segment = each_gene_name_collect
    if gene_count_in_segment==1:           
            tt_dataframe = pd.concat([tt_dataframe,pd.DataFrame({"left_col":[gene_name_in_segment],"right_col":[cell_name_in_segment],"distance":1,"same_ancestors":1,"subtree_number":0,"extract_source":"segment_direct"})]) 
    if tt_dataframe.shape[0]>0:
        tt_dataframe = tt_dataframe.assign(left_col=tt_dataframe['left_col'].str.split('\\/')).explode('left_col')
        tt_dataframe =tt_dataframe.assign(right_col=tt_dataframe['right_col'].str.split('\\/')).explode('right_col')
        temp_dataframe1 = tt_dataframe[tt_dataframe['left_col'].str.contains("|".join(gene_name_collect)) ]
        temp_dataframe1.columns = ['Marker', 'Org_Cell_type',"distance","same_ancestors","subtree_number","extract_source"]
        temp_dataframe2 = tt_dataframe[tt_dataframe['right_col'].str.contains("|".join(gene_name_collect)) ]
        temp_dataframe2.columns = ['Org_Cell_type', 'Marker',"distance","same_ancestors","subtree_number","extract_source"]
        temp_oututdata = pd.concat([temp_dataframe1,temp_dataframe2])
        temp_oututdata['decorate_Cell_type']="pass"
        for x in range(0,temp_oututdata.shape[0]):
            Org_Cell_type_fullname =  abbreviation2fullname(temp_oututdata[['Org_Cell_type']].iloc[x,0],abbreviations_dict)
            if  Org_Cell_type_fullname.strip()!="" and p.singular_noun(Org_Cell_type_fullname) :
                Org_Cell_type_fullname=p.singular_noun(Org_Cell_type_fullname)    
            if Org_Cell_type_fullname in cell_name_collect:
                find_index = find_indices(cell_name_collect,Org_Cell_type_fullname)
                if len( find_index)==1:
                    temp_oututdata.iat[x, temp_oututdata.columns.get_loc("decorate_Cell_type")] =  cell_name_collect[find_index[0]]
                elif len( find_index)>1:
                    temp_oututdata.iat[x, temp_oututdata.columns.get_loc("decorate_Cell_type")] =  Org_Cell_type_fullname
            for y in cell_name_collect:
                if y in Org_Cell_type_fullname:
                    temp_oututdata.iat[x, temp_oututdata.columns.get_loc("decorate_Cell_type")] =y
        temp_oututdata = temp_oututdata[temp_oututdata['decorate_Cell_type']!="pass"]   
        if "punct_src" in temp_oututdata.extract_source.tolist():
            temp_oututdata = temp_oututdata[ temp_oututdata['extract_source'].isin(["punct_src","segment_direct"])]
        temp_oututdata = temp_oututdata.drop_duplicates()
        output_dataframe_extarct_gene=pd.DataFrame()
        for row_index in range(0,temp_oututdata.shape[0]):
            for each_gene_name_collect in gene_name_collect:
                if re.search(re.escape(each_gene_name_collect),temp_oututdata[['Marker']].iloc[row_index,0])  and not re.search(re.escape(each_gene_name_collect)+"\s{0,1}\^{0,1}\-|\−",temp_oututdata[['Marker']].iloc[row_index,0]):
                    temp_append_dataframe = temp_oututdata.iloc[[row_index]]
                    temp_append_dataframe['Marker'].iloc[0] = each_gene_name_collect
                    output_dataframe_extarct_gene = pd.concat([temp_append_dataframe,output_dataframe_extarct_gene])
        output_dataframe_extarct_gene = output_dataframe_extarct_gene.drop_duplicates()
    else:
        output_dataframe_extarct_gene=tt_dataframe
    puzzle_dataframe = return_dict["puzzle_datafrme"]
    return({"temp_oututdata":output_dataframe_extarct_gene,"puzzle_dataframe":puzzle_dataframe})


def count_distance(text, w1, w2):
   index1 = None
   index2 = None
   distance = 2000000
   for idx, word in enumerate(text.lower().split(" ")):
      if word == w1.lower():
         if index2 is not None:
            distance = min(distance, abs(idx - index2) - 1)
         index1 = idx
      if word == w2.lower():
         if index1 is not None:
            distance = min(distance, abs(idx - index1) - 1)
         index2 = idx
   if index1 is not None and index2 is not None:
      return distance
   return 100000


def nearest_values_twolist(list1,list2):
    min_val = 1000000
    if len(list1)==0 or len(list2)==0:
        return(min_val,0)
    else:
        r1 = list1[0]
        r2 = list2[0]
        for row1 in list1:
            for row2 in list2:
                t = abs(row1 - row2)
                if t<min_val:
                    min_val = t
                    r1 = row1
                    r2 = row2
        return(r1,r2)

def closest_distance_words(text,w1,w2):
    try:
        ind1 = [w.start(0) for w in re.finditer(w1.lower()+' ', text.lower())] +[w.start(0) for w in re.finditer(' '+w1.lower()+' ', text.lower())]+[w.start(0) for w in re.finditer(' '+w1.lower()+'\.', text.lower())]
        ind2 = [w.start(0) for w in re.finditer(w2.lower()+' ', text.lower())] +[w.start(0) for w in re.finditer(' '+w2.lower()+' ', text.lower())]+[w.start(0) for w in re.finditer(' '+w2.lower()+'\.', text.lower())]
        i1,i2 = nearest_values_twolist(ind1,ind2)
        return(abs(i2-i1))
    except :
        return(1000000)


########### 额外补充数据预处理 & 储存
##### 1.  基因实体库
####grep -Eo 'gene_name \S+;' /mnt/icfs/bioinfodatabase/singlecell/10X/refdata-gex-GRCh38-2020-A/genes/genes.gtf | sed 's#gene_name "##' | sed 's#";##'| sort | uniq -d  > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/GTF_hg38_gene.txt
####grep -Eo 'gene_name \S+;' /mnt/icfs/bioinfodatabase/singlecell/10X/refdata-cellranger-mm10-3.0.0/genes/genes.gtf | sed 's#gene_name "##' | sed 's#";##'| sort | uniq -d  > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/GTF_mmu10_gene.txt
####awk -F"\t" '{print $7}' /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/MRK_List1.rpt > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/MGI_mmu_gene.txt
#### HGNC : https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit
#GTF_hg38_gene = pd.read_table("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/GTF_hg38_gene.txt",header=None,encoding='utf-8')
#HGNC_human_gene = pd.read_table("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/HGNC_human_gene.txt",header=0,encoding='utf-8')
#GTF_mmu10_gene = pd.read_table("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/GTF_mmu10_gene.txt",header=None,encoding='utf-8')
#MGI_mmu10_gene = pd.read_table("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/MGI_mmu_gene.txt",encoding='utf-8')
#不用了，会引入额外噪音 ： gene_collect = list(GTF_hg38_gene.iloc[:,0]) + list(GTF_mmu10_gene.iloc[:,0]) + list(MGI_mmu10_gene.iloc[:,0]) + list(HGNC_human_gene.iloc[:,0]) +  sum([str(x).split(",") for x in HGNC_human_gene.iloc[:,1]],[]) + sum([str(x).split(",") for x in HGNC_human_gene.iloc[:,2]],[])
#不用了，会引入额外噪音 ：  gene_collect = list(GTF_hg38_gene.iloc[:,0]) + list(GTF_mmu10_gene.iloc[:,0]) + list(HGNC_human_gene.iloc[:,0])  
#  gene_collect = list(GTF_hg38_gene.iloc[:,0]) + list(GTF_mmu10_gene.iloc[:,0]) 
### 不用了会引入额外噪音  part_gene = sum([str(x).split(",") for x in HGNC_human_gene.iloc[:,1]],[]) + sum([str(x).split(",") for x in HGNC_human_gene.iloc[:,2]],[])
### 不用了 会引入额外噪音  regex = re.compile("[A-Za-z].*[0-9].*")
### 不用了会引入额外噪音 part_gene =[i for i in part_gene if  regex.match(str(i))]
#### 不用了会引入额外噪音  gene_collect =set([str(x).strip() for x in gene_collect+part_gene if len(str(x).strip())>1 ] ) - set('nan') 
#gene_collect =set([str(x).strip() for x in gene_collect if len(str(x).strip())>1 ] ) - set('nan') 
#regex = re.compile("S?[0-9]{1,2}[A-Za-z]{0,1}$")
#gene_collect = [i for i in gene_collect if not regex.match(str(i))]
#regex = re.compile("[A-Za-z][0-9]{1,2}$")
#gene_collect = [i for i in gene_collect if not regex.match(str(i))]
#gene_collect = pd.DataFrame(gene_collect,columns=['gene_collect'])
#gene_collect['label']="original"  ### original :从数据库得到的gene list  . append : 根据语言库识别出的新gene
#gene_collect.to_csv(gene_collect_path,index=False,encoding='utf_8_sig')
##### 2.  细胞实体库
# grep "^CL:" /mnt/icfs/work/singlecellproject/personal/binbinfang/Cellmarker_NLP/marker_database/ontoProc/cellterm.txt > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cell_ontology_cellname.txt
# cell_ontology = pd.read_table("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cell_ontology_cellname.txt",encoding='utf-8',header=None)
##### 3.  组织实体库
# grep "^UBERON" /mnt/icfs/work/singlecellproject/personal/binbinfang/Cellmarker_NLP/marker_database/ontoProc/cellterm.txt > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cell_ontology_tissue.txt
##tissue_collect=pd.read_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/all_tissue.csv")
##tissue_pattern = [{"label": "tissue_from_ather", "pattern": x.lower()} for x in set(tissue_collect.iloc[:,0]) ]


##############  extra data  暂时不用的额外数据
#################### 组织库 ：不用了
#cut  -f4 /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/human_tissue_integrated_full.tsv | grep -v ","  > /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/human_tissue.txt
#cut  -f4 /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/mouse_tissue_integrated_full.tsv | grep -v ","  > /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/mouse_tissue.txt
#cut  -f4 /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/pig_tissue_integrated_full.tsv | grep -v ","  > /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/pig_tissue.txt
#cut  -f4 /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/rat_tissue_integrated_full.tsv | grep -v ","  > /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/rat_tissue.txt
#human_tissue_collect = pd.read_table("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/human_tissue.txt",header=None)
#mouse_tissue_collect = pd.read_table("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/mouse_tissue.txt",header=None)
#pig_tissue_collect = pd.read_table("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/pig_tissue.txt",header=None)
#rat_tissue_collect = pd.read_table("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/rat_tissue.txt",header=None)
#tissue_collect = pd.concat([human_tissue_collect, mouse_tissue_collect, pig_tissue_collect, rat_tissue_collect])
#tissue_collect.columns=["tissue_collect"]
#tissue_collect=tissue_collect.drop_duplicates(subset=['tissue_collect'])
#tissue_collect['tissue_collect'] = tissue_collect.tissue_collect.str.upper()
#tissue_collect.to_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/all_tissue.csv",index=False)
#tissue_collect=pd.read_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/all_tissue.csv")
#tissue_database = set(tissue_collect["tissue_collect"])
#################### 细胞库 ：不用了
#cell_ontology = pd.read_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/Cell Ontology terms list - resource_browser.csv",encoding='utf-8')
#cell_source_1 = pd.read_table("/mnt/icfs/work/singlecellproject/personal/binbinfang/Cellmarker_NLP/marker_database/ASs.txt",header=None,encoding='utf-8')
#cell_source_2 = pd.read_table("/mnt/icfs/work/singlecellproject/personal/binbinfang/Cellmarker_NLP/marker_database/CellTypes.txt",header=None,encoding = "ISO-8859-1")
#cell_contact = pd.concat([cell_source_1,cell_source_2,cell_ontology['Name'][0:525] ])
#cell_contact.columns=["cell_collect"]
#cell_contact=cell_contact.drop_duplicates(subset=['cell_collect'])
#cell_contact['label']="original"
#cell_contact = cell_contact[cell_contact['cell_collect']!="cell"]
#cell_contact = cell_contact[cell_contact['cell_collect']!="cells"]
#cell_contact.to_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cell_collect.csv",index=False,encoding='utf_8_sig')


################  3.  txt 2 segment
if "txt2segment" in do_part:
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    for file in tqdm(paper_id_file.iloc[:,0]):
        try:
            file =str(int(file))
            file_path1 = work_path+"/paper_process_result/"+file+"/"+file+".tsv"
            pmc_file_path = work_path+"/paper_process_result/"+file+"/PMC.tsv"
            output_file_path = work_path+"/paper_process_result/"+file+"/segment.txt"
            abbreviations_dict_path = work_path+"/paper_process_result/"+file+"/abbreviations_dict.pickle"
            if os.path.exists(pmc_file_path):
                txt2segment(file_path=pmc_file_path, output_file_path=output_file_path,abbreviations_dict_path=abbreviations_dict_path)
            else:
                txt2segment(file_path=file_path1, output_file_path=output_file_path,abbreviations_dict_path=abbreviations_dict_path)
        except :
            print(file+": txt2segment error")
    



################  4.  extract potentially  valuable  segment
if "extract_potentiall_valuable_segment" in do_part:
    #### import gene
    ####### main function
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    for each_paper_id in tqdm(paper_id_file.iloc[:,0]):
        try:
            paper_id = str(int(each_paper_id))
            open_file_path = work_path+"/paper_process_result/"+paper_id+"/segment.txt"
            output_file_path = work_path+"/paper_process_result/"+paper_id+"/"+"potentailly_valualbe_evidence.txt"
            output_fig_segment_file_path = work_path+"/paper_process_result/"+paper_id+"/"+"potentailly_valualbe_evidence_fig.txt"
            abbr_dict_path = work_path+"/paper_process_result/"+paper_id+"/abbreviations_dict.pickle"
            extract_potentiall_valuable_segment(open_file_path=open_file_path,output_file_path=output_file_path,output_fig_segment_file_path=output_fig_segment_file_path,abbr_dict_path=abbr_dict_path,paper_id=paper_id)
        except :
            print(each_paper_id+": extract_potentiall_valuable_segment error")
########### 5.  segment text clean   
if "add_segment_text_clean" in do_part:
    #### preprocess
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    punctuation_word = copy.deepcopy(string.punctuation)
    punctuation_word  = re.sub(r'\+|\-|\(|\)',"",punctuation_word)
    STOP_WORDS = STOP_WORDS - set(['on',"of","as","in","has","was","is","for","by","have","such"])
    #### main function
    for paper_id in tqdm(paper_id_file.iloc[:,0]):
        try:
            open_file_path = work_path+"/paper_process_result/"+str(int(paper_id))+"/potentailly_valualbe_evidence.txt"
            abbr_dict_path =work_path+"/paper_process_result/"+str(int(paper_id))+"/abbreviations_dict.pickle"
            add_segment_text_clean(open_file_path=open_file_path, abbr_dict_path=abbr_dict_path)
        except :
            print(paper_id+": add_segment_text_clean error")

########## 6. predict result 
if "add_predict_result" in do_part:
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    nlp_model= spacy.load(nlp_model_path )
    for paper_id in tqdm(paper_id_file.iloc[:,0]):
        try:
            input_filepath = work_path+"/paper_process_result/"+str(int(paper_id))+"/potentailly_valualbe_evidence_clean.txt"
            output_filepath = work_path+"/paper_process_result/"+str(int(paper_id))+"/potentailly_valualbe_evidence_AddPredict.txt"
            add_predict_result(input_filepath=input_filepath,nlp_model =nlp_model,output_filepath=output_filepath)
        except :
            print(paper_id+": potentailly_valualbe_evidence_clean error")

###### 7. extract cell-gene result
if "cell_marker_extract_main" in do_part:
    paper_id_file = pd.read_csv(paper_id_path,header=None)
    if paper_id_start!=None and paper_id_end!=None:
        paper_id_file =paper_id_file.iloc[int(paper_id_start) : int(paper_id_end),]
    #### main function
    for paper_id in tqdm(paper_id_file.iloc[:,0]):
        try:
            file_path=work_path+"/paper_process_result/"+str(int(paper_id))+"/potentailly_valualbe_evidence_AddPredict.txt"
            output_path=work_path+"/paper_process_result/"+str(int(paper_id))+"/summary.csv"
            abbreviations_dict_path_root = work_path+"/paper_process_result/"
            cell_marker_extract_main(file_path=file_path,output_path=output_path, abbreviations_dict_path_root=abbreviations_dict_path_root ,keyword = "Evidence")
        except :
            print(paper_id+": cell_marker_extract_main error")

'''



## train model 
#### step1 ： 构建训练集
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/train_S2_make_spacy_dataset.py --train_dataset_path /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/text_model/group_train_dataset.csv
#### step2 ： 训练模型
/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python -m spacy train /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/text_model/config.cfg  --verbose --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh/text_model/train_model_all_data/





########  以下 旧脚本，废弃不用
#### part1 train model 150篇   0 vs 1 
# head -n 1 /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/manual_label_dataset/dataset_1/28576768_train_dataset.csv > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv
# cat /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/manual_label_dataset/*/*_train_dataset.csv > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.tmp
# grep -v "^PMID" /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.tmp >> /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv

model_train_info_outpath ="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv"
train_info = pd.read_csv(model_train_info_outpath,encoding='ISO-8859-1')

train_info = train_info[train_info.Label.isin([0,1,2,3])]
train_info = train_info[~train_info.Segment_clean.str.startswith("Paper_")]

train_info = train_info[train_info.Label.isin([0,1])]
train_info['Label'] =train_info['Label'].astype(int).astype(str)

train = train_info.sample(frac = 0.8, random_state = 25)
test = train_info.drop(train.index)

train['tuples'] = train.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
train =train['tuples'].tolist()
test['tuples'] = test.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
test =test['tuples'].tolist()

train_docs = document(train)
doc_bin = DocBin(docs = train_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train.spacy")

test_docs = document(test)
doc_bin = DocBin(docs = test_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/test.spacy")

os.system("python3 -m spacy train /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/config.cfg --verbose --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_model_150_01/")

#### part1 train model 150篇   0 vs 1 和 2  暂时不做

model_train_info_outpath ="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv"
train_info = pd.read_csv(model_train_info_outpath,encoding='ISO-8859-1')

train_info = train_info[train_info.Label.isin([0,1,2])]

train = train_info.sample(frac = 0.8, random_state = 25)
test = train_info.drop(train.index)

train['tuples'] = train.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
train =train['tuples'].tolist()
test['tuples'] = test.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
test =test['tuples'].tolist()

train_docs = document(train)
doc_bin = DocBin(docs = train_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train.spacy")

test_docs = document(test)
doc_bin = DocBin(docs = test_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/test.spacy")

os.system("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python -m spacy train /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/config.cfg --verbose --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_model_150_01/")


########## part2 train test data
test_data =pd.read_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cellmarker_pdf_part/total_potentailly_valualbe_evidence_clean_add_label.csv")
test_data = test_data[test_data.Label.isin([True,False])]

nlp_model_150_01 = spacy.load("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_model_150_01/model-best")

def model_result_true(x,keyword):
    temp=nlp_model_150_01(str(x[keyword])).cats
    return(temp['True'])

test_data.insert(loc=test_data.shape[1], column='Predict_result', value=test_data.apply( model_result_true,keyword="Segment_clean",axis=1) )

test_data.to_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/cellmarker_pdf_part/total_potentailly_valualbe_evidence_clean_add_label_add_predict.csv",index=False,encoding='utf_8_sig')



######### 35篇预测150篇
nlp_model_35 = spacy.load("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated/model-best")

test_data =pd.read_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv",encoding='ISO-8859-1')
test_data = test_data[test_data.Label.isin([0,1,2,3])]
test_data = test_data[~test_data.Segment_clean.str.startswith("Paper_")]

test_data = test_data[test_data.Label.isin([0,1])]
test_data['Label'] =test_data['Label'].astype(int).astype(str)

def model_result_true(x,keyword):
    temp=nlp_model_35(str(x[keyword])).cats
    return(temp['True'])

test_data.insert(loc=test_data.shape[1], column='Predict_result', value=test_data.apply( model_result_true,keyword="Segment_clean",axis=1) )

test_data.to_csv("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info_add_predict.csv",index=False,encoding='utf_8_sig')





#### part1 train model 
# head -n 1 /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/manual_label_dataset/dataset_1/28576768_train_dataset.csv > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv
# cat /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/manual_label_dataset/*/*_train_dataset.csv > /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.tmp
# grep -v "^PMID" /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.tmp >> /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv

model_train_info_outpath ="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_info.csv"
train_info = pd.read_csv(model_train_info_outpath,encoding='ISO-8859-1')

train = train_info.sample(frac = 0.8, random_state = 25)
test = train_info.drop(train.index)

train['tuples'] = train.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
train =train['tuples'].tolist()
test['tuples'] = test.apply(lambda row: (str(row['Segment_clean']),str(row['Label'])), axis=1)
test =test['tuples'].tolist()

train_docs = document(train)
doc_bin = DocBin(docs = train_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train.spacy")

test_docs = document(test)
doc_bin = DocBin(docs = test_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/test.spacy")

os.system("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python -m spacy train /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/config.cfg --verbose --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/train_dataset/train_model_150/")





#############
nlp_model = spacy.load("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated/model-best")

predict_textcat(paper_id_path=paper_id_path,work_path=work_path,input_file_name="potentailly_valualbe_evidence.txt")

##########  supplment : train model 
model_train_info_outpath ="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/model_train_info.csv"
train_info = pd.read_csv(model_train_info_outpath,header=0,encoding='utf-8')


#temp = list(map(preprocessing, train_info["Evidence"]))
gene_collect = pd.read_csv(gene_collect_path)
use_gene_collect = set([str(x) for x in list(gene_collect['gene_collect'])])
temp =list()
for i in tqdm(range(0,train_info.shape[0])):
    temp.append(preprocessing(text=train_info["Evidence"][i],use_gene_collect=use_gene_collect,abbreviations_dict_path=work_path+"/paper_process_result/"+str(train_info["PMID"][i])+"/abbreviations_dict.pickle"))

train_info["text_clean"] = temp


train = train_info.sample(frac = 0.8, random_state = 25)
test = train_info.drop(train.index)

train['tuples'] = train.apply(lambda row: (str(row['text_clean']),str(row['Label'])), axis=1)
train =train['tuples'].tolist()
test['tuples'] = test.apply(lambda row: (str(row['text_clean']),str(row['Label'])), axis=1)
test =test['tuples'].tolist()


def document(data):
# 创建名为“text”的空列表
    text = []
    for doc, label in nlp_model_scibert.pipe(data, as_tuples = True):
        if (label=='True'):
            doc.cats['True'] = 1
            doc.cats['False'] = 0
        elif (label=='False'):
            doc.cats['True'] = 0
            doc.cats['False'] = 1
# 将文档添加到列表“text”中
        text.append(doc)
    return(text)

def document(data):
# 创建名为“text”的空列表
    text = []
    for each_data in data:
        try:
            doc = nlp_model_scibert(each_data[0])
            if (each_data[1]=='1'):
                doc.cats['True'] = 1
                doc.cats['False'] = 0
            elif (each_data[1]=='0'):
                doc.cats['True'] = 0
                doc.cats['False'] = 1
            text.append(doc)
        except Exception as e :
            print([e])
    return(text)

        




train_docs = document(train)
doc_bin = DocBin(docs = train_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP//train.spacy")

test_docs = document(test)
doc_bin = DocBin(docs = test_docs)
doc_bin.to_disk("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP//test.spacy")

os.system("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python -m spacy train /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/config.cfg --verbose --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated/")
os.system("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python -m spacy evaluate /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated//model-best/ --output /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/metrics.json /mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/test.spacy")


nlp_model = spacy.load("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated/model-best")

train_info.insert(loc=train_info.shape[1], column='Predict_result', value=train_info.apply( model_result,keyword="text_clean",axis=1) )
train_info.insert(loc=train_info.shape[1], column='Predict_result_true', value=train_info.apply( model_result_true,keyword="text_clean",axis=1) )
train_info.insert(loc=train_info.shape[1], column='Predict_result_false', value=train_info.apply( model_result_false,keyword="text_clean",axis=1) )
train_info.to_csv(model_train_info_outpath,index=False,encoding='utf_8_sig')
'''

'''
##### make train_label 
nlp_model = spacy.load("/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/output_updated/model-best")

paper_id_path= "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/manual_pubmedid.csv"
work_path = "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part"
paper_id_file = pd.read_csv(paper_id_path,header=None)

count=0
for paper_id in tqdm(paper_id_file.iloc[:,0]):
    paper_id_use = str(int(paper_id))
    if os.path.exists(work_path+"/paper_process_result/"+paper_id_use+ "/potentailly_valualbe_evidence_clean.txt"):
        train_data = pd.read_table(work_path+"/paper_process_result/"+paper_id_use+ "/potentailly_valualbe_evidence_clean.txt", delimiter="\t", quoting=csv.QUOTE_NONE)
        train_data =train_data[['PMID','context_background','Species', 'Disease',"Tissue_type","Judge_reliability(SDT)",'Weight', 'Gene_name','Cell_name', "Segment_clean",'Segment']]
        train_data.insert(loc=train_data.shape[1], column='Predict_prob', value=train_data.apply( model_result_true,keyword="Segment_clean",axis=1) )
        train_data = train_data.sort_values(by="Predict_prob",ascending=False)
        train_data['Label']=0
        train_data.to_csv(work_path+"/paper_process_result/"+paper_id_use+ "/"+paper_id_use+"_train_dataset.csv",index=False,encoding='utf_8_sig')
'''


