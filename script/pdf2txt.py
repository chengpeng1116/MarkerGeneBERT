### /mnt/icfs/personal/xiaolingzhang/miniconda3/envs/grobid/bin/python
import argparse
import scipdf
import re
import pandas as pd 
import os 
from tqdm import tqdm


def use_segment_extract(pdf_path,output_path):
    if os.path.exists(pdf_path):
        if os.path.exists(output_path):
            print(output_path+" have exist ")
        else :
            article_dict = scipdf.parse_pdf_to_dict(pdf_path) 
            output_all= open(output_path, 'w', encoding='utf-8')
            #### select title and abstract , marked ,output 
            paper_title = article_dict['title']
            print("Paper_title: "+paper_title ,file=output_all) 
            paper_abstract = article_dict['abstract']
            paper_abstract = re.sub("\n"," ",paper_abstract)
            print("Paper_abstract: "+paper_abstract ,file=output_all) 
            ####  select segment from section and fig caption
            result_start=0
            section_dict= dict()
            background_marker=0
            disccuion_marker=0
            result_marker=0
            method_marker=0
            for section_number in range(0,len(article_dict['sections'])):
                if article_dict['sections'][section_number]['heading'].lower() in ['background',"introduction"] or article_dict['sections'][section_number]['heading'].lower().startswith("introduction") or article_dict['sections'][section_number]['text'].lower().startswith("introduction") or article_dict['sections'][section_number]['text'].lower().startswith("background"):
                    section_dict[section_number] = "background"
                    background_marker=1
                else:
                    background_marker=2
                    if article_dict['sections'][section_number]['heading'].lower() in ['result',"results"] or article_dict['sections'][section_number]['heading'].lower().startswith("result") or article_dict['sections'][section_number]['text'].lower().startswith("result"):
                        result_marker=1
                        section_dict[section_number] = "result"   
                    elif  article_dict['sections'][section_number]['heading'].lower() in ['discussion',"discussions"] or article_dict['sections'][section_number]['heading'].lower().startswith("discussion") or article_dict['sections'][section_number]['text'].lower().startswith("discussion"):
                        result_marker=2
                        disccuion_marker=1
                        section_dict[section_number] = "discussion"  
                    elif  article_dict['sections'][section_number]['heading'].lower().startswith("method") or article_dict['sections'][section_number]['heading'].lower().startswith("material") or article_dict['sections'][section_number]['text'].lower().startswith("material") or article_dict['sections'][section_number]['text'].lower().startswith("method"):
                        method_marker=1
                        section_dict[section_number] = "method"  
                    elif background_marker!=1 and  result_marker==0:
                        section_dict[section_number] = "other_background"
                    elif result_marker==1 and  disccuion_marker==0:
                        section_dict[section_number] = "result"
                    elif background_marker==2 and result_marker==2 and disccuion_marker==2 :
                        section_dict[section_number] = "method"
                    else:
                        section_dict[section_number] = "NULL"
            section_count=1
            method_count=1
            result_segment=""
            back_segment=""
            method_segment=""
            other_back_segment=""
            for section_number in range(0,len(article_dict['sections'])):                   
                if section_dict[section_number]=="background":
                    back_segment+=article_dict['sections'][section_number]['text']+" @@ "
                if section_dict[section_number]=="result":
                    result_segment+=article_dict['sections'][section_number]['text']+" "
                if section_dict[section_number]=="method":
                    if method_count<=4:
                        if article_dict['sections'][section_number]['text']!="":       
                            method_segment+=article_dict['sections'][section_number]['text']+" "
                            method_count+=1
                if section_dict[section_number]=="other_background":
                    other_back_segment+=article_dict['sections'][section_number]['text']+" @@ "
            back_segment = "@@ "+re.sub("\n"," @@ ",back_segment)
            other_back_segment = "@@ " + re.sub("\n"," @@ ",other_back_segment)
            print("Paper_intro: "+back_segment ,file=output_all) 
            #print("Paper_other_back: "+other_back_segment ,file=output_all) 
            if result_segment!="":
                section_1_segment =  result_segment.split("\n")[0]
                section_2_segment =  result_segment.split("\n")[1] 
                print("Paper_section1: "+section_1_segment ,file=output_all) 
                print("Paper_section2: "+section_2_segment ,file=output_all) 
            result_segment = re.sub("\n"," ",result_segment)
            print(result_segment ,file=output_all) 
            method_segment = re.sub("\n"," ",method_segment)   
            print("Paper_method: "+method_segment ,file=output_all) 
            for fig_number in range(0,len(article_dict['figures'])):
                fig_caption =  article_dict['figures'][fig_number]["figure_caption"]
                fig_caption = re.sub("\n"," ",fig_caption)
                print(fig_caption ,file=output_all) 
                    #### 提取图片 ，目前有报错，类似 需要载入 字体等信息，先搁置
                        # scipdf.parse_figures(pdf_path + filename + '.pdf', output_folder='') 
            output_all.close()
            print(pdf_path+" 2 txt done ")
    else:
        print("file do not exist : "+pdf_path)

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--Analysis_path',dest='Analysis_path', type = str, default = None)


args = parser.parse_args()
Analysis_path = args.Analysis_path 

Analysis_file  =  os.listdir(Analysis_path)
for each_Analysis_file in tqdm(Analysis_file):
    each_Analysis_file_path = Analysis_path +"/" + each_Analysis_file
    if os.path.isdir(each_Analysis_file_path) and os.path.exists(each_Analysis_file_path+"/"+each_Analysis_file+".pdf") and not os.path.exists(each_Analysis_file_path+"/PMC.tsv"):
        try:
            use_segment_extract(pdf_path=each_Analysis_file_path+"/"+each_Analysis_file+".pdf",output_path=each_Analysis_file_path+"/"+each_Analysis_file+".tsv" ) 
        except Exception as e :
            temp = open( each_Analysis_file_path+"/scipdf_extract_error","w")
            print(e,file=temp)
            temp.close()











