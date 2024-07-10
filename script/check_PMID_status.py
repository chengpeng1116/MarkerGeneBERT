# /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python

import argparse
import pandas  as pd 
import os 
from tqdm import tqdm

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--check_path', type = str, default = None)
parser.add_argument('--output_path', type = str, default = None)

args = parser.parse_args()
check_path=args.check_path
output_path=args.output_path

#check_path="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/paper_process_result"
#output_path="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/DB_PMID_status.csv"

PMID_list= os.listdir(check_path)
PMID_list = [x for x in PMID_list if os.path.isdir(check_path+"/"+x)  ] +["temp"]

output_file = pd.DataFrame({"PMID":PMID_list,"Download_status":"","Download_type":"","Analysis_status":""})

for i in tqdm(range(0,len(PMID_list))):
    PMID = PMID_list[i]
    if os.path.exists(check_path+"/"+PMID+"/abstract.tsv"):
        output_file.loc[i, "Download_status"]="Downloaded"
        if os.path.exists(check_path+"/"+PMID+"/PMC.tsv"):
            output_file.loc[i, "Download_type"]="XML"
        elif os.path.exists(check_path+"/"+PMID+"/"+PMID+".pdf"):
            output_file.loc[i, "Download_type"]="PDF"
        else:
            output_file.loc[i, "Download_type"]="PDF_download_fail"
    else:
        output_file.loc[i, "Download_status"]="Not_downloaded"
    if os.path.exists(check_path+"/"+PMID+"/summary.csv"):
       output_file.loc[i, "Analysis_status"]="Analyzed"
    elif os.path.exists(check_path+"/"+PMID+"/scipdf_extract_error"):
        if not os.path.exists(check_path+"/"+PMID+"/"+PMID+".tsv") :
            output_file.loc[i, "Analysis_status"]="Not_Analyzed"
        else:
            output_file.loc[i, "Analysis_status"]="scipdf_error"
    else:
        if not output_file.loc[i, "Download_type"]=="Not_downloaded":
            output_file.loc[i, "Analysis_status"]="Not_Analyzed"


#### 把 未分析的的文献 PMID 均匀的分散开，用于提高后续并行更新时的速度
df1 = output_file[output_file.Analysis_status == "Analyzed"]
df2 = output_file[output_file.Analysis_status == "Not_Analyzed"]
df3 = output_file[output_file.Analysis_status == "scipdf_error"]

spacing = int(len(df1) / (len(df2) - 1))
df2ind = range(0, (len(df2) - 1) * (spacing + 1) + 1, spacing + 1)

df1ind = [range(i + 1, i + 1 + spacing) for i in df2ind]
df1ind = [i for l in df1ind for i in l][:len(df1)]

if df1ind!=[]:
    df1.index = df1ind
if df2ind!=[]:
    df2.index = df2ind


output_file = pd.concat([df1, df2,df3]).sort_index()


output_file.to_csv(output_path,index=False)
