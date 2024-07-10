import argparse
import pandas  as pd 
import os 

import argparse
import pandas  as pd 
import os 

parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--update_pmid', type = str, default = None)
parser.add_argument('--DB_PMID_status', type = str, default = None)
parser.add_argument('--output_path', type = str, default = None)

args = parser.parse_args()
update_pmid=args.update_pmid
DB_PMID_status=args.DB_PMID_status
output_path=args.output_path

#DB_PMID_status="/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/manual_part/DB_PMID_status.csv"

DB_PMID_status_data  = pd.read_csv(DB_PMID_status)
DB_PMID_status_data_use = DB_PMID_status_data[DB_PMID_status_data["Download_status"].isin( ["Downloaded","Not_downloaded"]) ]


update_pmid_data = pd.read_csv(update_pmid)
temp = update_pmid_data.columns.tolist()
temp[0]="PMID"
update_pmid_data.columns=temp
update_pmid_data['PMID'] = update_pmid_data['PMID'].astype(str)


output_data = update_pmid_data[~update_pmid_data["PMID"].isin( DB_PMID_status_data_use["PMID"]) ]


output_data.to_csv(output_path,index=False)