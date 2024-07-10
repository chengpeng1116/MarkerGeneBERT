import sys
import argparse
import os
import re
from tqdm import tqdm

### /mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python

main_sh_path = os.getcwd()
root_sh_dir = sys.path[0] 
# root_sh_dir = "/mnt/rorke/personal/pengcheng/project_collect/Cellmarker_NLP/pipeline_test/pipeline_sh"
def make_dir(temp_dir):
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


parser = argparse.ArgumentParser(description=' ')
parser.add_argument('--cfg', type = str, default = None)


args = parser.parse_args()

cfg = args.cfg 
config_path = main_sh_path + "/" + args.cfg


config = {}
for config_row in open(config_path):
    if not config_row.startswith("#"):
        if re.match("\s*\S*\s*=+?.*", config_row) != None:
            break_index = config_row.find("=")      
            key = config_row[:break_index].lstrip().strip()
            value = config_row[break_index + 1:]
            if re.match(".*\s+#.*", value) != None:    
                value = re.sub("\s+#.*", "", value)
            value = value.lstrip().strip()
            config[key] = value   


Work_path =  main_sh_path+ "/" 
DB_path =  config['DB_path']
if not  "update_pmid" in config.keys():
    update_pmid =  Work_path+"/update_pmid.csv"
else:
    update_pmid =  config['update_pmid']

Analysis_path = Work_path+"/paper_process_result"
make_dir(Analysis_path)

DB_path_file  =  os.listdir(DB_path)
check_status = 0
for each_DB_path_file in tqdm(DB_path_file):
    each_DB_path_file_path = DB_path +"/" + each_DB_path_file
    if os.path.isdir(each_DB_path_file_path) and each_DB_path_file!="temp_PDF":
        analysis_path_file = Analysis_path + "/" + each_DB_path_file
        make_dir(analysis_path_file)
        if  os.path.exists(each_DB_path_file_path+"/"+each_DB_path_file+".pdf") and not os.path.exists(analysis_path_file+"/"+each_DB_path_file+".pdf") :
            #if not os.path.exists(analysis_path_file+"/"+each_DB_path_file+".pdf"):
            #    os.system("ln -s "+each_DB_path_file_path+"/"+each_DB_path_file+".pdf "+ analysis_path_file )
            #elif os.stat(each_DB_path_file_path+"/"+each_DB_path_file+".pdf").st_size != os.stat(analysis_path_file+"/"+each_DB_path_file+".pdf").st_size:
            os.system("ln -sf "+each_DB_path_file_path+"/"+each_DB_path_file+".pdf "+ analysis_path_file )
            check_status =1
        if  os.path.exists(each_DB_path_file_path+"/abstract.tsv") and not os.path.exists(analysis_path_file+"/abstract.tsv") :
            #if not os.path.exists(analysis_path_file+"/abstract.tsv"):
                #os.system("ln -s "+each_DB_path_file_path+"/abstract.tsv "+ analysis_path_file )
            #elif os.stat(each_DB_path_file_path+"/abstract.tsv").st_size != os.stat(analysis_path_file+"/abstract.tsv").st_size:
            os.system("ln -sf "+each_DB_path_file_path+"/abstract.tsv "+ analysis_path_file )     
            check_status =1
        if  os.path.exists(each_DB_path_file_path+"/PMC.tsv") and not os.path.exists(analysis_path_file+"/PMC.tsv") :
            #if not os.path.exists(analysis_path_file+"/PMC.tsv"):
                #os.system("ln -s "+each_DB_path_file_path+"/PMC.tsv "+ analysis_path_file )
            #elif os.stat(each_DB_path_file_path+"/PMC.tsv").st_size != os.stat(analysis_path_file+"/PMC.tsv").st_size:
            os.system("ln -sf "+each_DB_path_file_path+"/PMC.tsv "+ analysis_path_file )       
            check_status =1

step_sh = open(Work_path + "/step1_check_PMID_status.sh", "w")
print("# 用于检测分析路径下各文献分析状态，并将分析状态结果保存至DB_PMID_status.csv" , file= step_sh)
order = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+root_sh_dir+"/check_PMID_status.py  --check_path " +Analysis_path +" --output_path " +Work_path+"/DB_PMID_status.csv"
print(order, file= step_sh)
step_sh.close()

if check_status==1 or not os.path.exists(Work_path+"/DB_PMID_status.csv") :
    os.system("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+root_sh_dir+"/check_PMID_status.py  --check_path " +Analysis_path +" --output_path " +Work_path+"/DB_PMID_status.csv")


step_sh = open(Work_path + "/step2_update_DB.sh", "w")
print("#根据关键词从PUBMED得到PMID列表，保存于update_pmid.csv ", file= step_sh)
order = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript "+root_sh_dir+"/batch_download_PMID.R " + ' --keyword '+"'" + '"Animals"[MeSH Terms] AND ("Single-Cell Analysis"[MeSH Terms] OR ("single-cell" AND "expression"))' +"'" +" --mindate 2020  --maxdate 2020  --output_file_path "+ update_pmid+" &&\\\n"
print(order, file= step_sh)
print("#将update_pmid.csv文件中的PMID与DB_PMID_status.csv进行比较，输出本数据库未收录的PMID列表，保存于文件update_PMID_checked.csv", file= step_sh)
order = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+root_sh_dir+"/update_DB.py  --update_pmid " +update_pmid +" --DB_PMID_status " +Work_path+"/DB_PMID_status.csv --output_path "+ Work_path+"/update_PMID_checked.csv"
print(order, file= step_sh)
step_sh.close()


step_sh = open(Work_path + "/step3_download_PMID.sh", "w")
print("#使用R脚本从tidyPMC下载XML格式的文献。下载成功则更新至数据库目录下，下载失败则将PMID保存至update_PMID_checked_PMC_download_fail.csv", file= step_sh)
order = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript "+root_sh_dir+"/extract_paper_info.R " + Work_path+"/update_PMID_checked.csv" +" " + DB_path+ " " + Work_path+"/update_PMID_checked_PMC_download_fail.csv &&\\\n"
print(order, file= step_sh)
temp_pdf_path = DB_path+"/temp_pdf/"
make_dir(temp_pdf_path)
print("#使用scihub，将update_PMID_checked_PMC_download_fail.csv中的文献利用scihub下载PDF格式的文献", file= step_sh)
order = "##下载量小时用 ： /mnt/icfs/personal/xiaolingzhang/.local/bin/scihub -ns -ow N -s "+ Work_path+"/update_PMID_checked_PMC_download_fail.csv" +"  -O "+ temp_pdf_path +"  -u https://sci-hub.ren"
print(order, file= step_sh)

make_dir(Work_path+"/sub_step3_pdf_download/")
string_s = '''
import pandas as pd 
import re 
import os 
from tqdm import tqdm

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


tt = pd.read_csv("'''+ Work_path+'/update_PMID_checked_PMC_download_fail.csv",header=None) '+'''

tt.columns=["PMID"]
temp = list(set(tt['PMID'].tolist()))
temp = [str(x ) for x in temp]
split_order_list= list(split(range(len(temp)), 10))
temp_use = [ temp[x]  for x in split_order_list[count_s] ]
download_pdf =  os.listdir("'''+root_sh_dir+'''/DB_source/temp_pdf/")
download_pdf = [re.sub("\\.pdf","",x) for x in download_pdf ]
need_update = list(set(temp_use)-set(download_pdf))
need_update =  [int(x) for x in need_update]
fail_pmid = []
for each_need_update in tqdm(need_update):
    i=0
    fail_pmid.append(each_need_update)
    while i <=4:
        i+=1
        try:
            os.system("/mnt/icfs/personal/xiaolingzhang/.local/bin/scihub -ns -ow N -s "+ str(each_need_update) +"  -O '''+root_sh_dir+'''/DB_source/temp_pdf/  -u https://sci-hub.ren") 
            if os.path.exists(root_sh_dir+"/DB_source/temp_pdf/"+str(each_need_update)+".pdf"):
                i=10
                no_use = fail_pmid.pop()
        except:
            i+=1


fail_pmid = pd.DataFrame({"PMID":fail_pmid})
fail_pmid.to_csv("'''+Work_path+'''/sub_step3_pdf_download/pdf_fail_"+str(count_s)+".csv",index=False   )

'''
print("#cd  "+Work_path+"/sub_step3_pdf_download/",file= step_sh)
for i in range(0,10):
    output_file = open(Work_path+"/sub_step3_pdf_download/run_pdf_"+str(i)+".sh","w")
    print("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+Work_path+"/sub_step3_pdf_download/run_pdf_"+str(i)+".py",file=output_file)
    output_file.close()
    output_file1 = open(Work_path+"/sub_step3_pdf_download/run_pdf_"+str(i)+".py","w")
    print("count_s="+str(i),file=output_file1)
    print(string_s ,file=output_file1)
    output_file1.close()
    string_temp = "# qsub " + Work_path+"/sub_step3_pdf_download/run_pdf_"+str(i)+".sh  "
    print(string_temp,file= step_sh)

print("#将PDF格式的文献转存至数据库目录下", file= step_sh)
order = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/R4.1/bin/Rscript "+root_sh_dir+"/mv_pdf_download_fromTemp.R " + temp_pdf_path+ " " + DB_path
print(order, file= step_sh)
step_sh.close()


step_sh = open(Work_path + "/step4_pdf2txt_a04.sh", "w")
print("#将分析路径下的文献PDF解析为txt格式"  ,file=step_sh ) 
print("/mnt/icfs/personal/xiaolingzhang/miniconda3/envs/grobid/bin/python "+root_sh_dir+"/pdf2txt.py --Analysis_path " +Analysis_path ,file=step_sh ) 
step_sh.close()


step_sh = open(Work_path + "/step5_cellmarker_extract.sh", "w")
print("#将分析路径下的所有文献进行cell-marker的提取",file=step_sh )
print("cd "+Work_path + "/sub_step5_cellmarker_extract/",file=step_sh ) 
make_dir(Work_path + "/sub_step5_cellmarker_extract/")



if os.path.exists(Work_path+"/DB_PMID_status.csv"):
    count = len(open(Work_path+"/DB_PMID_status.csv",'r').readlines())
    if count < 50:
        sub_step_sh = open(Work_path + "/sub_step5_cellmarker_extract/s5_run.sh", "w")
        print("/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+root_sh_dir+"/std_all_part.py --work_path "+ Work_path + ' --paper_id_path '+ Work_path+"/DB_PMID_status.csv" ,file=sub_step_sh ) 
        sub_step_sh.close()
        print("qsub "+Work_path + "/sub_step5_cellmarker_extract/s5_run.sh",file=step_sh ) 
    else:
        split_order_list= list(split(range(count), 20))
        for i in range(0,len(split_order_list)):
            paper_id_start = split_order_list[i][0]
            paper_id_end = split_order_list[i][-1]+1
            order_string = "/mnt/icfs/work/singlecelldevelopment/miniconda3/envs/Cellmarker_NLP/bin/python "+root_sh_dir+"/std_all_part.py --work_path "+ Work_path + ' --paper_id_path '+ Work_path+"/DB_PMID_status.csv --paper_id_start " +str(paper_id_start) +" --paper_id_end "+str(paper_id_end)
            sub_step_sh = open(Work_path + "/sub_step5_cellmarker_extract/s5_run_"+str(i+1)+".sh", "w")
            print(order_string ,file=sub_step_sh)
            sub_step_sh.close()
            print("qsub "+Work_path + "/sub_step5_cellmarker_extract/s5_run_"+str(i+1)+".sh",file=step_sh ) 
else:
    print("DB_PMID_status.csv mot exist" ) 
step_sh.close()

step_sh = open(Work_path + "/step6_summary_cellmarker_from_DB.sh", "w")
print("#将分析路径下的所有文献进行cell-marker的提取",file=step_sh ) 
print("head -n 1 "+ Work_path+"/paper_process_result/28576768/cellmarker_result.csv > " +Work_path +'/cellmarker_result_all.csv', file=step_sh )
print("cat "+ Work_path+'/paper_process_result/*/cellmarker_result.csv | grep -v "PMID"  >> ' +Work_path +'/cellmarker_result_all.csv', file=step_sh )
step_sh.close()









