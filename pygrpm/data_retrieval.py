import os
import ast
import sys
import json
import nbib
import pandas as pd
import requests as rq
import matplotlib.pyplot as plt
from datetime import datetime
from bs4 import BeautifulSoup
from io import StringIO
from Bio import Entrez
import time



# Define Functions  ================
def transform_string(string):
    string = string.replace("\n",", ")
    string = string.replace("\'","\"")
    string = string.replace('\"\"', '\"')
    string = string.replace('p.\"','p.')
    string = string.replace('c.\"','c.')
    string = string.replace('g.\"','g.')
    string = string.replace('\">','>')
    string = string.replace('.C\"204','.C204') #<= if stucked, look for bugs like this into text
    return string

def query_build(pmid_list):
    query = "+OR+".join(pmid_list)
    return query

def build_pubmed_query(pmid_l, limit = 1300):
    query = []

    if len(pmid_l)<=limit:
        pmid_l01 = pmid_l
        query = [query_build(pmid_l01)]

    if limit<len(pmid_l)<=limit*2:
        j = len(pmid_l)//2
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02)]

    if limit*2<len(pmid_l)<=limit*3:
        j = len(pmid_l)//3
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:j*2]
        pmid_l03 = pmid_l[j*2:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02),
                 query_build(pmid_l03)]

    if limit*3<len(pmid_l)<=limit*4:
        j = len(pmid_l)//4
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:j*2]
        pmid_l03 = pmid_l[j*2:j*3]
        pmid_l04 = pmid_l[j*3:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02),
                 query_build(pmid_l03),
                 query_build(pmid_l04)]

    if limit*4<len(pmid_l)<=limit*5:
        j = len(pmid_l)//5
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:j*2]
        pmid_l03 = pmid_l[j*2:j*3]
        pmid_l04 = pmid_l[j*3:j*4]
        pmid_l05 = pmid_l[j*4:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02),
                 query_build(pmid_l03),
                 query_build(pmid_l04),
                 query_build(pmid_l05)]

    if limit*5<len(pmid_l)<=limit*6:
        j = len(pmid_l)//6
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:j*2]
        pmid_l03 = pmid_l[j*2:j*3]
        pmid_l04 = pmid_l[j*3:j*4]
        pmid_l05 = pmid_l[j*4:j*5]
        pmid_l06 = pmid_l[j*5:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02),
                 query_build(pmid_l03),
                 query_build(pmid_l04),
                 query_build(pmid_l05),
                 query_build(pmid_l06)]

    if limit*6<len(pmid_l)<=limit*7:
        j = len(pmid_l)//7
        pmid_l01 = pmid_l[:j]
        pmid_l02 = pmid_l[j:j*2]
        pmid_l03 = pmid_l[j*2:j*3]
        pmid_l04 = pmid_l[j*3:j*4]
        pmid_l05 = pmid_l[j*4:j*5]
        pmid_l06 = pmid_l[j*5:j*6]
        pmid_l07 = pmid_l[j*6:]
        query = [query_build(pmid_l01),
                 query_build(pmid_l02),
                 query_build(pmid_l03),
                 query_build(pmid_l04),
                 query_build(pmid_l05),
                 query_build(pmid_l06),
                 query_build(pmid_l07)]
    return query

def GetPubmedData(page, pmid):
    page = str(page)
    url = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + pmid + '&format=pubmed&size=200&page='+ page
    output = rq.get(url)
    html = output.text
    soup = BeautifulSoup(html, features="html.parser")
    for script in soup(["script", "style"]):
        script.extract()
    text = soup.get_text()
    postString = text.split("\n\n\n\n\n\n\n\n\n\n",2)[2]
    nbib01 = postString.replace('\n\n','')
    return nbib01

def save_checkpoint(complete_df, complete_nbibtable, db_path):
    complete_df = complete_df.reindex(columns=['gene','rsid', 'pmid', 'mesh', 'qualifier', 'major'])
    #complete_df.to_csv(db_path+'/grpm_table_output.csv')
    complete_df.to_parquet(db_path+'/grpm_dataset.parquet')

    grpm_nbib = complete_nbibtable.reindex(columns=['gene','pubmed_id', 'citation_owner', 'nlm_status', 'last_revision_date', 'electronic_issn', 'linking_issn', 'journal_volume', 'journal_issue', 'publication_date', 'title', 'abstract', 'authors', 'language', 'grants', 'publication_types', 'electronic_publication_date', 'place_of_publication', 'journal_abbreviated', 'journal', 'nlm_journal_id', 'descriptors', 'pmcid', 'keywords', 'conflict_of_interest', 'received_time', 'revised_time', 'accepted_time', 'pubmed_time', 'medline_time', 'entrez_time', 'pii', 'doi', 'publication_status', 'print_issn', 'pages'])
    #grpm_nbib.to_csv(db_path+'/complete_nbibtable.csv')
    grpm_nbib.to_parquet(db_path+'/grpm_nbib.parquet')
    print("saved checkpoint")

def plot_retrieval_data(df_count_sort, timestamp, db_path, gene):
    x1 = df_count_sort['mesh'].head(30)
    y1 = df_count_sort['pmid-unique'].head(30)
    plt.figure(figsize=(5, 8))
    plt.title('Scatter Plot: '+gene+' pmid-mesh (total)', loc='center',pad=10)
    plt.scatter(y1, x1)
    plt.gca().invert_yaxis()
    #plt.yticks(rotation=90)
    plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
    plt.xlabel('pmid count', position=(0.5, 1.08))
    ax = plt.gca()
    ax.xaxis.set_label_position('top')
    plt.savefig(db_path+'/'+gene+'_mesh_plot_'+timestamp+'_total.png',dpi=120, bbox_inches = "tight")
    plt.close()

def get_studytype_data(ref, gene, save_path):
    dfbib = pd.DataFrame(ref)
    dfbib.pubmed_id = dfbib.pubmed_id.astype('str')
    if 'publication_types' in dfbib.columns and len(dfbib['publication_types'])>1:
        dfbib_studyty = dfbib[['pubmed_id','publication_types']].dropna().reset_index(drop=True)

        #PMID-Studytype table build:
        studytype_list = []
        for i in range(len(dfbib_studyty)):
            for studytype in dfbib_studyty['publication_types'][i]:
                out = dfbib_studyty['pubmed_id'][i], studytype
                studytype_list.append(out)
        studytype_df = pd.DataFrame(studytype_list).rename(columns={0: 'pmid',1:'study_type'})
        mask_st = studytype_df['study_type'].str.contains('Research Support|Journal Article')
        studytype_df_less = studytype_df[~mask_st].reset_index(drop=True)

        mask_lessing = studytype_df['pmid'].isin(studytype_df_less['pmid'])
        studytype_df_diff = studytype_df[~mask_lessing].reset_index(drop=True)
        studytype_df_diff['study_type2'] = 'Unknown'
        studytype_df_diff = studytype_df_diff[['pmid','study_type2']].rename(columns={'study_type2':'study_type'}).drop_duplicates().reset_index(drop=True)
        studytype_df_concat = pd.concat([studytype_df_less, studytype_df_diff], ignore_index=True)

        studytype_df_concat.to_csv(save_path+'/'+gene+'_lit1pmid_studytype.csv')
    else:
        print(gene+' no publication_types in nbib')
        pass