#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[ ]:


import os
import glob
import pyperclip
import pandas as pd
import requests as rq
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

#from Bio import Entrez
#Entrez.email = "your_email@example.com"


# In[ ]:


#Only for for Google Colab
import sys
if 'google.colab' in sys.modules:
    from google.colab import drive
    drive.mount('/content/drive')
    os.chdir('/content/drive/MyDrive/GRPM db/GRPM db - Vitam')


# In[ ]:


#If needed, import openai
import openai
openai.api_key = "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"
def chatWithGPT(prompt):
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "user", "content": prompt}])
    return print(completion.choices[0].message.content)

def fixMyCode(code):
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "user", "content": "find error in my python script below and fix it: " + code}])
    return print(completion.choices[0].message.content)

# Initialize the conversation (text-davinci-003)----------------
conversation = [
    {"role": "assistant", "content": "You are a helpful assistant expert in science and informatics."},
]
def expand_conversation(message):
    conversation.append({'role': 'user', 'content': message})

def build_prompt(conversation):
    prompt = ""
    for message in conversation:
        role = message['role']
        content = message['content']
        prompt += f'{role}: {content}\n'
    return prompt

def send_message(message):
    expand_conversation(message)
    response = openai.Completion.create(
        engine='text-davinci-003',
        prompt = build_prompt(conversation),
        temperature=0.7,
        max_tokens=2000,
        n=1,
        stop=None
    )
    conversation.append({'role': 'assistant', 'content': response.choices[0].text.strip()})
    # Print the assistant's response
    print(response.choices[0].text.strip())

# Initialize the conversation (gpt-3.5-turbo)-----------------

# Conversation history
conversation_gpt = [
    {"role": "assistant", "content": "You are a helpful assistant expert in science and informatics."},
]

def expand_conversation_gpt(message):
    conversation_gpt.append({"role": "user", "content": message})

def build_messages(conversation):
    messages = []
    for message in conversation:
        messages.append({"role": message["role"], "content": message["content"]})
    return messages

def send_message_gpt(message):
    expand_conversation_gpt(message)
    messages = build_messages(conversation_gpt)
    response = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=messages,
        max_tokens= 2000
    )
    # Add the assistant's reply to the conversation
    conversation_gpt.append({"role": "assistant", "content": response.choices[0].message.content})

    # Print the assistant's response
    print(response.choices[0].message.content)


# In[ ]:


# Example usage
send_message_gpt("explain me nutrigenetics")


# In[ ]:


message = "explain me nutrigenetics" #input()
send_message(message)


# # Load Human Gene list
# protein_coding_genes = pd.read_csv('human_genes_repo/H_GENES_proteincoding_genes.csv')
# IG_TR_genes          = pd.read_csv('human_genes_repo/H_GENES_IGTR_genes.csv')
# RNA_genes            = pd.read_csv('human_genes_repo/H_GENES_RNA_genes.csv')
# pseudo_genes         = pd.read_csv('human_genes_repo/H_GENES_pseudo_genes.csv')
# misc_genes           = pd.read_csv('human_genes_repo/H_GENES_misc_genes.csv')
# 
# #First: protein coding genes:
# protein_coding_genes_list = protein_coding_genes['Gene name'].dropna().tolist()
# RNA_genes_list = RNA_genes['Gene name'].dropna().tolist()
# pseudo_genes_list = pseudo_genes['Gene name'].dropna().tolist()

# GRPM_report = pd.read_csv('grpm_db/GRPM_report.csv', index_col=0).T.reset_index()
# GRPM_report

# # Import GRPM Report

# In[ ]:


#set surce database:-----------------------
db_tag = 'pcg'
db_name = 'grpm_db_'+db_tag
# pcg    = protein coding genes = grpm_db
# rna    = rna genes            = grpm_db_rna
# pseudo = pseudogenes          = grpm_db_pseudo
# pseudo = pseudogenes          = grpm_db_pseudo
#------------------------------------------

# IMport GRPM_report
GRPM_report = pd.read_csv(db_name+'/GRPM_report.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})

data_types = {
    'gene': str,
    'ncbi_dbsnp': int,
    'lit2_variant': int,
    'lit2_variant_norsid': int,
    'lit2_rsid': int,
    'lit2_rsid_plus1': int,

    'lit1_rsid': int,
    'lit1_rsid_pmid_plus1': int,
    'lit1_pmid': int,
    'lit1_pmid_pmid_plus1': int,
    'pubmed_pmid_query': int,
    'nbib_objects': int,
    'nbib_objects_withdescriptors': int,
    'pubmed_pmid': int,
    'pubmed_pmid_withmesh': int,
    'pubmed_pmidmesh': int,
    'pubmed_mesh_qualifier_major': int,
    'pubmed_mesh': int,
    'rsid_pmid10': int,
    'rsid_pmid50': int,
    'rsid_pmid100': str,
    'top10mesh_all': str,
    'top10rsid_all': str,
    'pubmed_runtime': str,
    'total_runtime': str,
    'time_stamp': str
    }
GRPM_report = GRPM_report.astype(data_types)

GRPM_report.T
GRPM_report.dtypes
GRPM_report


# # Import GRPMX Reports

# In[ ]:


#Check survey:
folder_path = "ref-mesh-archive"  # Replace with the actual folder path

# Create a file path pattern to match CSV files
file_pattern = os.path.join(folder_path, "*.csv")

# Use glob to get a list of file paths matching the pattern
csv_files = glob.glob(file_pattern)
csv_files_name = []
# Print the list of CSV files
for file in csv_files:
    file_name = os.path.basename(file)
    csv_files_name.append(file_name)

print('Available reference mesh lists:')
csv_files_df = pd.Series(csv_files_name)
file_list = csv_files_df[csv_files_df.str.contains('ref_mesh_')]
file_names = [x.split('mesh_')[1].split('.csv')[0] for x in file_list]
file_names


# chatWithGPT("in python how to extract from a list everiting between '_' and '.csv'?")

# In[ ]:


# Import Survey
db_tag = 'pcg'
# pcg    = protein coding genes = grpm_db_pcg
# rna    = rna genes            = grpm_db_rna
# pseudo = pseudogenes          = grpm_db_pseudo

# Create an empty list to store folder names
folder_names = []
current_dir = os.getcwd()
# Iterate over the directories in the workspace
for root, dirs, files in os.walk(current_dir):
    for dir_name in dirs:
        # Check if the folder name contains the string 'survey'
        if 'survey' in dir_name:
            folder_names.append(dir_name)

# Create a pandas Series from the list of folder names
folder_series = pd.Series(folder_names)
folder_series = folder_series.str.replace('grpm_survey_'+db_tag+'_','')
# Print the resulting Series
print(folder_series)


# In[ ]:


db_tag = 'pcg'

# choose directories to analyze:-------------
first, last = 0, 28
#first, last = 0, 3
#--------------------------------------------
#
#dir_names = ['directory_01', 'directory_02', 'directory_03', 'directory_04', 'directory_05', 'directory_06','directory_07', 'directory_08', 'directory_09', 'directory_10', 'directory_11', 'directory_12','directory_13', 'directory_14', 'directory_15', 'directory_16', 'directory_17', 'directory_18','directory_19', 'directory_20', 'directory_21', 'directory_22', 'directory_23', 'directory_24','directory_25', 'directory_26', 'directory_27', 'directory_28']
#
#for folder, name in zip(folder_series, dir_names[:len(folder_series)]):
#    director =  'grpm_survey_'+ db_tag + '_'+ folder
#    globals()[name] = director

# define survey directories---------------
if 1 > 2:
    directory_01 = 'grpm_survey_'+db_tag+'_'+'nutri'
    directory_02 = 'grpm_survey_'+db_tag+'_'+'neuro'
    directory_03 = 'grpm_survey_'+db_tag+'_'+'infect'
    directory_04 = 'grpm_survey_'+db_tag+'_'+'repro'
    directory_05 = 'grpm_survey_'+db_tag+'_'+'dmt2_ms'
    directory_06 = 'grpm_survey_'+db_tag+'_'+'vitam'
    directory_07 = 'grpm_survey_'+db_tag+'_'+'eat_taste'
    directory_08 = 'grpm_random_'+db_tag+'_'+'03'
    directory_09 = 'grpm_random_'+db_tag+'_'+'04'
    directory_10 = 'grpm_random_'+db_tag+'_'+'05'
else:
    directory_01 = 'grpm_survey_'+db_tag+'_'+'nutri'
    directory_02 = 'grpm_survey_'+db_tag+'_'+'ob_bmi'
    directory_03 = 'grpm_survey_'+db_tag+'_'+'dmt2_ms'
    directory_04 = 'grpm_survey_'+db_tag+'_'+'cvd'
    directory_05 = 'grpm_survey_'+db_tag+'_'+'vitam'
    directory_06 = 'grpm_survey_'+db_tag+'_'+'eat_taste'
    directory_07 = 'grpm_survey_'+db_tag+'_'+'intol'
    directory_08 = 'grpm_survey_'+db_tag+'_'+'aller'
    directory_09 = 'grpm_survey_'+db_tag+'_'+'oxi_stress'
    directory_10 = 'grpm_survey_'+db_tag+'_'+'xeno'

#define random directories
directory_11 = 'grpm_random_'+db_tag+'_'+'04'
directory_12 = 'grpm_random_'+db_tag+'_'+'05'
directory_13 = 'grpm_random_'+db_tag+'_'+'06'
directory_14 = 'grpm_random_'+db_tag+'_'+'07'
directory_15 = 'grpm_random_'+db_tag+'_'+'08'
directory_16 = 'grpm_random_'+db_tag+'_'+'09'
directory_17 = 'grpm_random_'+db_tag+'_'+'00'
directory_18 = 'grpm_random_'+db_tag+'_'+'norep1'
directory_19 = 'grpm_random_'+db_tag+'_'+'norep2'
directory_20 = 'grpm_random_'+db_tag+'_'+'norep3'
directory_21 = 'grpm_random_'+db_tag+'_'+'norep4'
directory_22 = 'grpm_random_'+db_tag+'_'+'norep5'
directory_23 = 'grpm_random_'+db_tag+'_'+'norep6'
directory_24 = 'grpm_random_'+db_tag+'_'+'norep7'
directory_25 = 'grpm_random_'+db_tag+'_'+'norep8'
directory_26 = 'grpm_random_'+db_tag+'_'+'norep9'
directory_27 = 'grpm_random_'+db_tag+'_'+'norep10'
directory_28 = 'grpm_random_'+db_tag+'_'+'norep11'

#set directories for this comapration:
directories = [directory_01,directory_02,directory_03,directory_04,directory_05,directory_06,directory_07,directory_08,directory_09,directory_10,directory_11,directory_12,directory_13,directory_14,directory_15,directory_16,directory_17,directory_18,directory_19, directory_20,directory_21,directory_22,directory_23,directory_24,directory_25,directory_26,directory_27,directory_28 ]
directories = directories[first:last]
survey_l = pd.Series(directories)[pd.Series(directories).str.contains('survey')]


# import GRPMX reports------------------
repo_list_names = ['GRPMX_report_01', 'GRPMX_report_02', 'GRPMX_report_03', 'GRPMX_report_04', 'GRPMX_report_05', 'GRPMX_report_06', 'GRPMX_report_07', 'GRPMX_report_08', 'GRPMX_report_09', 'GRPMX_report_10', 'GRPMX_report_11', 'GRPMX_report_12', 'GRPMX_report_13', 'GRPMX_report_14', 'GRPMX_report_15', 'GRPMX_report_16', 'GRPMX_report_17', 'GRPMX_report_18', 'GRPMX_report_19', 'GRPMX_report_20','GRPMX_report_21','GRPMX_report_22','GRPMX_report_23','GRPMX_report_24','GRPMX_report_25','GRPMX_report_26','GRPMX_report_27','GRPMX_report_28' ]
repo_list_names = repo_list_names[first:last]

repo_list = []
for directory, name in zip(directories, repo_list_names):
    dataframe = pd.read_csv(directory +'/GRPMX_report.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
    globals()[name] = dataframe
    repo_list.append(dataframe)

    print(name,'df loaded')

for i, repo in enumerate(repo_list):
  data_types = {
      'gene': str,
      'reference_mesh': int,
      'starting_pmidmesh': int,
      'starting_pmid': int,
      'starting_mesh': int,
      'starting_rsid': int,
      'matching_pmidmesh': int,
      'matching_pmids': int,
      'matching_mesh': int,
      'matching_rsid': int,
      'dropped_rsid': int,
      'matching_mesh_ratio': float,
      'matching_pmids_ratio': float,
      'matching_pmidmesh_ratio': float,
      'matching_rsid_ratio': float,
      'matching_rsid_pmid10': str,
      'matching_rsid_pmid100': str,
      'matching_top10mesh': str,
      'matching_top10rsid': str
  }
  repo_list[i] = repo.astype(data_types)

master_df = pd.DataFrame({'directories': directories, 'repo_list': repo_list_names})


# In[ ]:


len(master_df)
send_message_gpt('how to set on a seborn heatmap the max heat value for every row?')


# In[ ]:


send_message('how to unite two lists?')


# ### add labels to master df

# In[ ]:


add_label = True
if add_label == True:
    titles = [
        "Genaral Nutrition",
        "Obesity, Weight Control and Compulsive Eating",
        "Cardiovascular Health and Lipid Metabolism",
        "Diabetes Mellitus Type II and Metabolic Syndrome",
        "Vitamin and Micronutrients Metabolism and Deficiency-Related Diseases",
        "Eating Behavior and Taste Sensation",
        "Food Intolerances",
        "Food Allergies",
        "Diet-induced Oxidative Stress",
        "Xenobiotics Metabolism"]
    add_rand = ['random'] * (len(master_df)- len(titles))
    titles.extend(add_rand)
    master_df['label'] = titles
    #new_order = ['label', 'directories', 'repo_list']
    #master_df  = master_df[new_order]
    last_column = master_df.columns[-1]
    last_column_data = master_df.pop(last_column)
    master_df.insert(0, last_column, last_column_data)
master_df

master_df


# In[ ]:


repo_list[0].dtypes
#enumerate(repo_list)


# #check
# type(GRPMX_report_01.reference_mesh[0])
# display(GRPMX_report_20.matching_pmids)
# display(repo_list[19].matching_pmids)

# #controllo
# range(len(directories)-1)
# GRPMX_report_3 = pd.read_csv(directory_3+'/GRPMX_report.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
# GRPMX_report_3

# # import GRPMX data

# In[ ]:


# import GRPMX data
grpm_out_names = ['GRPMX_01', 'GRPMX_02', 'GRPMX_03', 'GRPMX_04', 'GRPMX_05', 'GRPMX_06', 'GRPMX_07', 'GRPMX_08', 'GRPMX_09', 'GRPMX_10', 'GRPMX_11', 'GRPMX_12', 'GRPMX_13', 'GRPMX_14', 'GRPMX_15', 'GRPMX_16', 'GRPMX_17', 'GRPMX_18', 'GRPMX_19', 'GRPMX_20','GRPMX_21', 'GRPMX_22','GRPMX_23','GRPMX_24','GRPMX_25','GRPMX_26','GRPMX_27','GRPMX_28' ]
grpm_out_names = grpm_out_names[first:last]

grpm_out_list = []
for directory, name in zip(directories, grpm_out_names):
    dataframe = pd.read_csv(directory +'/grpmx_filtered_output.csv',index_col=0)
    globals()[name] = dataframe
    grpm_out_list.append(dataframe)
    print(name,'df loaded')

master_df['grpm_out_list'] = grpm_out_names

for i in range(len(directories)):
    master_df.at[i, 'grpm_genes'] = grpm_out_list[i].gene.nunique()
    master_df.at[i, 'grpm_pmid'] = grpm_out_list[i].pmids.nunique()
    master_df.at[i, 'grpm_mesh'] = grpm_out_list[i].mesh.nunique()

print('grpm_out list len:',len(grpm_out_list))
#for i in range(len(directories)):
#    print(master_df.directories.iloc[i], 'gene number: ', grpm_out_list[i].gene.nunique())
master_df[master_df.directories.str.contains('survey')]


# In[ ]:


#choose mesh to get rid
check = 1
print(master_df.label[check])
grpm_out_list[check].groupby('mesh').describe()


# In[ ]:


#PURGE UNWANTED MESH from survey.csv
grpm_out_list[check].mesh


# In[ ]:


#PURGE UNWANTED MESH from survey.csv
id = check
save_clean = False
tag = '_clean'
exact_match = False

dff = grpm_out_list[check]

#clean mesh ref
get_rid_list = ['Polymor','Genetic']#['Fistula', 'Neoplasm','Labor', 'Chick']
get_rid_ex = ['Polymorphism, Single Nucleotide', 'Genetic Predisposition to Disease']
dff_rid = grpm_out_list[id].copy()

if exact_match == False:
    for get_rid in get_rid_list:
        mask = dff['mesh'].str.contains(get_rid)
        dff_rid = dff_rid[-mask]
else:
    mask = dff['mesh'].isin(get_rid_ex)
    dff_rid = dff_rid[-mask]

if save_clean == True:
    dff_rid.to_csv(directories[id]+'/grpmx_filtered_output'+tag+'.csv')
dff_rid


# In[ ]:


#cleaned = pd.read_csv(directories[check]+'/grpmx_filtered_output'+tag+'.csv', index_col=0)


# In[ ]:


master_df


# #check
# display(pd.read_csv(directory_20 +'/grpmx_filtered_output.csv',index_col=0))
# display(grpm_out_list[19])

# # Add "interest index" to report:

# In[ ]:


repo_list[0].columns


# In[ ]:


## Add "interest index" to report:
replace_all = False

#--------------------------------------
def add_interest_index(df):
    max_match_pmids = int(df['matching_pmids'].max())
    # duplicate df
    df_int = df#.copy()
    df_int['matching_pmids_index'] = round((df_int['matching_pmids'] / max_match_pmids), 5) #matching_pmids = Pm(g)
    df_int['interest_index'] = round(df_int['matching_pmids_index'] * df_int['matching_pmids_ratio'], 5)
    df_int['matching_mesh_ref_ratio'] = round(df_int['matching_mesh'] / df_int['reference_mesh'], 5)
    new_col_order = ['gene', 'interest_index', 'reference_mesh', 'starting_pmidmesh',
                     'starting_pmid', 'starting_mesh', 'starting_rsid',
                     'matching_pmids_index','matching_pmidmesh', 'matching_pmids',
                     'matching_mesh', 'matching_rsid',
                     'dropped_rsid', 'matching_mesh_ratio',
                     'matching_mesh_ref_ratio', 'matching_pmids_ratio',
                     'matching_pmidmesh_ratio', 'matching_rsid_ratio',
                     'matching_rsid_pmid10', 'matching_rsid_pmid100', 'matching_top10mesh',
                     'matching_top10rsid']
    df_int = df[new_col_order]
    return df_int
#-------------------------------------

##GRPMX_report_int

report_int_names = ['GRPMX_report_int_01','GRPMX_report_int_02','GRPMX_report_int_03','GRPMX_report_int_04','GRPMX_report_int_05','GRPMX_report_int_06','GRPMX_report_int_07','GRPMX_report_int_08','GRPMX_report_int_09','GRPMX_report_int_10','GRPMX_report_int_11','GRPMX_report_int_12','GRPMX_report_int_13','GRPMX_report_int_14','GRPMX_report_int_15','GRPMX_report_int_16','GRPMX_report_int_17','GRPMX_report_int_18','GRPMX_report_int_19','GRPMX_report_int_20','GRPMX_report_int_21','GRPMX_report_int_22','GRPMX_report_int_23','GRPMX_report_int_24','GRPMX_report_int_25','GRPMX_report_int_26','GRPMX_report_int_27','GRPMX_report_int_28' ]
report_int_names = report_int_names[first:last]
master_df['repo_int_list'] = report_int_names

repo_int_list = []
for repo, name, directory in zip(repo_list, report_int_names, directories):
    if os.path.exists(directory + '/GRPMX_report_int.csv') and replace_all == False:
        dataframe = pd.read_csv(directory + '/GRPMX_report_int.csv')
        data_types = {
            'gene': str,
            'reference_mesh': int,
            'starting_pmidmesh': int,
            'starting_pmid': int,
            'starting_mesh': int,
            'starting_rsid': int,
            'matching_pmidmesh': int,
            'matching_pmids': int,
            'matching_mesh': int,
            'matching_rsid': int,
            'dropped_rsid': int,
            'interest_index': float,
            'matching_mesh_ratio': float,
            'matching_mesh_ref_ratio': float,
            'matching_pmids_ratio': float,
            'matching_pmidmesh_ratio': float,
            'matching_rsid_ratio': float,
            'matching_pmids_index': float,
            'matching_rsid_pmid10' : str,
            'matching_rsid_pmid100': str,
            'matching_top10mesh': str,
            'matching_top10rsid': str
            }
        dataframe = dataframe.astype(data_types)
        globals()[name] = dataframe
        repo_int_list.append(dataframe)
        print(name,'df loaded')

    else:
        dataframe = add_interest_index(repo)
        globals()[name] = dataframe
        dataframe.to_csv(directory+'/GRPMX_report_int.csv', index=0)
        repo_int_list.append(dataframe)
        print(name,'df created')

if 1 > 2:
    for repo in repo_int_list:
        repo.drop(['starting_pmidmesh', 'starting_pmid','starting_mesh', 'starting_rsid'],axis=1, inplace=True)

master_df


# In[ ]:


i = 24
type(repo_int_list[i].matching_pmids[i])
GRPMX_report_int_05


# #check
# display(GRPMX_report_int_18)
# display(repo_int_list[17])

# # Analyse GRPMX Repo Data

# ## Sort for interest index

# In[ ]:


#Sort for interest index -------------
print('GRPMX table sort by interest index',directory,':')

def sortby_int(df):
    df_int_sort = df.sort_values(by='interest_index',ascending=False).reset_index(drop=True)
    return df_int_sort

repo_int_sort_name = ['GRPMX_report_int_01_sort','GRPMX_report_int_02_sort','GRPMX_report_int_03_sort','GRPMX_report_int_04_sort','GRPMX_report_int_05_sort','GRPMX_report_int_06_sort','GRPMX_report_int_07_sort','GRPMX_report_int_08_sort','GRPMX_report_int_09_sort','GRPMX_report_int_10_sort','GRPMX_report_int_11_sort','GRPMX_report_int_12_sort','GRPMX_report_int_13_sort','GRPMX_report_int_14_sort','GRPMX_report_int_15_sort','GRPMX_report_int_16_sort','GRPMX_report_int_17_sort','GRPMX_report_int_18_sort','GRPMX_report_int_19_sort','GRPMX_report_int_20_sort','GRPMX_report_int_21_sort','GRPMX_report_int_22_sort','GRPMX_report_int_23_sort','GRPMX_report_int_24_sort','GRPMX_report_int_25_sort','GRPMX_report_int_26_sort','GRPMX_report_int_27_sort','GRPMX_report_int_28_sort' ]
repo_int_sort_name = repo_int_sort_name[first:last]
master_df['repo_int_sort_list'] = repo_int_sort_name

repo_int_sort_list = []
for directory, name, repo in zip(directories, repo_int_sort_name, repo_int_list):
    dataframe = sortby_int(repo)
    globals()[name] = dataframe
    repo_int_sort_list.append(dataframe)
    print(name,'df sorted')
print('')
for i in range(len(directories)):
    print(master_df.directories.iloc[i], 'gene number: ', len(grpm_out_list[i].gene))
master_df.T


# In[ ]:


master_df.directories.iloc[0]
repo_int_sort_list[0]


# #check
# display(GRPMX_report_int_01_sort)#[['gene','matching_pmids_ratio','interest_index']]
# repo_int_sort_list[0]

# ## Create datasheet with all Reports + grpm data

# In[ ]:


repo_int_sort_list[0]


# In[ ]:


#Create datasheet-----------------------------------------
value   = 'interest_index'
value_2 = 'matching_mesh'
value_x = 'matching_mesh_ref_ratio'
value_3 = 'matching_pmids'
value_4 = 'matching_pmids_ratio'
out_1 = 'gene'
out_2 = 'pmids'
out_3 = 'mesh'

cifr = 2
f_gene, l_gen2 = 0, 200 # <== set gene-range
comp_df_trial = pd.DataFrame()
comp_df_trial = pd.DataFrame({'dir':directories,
                              'int_mean'           :0,
                              'mesh_mean'          :0,
                              'mesh_ref_ratio_mean':0,
                              'pmids_mean'         :0,
                              'pmids_ratio_mean'   :0,
                              out_1: 0,out_2: 0,out_3: 0 })

for i in range(len(comp_df_trial)):
    repost = repo_int_sort_list[i]
    grpm_out = grpm_out_list[i]
    comp_df_trial.loc[i, 'int_mean' ] = round(repost[value  ].iloc[  f_gene:l_gen2].mean(), cifr+2)
    comp_df_trial.loc[i, 'mesh_mean'] = round(repost[value_2].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial.loc[i, 'mesh_ref_ratio_mean'] = round(repost[value_x].iloc[f_gene:l_gen2].mean(), cifr+2)
    comp_df_trial.loc[i, 'pmids_mean'         ] = round(repost[value_3].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial.loc[i, 'pmids_ratio_mean'   ] = round(repost[value_4].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial.loc[i, out_1  ] = grpm_out[out_1].nunique()
    comp_df_trial.loc[i, out_2  ] = grpm_out[out_2].nunique()
    comp_df_trial.loc[i, out_3  ] = grpm_out[out_3].nunique()

#comp_df_trial[[out_1,out_2,out_3]] = comp_df_trial[[out_1,out_2,out_3]].astype(int)
print('Survey datasheet - gene mean values for InterestIndex Top',l_gen2,'genes')
comp_df_trial#.to_csv('nutri+neuro+17random.csv')
#pyperclip.copy(comp_df_trial.to_csv())


# In[ ]:


repo_int_sort_list[5][['gene','interest_index','matching_pmids','matching_mesh','matching_mesh_ref_ratio','reference_mesh']]


# In[ ]:


#Create datasheet-----------------------------------------
value = 'interest_index'
value_2 = 'matching_mesh'
value_x = 'matching_mesh_ref_ratio'
value_3 = 'matching_pmids'
value_4 = 'matching_pmids_ratio'
out_1 = 'gene'
out_2 = 'pmids'
out_3 = 'mesh'

cifr = 2
f_gene, l_gen2 = 0, len(repo_list[0]) # <== set gene-range

comp_df_trial_all = pd.DataFrame({'dir':directories,
                                  'int_mean'           :0,
                                  'mesh_mean'          :0,
                                  'mesh_ref_ratio_mean':0,
                                  'pmids_mean'         :0,
                                  'pmids_ratio_mean'   :0,
                                  out_1: 0,out_2: 0,out_3: 0 })

for i in range(len(comp_df_trial)):
    repost = repo_int_sort_list[i]
    grpm_out = grpm_out_list[i]
    comp_df_trial_all.loc[i, 'int_mean' ] = round(repost[value].iloc[  f_gene:l_gen2].mean(), cifr+2)
    comp_df_trial_all.loc[i, 'mesh_mean'] = round(repost[value_2].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial_all.loc[i, 'mesh_ref_ratio_mean'] = round(repost[value_x].iloc[f_gene:l_gen2].mean(), cifr+2)
    comp_df_trial_all.loc[i, 'pmids_mean'] = round(repost[value_3].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial_all.loc[i, 'pmids_ratio_mean'] = round(repost[value_4].iloc[f_gene:l_gen2].mean(), cifr)
    comp_df_trial_all.loc[i, out_1  ] = grpm_out[out_1].nunique()
    comp_df_trial_all.loc[i, out_2  ] = grpm_out[out_2].nunique()
    comp_df_trial_all.loc[i, out_3  ] = grpm_out[out_3].nunique()

#comp_df_trial[[out_1,out_2,out_3]] = comp_df_trial[[out_1,out_2,out_3]].astype(int)
print('Survey datasheet - gene mean values for InterestIndex Top',l_gen2,'genes')
comp_df_trial_all#.to_csv('nutri+neuro+17random.csv')


# ## Calculate percentile95

# #Assuming you have a DataFrame named 'df' with a column named 'column_name'
# data = repo_int_sort_list[0]['interest_index']
# 
# #Calculate the 95th percentile using numpy
# percentile_95 = np.percentile(data, 95)
# 
# print("95th percentile:", percentile_95)
# 

# ### Filter grpmx output for quantile 0.95

# for i in range(len(directories)):
#     print(master_df.directories.iloc[i], 'interest index 95 quantile', round(repo_int_sort_list[i].interest_index.quantile(0.95),5))

# #load my GRPMx Data
# dir = 0
# df_grpmx = pd.read_csv(directories[dir]+'/grpmx_filtered_output.csv', index_col=0)
# df_grpmx_repo = pd.read_csv(directories[dir]+'/GRPMX_report_int.csv')
# 
# #filter for 0.95 quantile
# df_grpmx_repo_95 = df_grpmx_repo[df_grpmx_repo.interest_index >= df_grpmx_repo.interest_index.quantile(0.95)]
# df_grpmx_95 = df_grpmx[df_grpmx.gene.isin(df_grpmx_repo_95.gene)]
# print('df_grpmx_95 genes:', df_grpmx_95.gene.nunique())
# df_grpmx_95

# ### Calculate quantile 0.95 for all datasets

# In[ ]:


for i in range(len(comp_df_trial_all)):
    repo = repo_int_sort_list[i]
    grpm_out = grpm_out_list[i]
    percentile95 = repo.interest_index.quantile(0.95)

    comp_df_trial_all.loc[i, 'percentile_95'] = percentile95
    dummy = repost[repost.interest_index > percentile95]
    comp_df_trial_all.loc[i, 'genes_95'] = len(dummy)
    comp_df_trial_all.loc[i, 'int_min_95'] = dummy.interest_index.min()

print(comp_df_trial_all[['dir','percentile_95','genes_95', 'int_min_95']])
#memo: la soglia deve essere basata sull interest index non sul percentile
#comp_df_trial_all.to_csv('report_df.csv', index=0)


# In[ ]:


add_label = True
if add_label == True:
    #titles = master_df['label']
    comp_df_trial_all['label'] = master_df['label']
    columns = comp_df_trial_all.columns.tolist()
    columns = [columns[-1]] + columns[:-1]
    comp_df_trial_all = comp_df_trial_all[columns]
comp_df_trial_all


# In[ ]:


pd.Series(directories)[pd.Series(directories).str.contains('survey')]


# In[ ]:


directories[:len(survey_l)]


# In[ ]:


## Survey Density Rugged Plot

f_gene, l_gene = 0, 200 # set top gene survey_l

color_labels = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black', '#D2B48C', '#A52A2A', '#800080']
color_labels = color_labels[:len(survey_l)]
max_y = 0

for n, i, col in zip(range(len(survey_l)), directories[:len(survey_l)], color_labels):
    dat = repo_int_sort_list[n].interest_index.iloc[f_gene:l_gene]

    sns.set_style('whitegrid')
    sns.kdeplot(np.array(dat), color = col, bw_method=0.5)
    sns.rugplot(np.array(dat), color = col, label=directories[n])
    # Add vertical line for the 95th percentile
    percentile_95 = dat.quantile(0.95)
    plt.axvline(x = percentile_95, color= col, linestyle='--')#, label='95th percentile')

    #density = sns.kdeplot(np.array(dat)).get_lines()[0].get_data()
    #if max(density[1]) > max_y:
    #    max_y = max(density[1])

#plt.yscale('log')
plt.title('GRPMX interest index plot: '+directories[n])
plt.xlim(0)
plt.ylim(0, 30)
#plt.yscale('log')
plt.legend()
print('GRPMX interest index plot:', directories[n])
plt.show()


# In[ ]:


max_y


# In[ ]:


print('Statistics:')
f_gene, l_gen2 = 0, len(repo_list[0])
for i in range(len(directories)):
    repost = repo_int_sort_list[i]
    grpm_out = grpm_out_list[i]

    dir = directories[i]
    print(dir   ,'int_ind mean:',round(repost.interest_index.iloc[first:last].mean(), 6))
#GRPMX_report_int_7_sort.gene


# In[ ]:


#Sort for interest index -------------
for repo, dir in zip(repo_int_sort_list,directories):
    print('GRPMX table sort by interest index',dir,':')
    print(' ',repo.gene.iloc[0])
    print(' ',repo.gene.iloc[1])
    print(' ',repo.gene.iloc[2])
    print(' ',repo.gene.iloc[3])


# #show GRPM database
# GRPM_report_sort = GRPM_report.sort_values(by='pubmed_pmid', ascending = False).reset_index(drop=True)
# print('GRPM db table sort by pmid number:')
# GRPM_report_sort

# ### Density plot: GRPMX repo

# In[ ]:


## Survey Desnity Rugged Plot
f_gene, l_gen2 = 0, 100#, 20000
survey_l = len(pd.Series(directories)[pd.Series(directories).str.contains('survey')])

datas_names = ['data_10',
               'data_20',
               'data_30',
               'data_40',
               'data_50',
               'data_60',
               'data_70',
               'data_80',
               'data_90',
               'data_10',]
datas = []
for name, repo in zip(datas_names[first:last],repo_int_sort_list):
    dataframe = repo.interest_index.iloc[f_gene:l_gen2]
    globals()[name] = dataframe
    datas.append(dataframe)

color_labels = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black', '#D2B48C', '#A52A2A', '#800080']
#color_labels = ['red', 'green', 'blue', 'cyan', 'magenta']
color_labels = color_labels[0:survey_l]
start, end = 0, survey_l

plt.figure(figsize=(10, 4))
#plot randoms in grey
for dat in datas[survey_l+1:]:
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(dat), bw_method=0.5, color='#808080')
    sns.rugplot(np.array(dat), color='#808080')

#plot survey above in color
for dat, col, dir  in zip(datas[0:survey_l],
                    color_labels,
                    directories[0:survey_l]
                    ):
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(dat), bw_method=0.5, color=col, label= dir )
    sns.rugplot(np.array(dat), color=col)

#plt.xscale('log')
plt.title('Ig = (Pm/P)g')
plt.title('Interest index = a^2/bc')
#plt.xscale('log')
plt.ylim(0, 20)
plt.xlim(0)
plt.legend()
print('GRPMX interest index plot:', directory)
frame = {'directory': directories[start:len(color_labels)], 'color code':color_labels}
print(pd.DataFrame(frame))
plt.show()


# In[ ]:


master_df


# ### Density plot: GRPMX repo - General

# In[ ]:


#settings-----
aggregate_randoms = True
survey_l = len(pd.Series(directories)[pd.Series(directories).str.contains('survey')])
#-------------

#full plot
datas_name = ['data_01','data_02','data_03','data_04','data_05','data_06','data_07','data_08','data_09','data_10', 'data_11', 'data_12', 'data_13', 'data_14', 'data_15', 'data_16', 'data_17', 'data_18', 'data_19','data_20','data_21','data_22','data_23','data_24','data_25','data_26','data_27','data_28']
start, end = 0, len(datas_name)

# gene list slice:
f_gene, l_gen2 = 0, 200#, 20000
def slice_int(df):
    df_int_slice = df.interest_index.iloc[f_gene:l_gen2]
    return df_int_slice

# Choose dataset to comapre, slice datas list:
datas = []
for directory, name, repo in zip(directories[start:end], datas_name[start:end], repo_int_sort_list[start:end]):
    dataframe = slice_int(repo)
    globals()[name] = dataframe
    datas.append(dataframe)

# set survey series from master_df
survey_l = len(master_df[master_df.directories.str.contains('survey')])
color_labels = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black', '#D2B48C', '#A52A2A', '#D8BFD8']
color_labels = color_labels[:survey_l]

if len(datas[survey_l+1:end]) != 0:
    cluster_rand = pd.concat(datas[len(color_labels)+1:end], axis=1)

# Desity Rugged Plot--------------------------
plt.figure(figsize=(10, 4))

#plot randoms in grey
if len(datas[survey_l+1:end]) != 0:
    if aggregate_randoms == True:
        dat = cluster_rand.mean(axis=1)
        sns.set_style('whitegrid')
        sns.kdeplot(np.array(dat), bw_method=0.5, color='#808080')
        sns.rugplot(np.array(dat), color='#808080', label='random mean')
    else:
        for dat in datas[survey_l+1:end]:
            sns.set_style('whitegrid')
            sns.kdeplot(np.array(dat), bw_method=0.5, color='#808080')
            sns.rugplot(np.array(dat), color='#808080')

#plot survey above in color
for dat, col, dir in zip(datas[start:survey_l],
                         color_labels,
                         directories[start:survey_l]):
    sns.set_style('whitegrid')
    sns.kdeplot(np.array(dat), bw_method=0.5, color=col)
    sns.rugplot(np.array(dat), color=col, label=dir)

#plt.xscale('log')
plt.title('Interest index plot for survey VS random, top '+str(l_gen2)+' genes')
plt.xscale('log')
plt.ylim(0, 20)
plt.xlim(10 ** -2.3)
plt.legend()

print('GRPMX interest index plot:')
#frame = {'directory': directories[start:end], 'color code':color_labels[start:end]}
frame = {'directory': directories[start:survey_l], 'color code':color_labels}
print(pd.DataFrame(frame))
plt.show()


# In[ ]:


master_df


# ## Plotting comparison over Master reference mesh

# In[ ]:


master_df


# In[ ]:


genes_l = 25 #set master top genes
n = 0 #choose master gene ref from master_df

master_genes_repo = master_df.directories.iloc[n]
master_genes = repo_int_sort_list[n][:genes_l].gene
len(master_genes)

#slice all report for master_genes
#pyperclip.copy(str(master_df.repo_int_sort_list.to_list()))

def master_sampler(df):
    df_int_master_sample = df[df.gene.isin(master_genes)]
    return df_int_master_sample

master_list_names = ['GRPMX_report_int_01_master','GRPMX_report_int_02_master','GRPMX_report_int_03_master','GRPMX_report_int_04_master','GRPMX_report_int_05_master','GRPMX_report_int_06_master','GRPMX_report_int_07_master','GRPMX_report_int_08_master','GRPMX_report_int_09_master','GRPMX_report_int_10_master','GRPMX_report_int_11_master', 'GRPMX_report_int_12_master', 'GRPMX_report_int_13_master', 'GRPMX_report_int_14_master', 'GRPMX_report_int_15_master', 'GRPMX_report_int_16_master', 'GRPMX_report_int_17_master', 'GRPMX_report_int_18_master', 'GRPMX_report_int_19_master', 'GRPMX_report_int_20_master','GRPMX_report_int_21_master','GRPMX_report_int_22_master','GRPMX_report_int_23_master','GRPMX_report_int_24_master','GRPMX_report_int_25_master','GRPMX_report_int_26_master','GRPMX_report_int_27_master','GRPMX_report_int_28_master']
master_list_names = master_list_names[first:last]
master_df['master_list_names'] = master_list_names
print('master gene ref:', master_genes_repo)

# Choose dataset to compare, slice datas list:
start, end = 0, 22
master_list = []
for directory, name, repo in zip(directories[start:end], master_list_names[start:end], repo_int_sort_list[start:end]):
    #if directory == master_df.directories.iloc[n]:
    #    pass
    #else:
        dataframe = master_sampler(repo)
        globals()[name] = dataframe
        master_list.append(dataframe)
        print(name,'df mastered')

#---------------------------------------------------

# set referece dataframe
n = n
referece = master_list[n].copy()
referece['gene_id'] = range(1, len(referece) + 1)
gene_id_df = referece[['gene','gene_id']]

# merge gene id into all dataframes
def merge_id(df):
    df_int_master_merge = pd.merge(df, gene_id_df, on='gene', how='left')
    return df_int_master_merge
print('')
# Choose dataset to compare, slice datas list:
new_master_list = []
for name, repo in zip(master_list_names, master_list):
    dataframe = merge_id(repo)
    globals()[name] = dataframe
    new_master_list.append(dataframe)
    print(name,'df mastered with gen id')


# In[ ]:


#Scatter Plotting
plot_random = False

# select special_list from master_df
use_special_list = True
select_dirs = [1,2,3]

# set survey series from master_df
survey_ser = master_df[master_df.directories.str.contains('survey')]
special_list = survey_ser
survey_l = len(survey_ser)
random_ser = master_df[master_df.directories.str.contains('random')]
random_ser_num = len(random_ser)

color_labels = ['green', 'blue', 'cyan', 'magenta', 'yellow', 'black', '#D2B48C', '#A52A2A', '#D8BFD8']
color_labels = color_labels[:survey_l]

dir_select = []
new_master_list_select = []
if use_special_list == True:
    special_list = select_dirs
    master = [new_master_list[index] for index in special_list]
    dire = [master_df.directories[index] for index in special_list]
    new_master_list_select.extend(master)
    dir_select.extend(dire)
else:
    new_master_list_select = new_master_list[2:survey_l]
    dir_select = directories[2:survey_l]


# Desity Rugged Plot
plt.figure(figsize=(15, 8))

if plot_random == True:
    #plot randoms in grey
    start, end = len(directories)-random_ser_num, len(directories)
    for dat in new_master_list[start:end]:
        # Create the dot plot using Seaborn
        x= dat.gene_id
        y = dat.interest_index
        # Plot the first series
        plt.scatter(x, y, color='#808080')

#plot surveies in color:
for dat, col, dir in zip(new_master_list_select,
                         color_labels[:len(new_master_list_select)],
                         dir_select):
    x = dat.gene_id
    y = dat.interest_index
    plt.scatter(x, y,  color=col, label=dir)

# plot master series in red:
x = new_master_list[n].gene_id
y = new_master_list[n].interest_index
plt.scatter(x, y,  color='red', label=master_df.directories[n])

# Set plot title and labels
plt.title('Interest index plot for survey VS random, top '+str(genes_l)+' genes')
plt.yscale('log')
plt.xlabel('gene id')
plt.ylabel('Interest Index')
plt.legend()

# Show the plot
print('GRPMX interest index plot:')
frame = {'directory': dir_select[:len(special_list)], 'color code':color_labels[:len(special_list)]}
print(pd.DataFrame(frame))
plt.show()


# ## Boxplotting

# In[ ]:


# group all random interest index data:
new_df = pd.DataFrame(columns= ['gene_id', 'interest_index'])

survey_ser = master_df[master_df.directories.str.contains('survey')]
survey_ser_num = len(survey_ser)

#start, end = len(survey_ser)+1, len(directories)
start, end = len(directories)-random_ser_num, len(directories)
for dat in new_master_list[start:end]:
    new_df = pd.concat([new_df, dat[['gene_id', 'interest_index']]],axis=0)

#Boxplotting:----------------
#set num of genes to plot
slice = 25
new_df_slice = new_df[new_df.gene_id <= slice]

if random_ser_num != 0:
    # Data for boxplot
    x_box = new_df_slice.gene_id
    y_box = new_df_slice.interest_index
    x_box_cat = pd.Categorical(x_box) # Convert x_box to categorical data

# Data for scatter plot
n = n
x_scatter = new_master_list[n][new_master_list[n].gene_id <= slice].gene_id
y_scatter = new_master_list[n][new_master_list[n].gene_id <= slice].interest_index
x_scatter_cat = pd.Categorical(x_scatter) # Convert x_scatter to categorical data


plt.figure(figsize=(12, 8))
if random_ser_num != 0:
    # Plot boxplot
    sns.boxplot(x=x_box_cat.codes, y=y_box)

# Plot scatter plot
plt.scatter(x=x_scatter_cat.codes, y=y_scatter, color='red', label= directories[n])

# Set x-axis tick labels
#plt.xticks(range(len(x_box)), x_box)
plt.title('Interest index plot for survey VS random, top '+str(genes_l)+' genes')
plt.yscale('log')
plt.legend()
# Show the plot
plt.show()


# #Create a DataFrame with sample data
# 
# #set num of genes to plot
# slice = 25
# new_df_slice = new_df[new_df.gene_id <= slice]
# 
# plt.figure(figsize=(12, 8))
# # Create a boxplot using seaborn
# sns.boxplot(x='gene_id', y='interest_index', data=new_df_slice)
# 
# #set reference survey dir
# n = 0
# x = new_master_list[n][new_master_list[n].gene_id <= slice].gene_id
# y = new_master_list[n][new_master_list[n].gene_id <= slice].interest_index
# plt.scatter(x, y,  color='red', label= directories[n])
# #sns.scatterplot(x=x, y=y)
# 
# #xr = new_master_list[1].gene_id
# #yr = new_master_list[1].interest_index
# #plt.scatter(xr, yr,  color='green')
# plt.yscale('log')
# plt.legend()
# #Display the plot
# plt.show()

# slice = 25
# new_df_slice = new_df[new_df.gene_id <= slice]
# pd.Categorical(new_df_slice.gene_id).codes

# Condivido un sunto dei risultati ottenuti fin ora in foma di scatter plot.
# 
# I risultati della comparazione di Diverse liste di Mesh sul Database GRPM creato sulla base dei dati LitVar e Pubmed:
# In rosso ho plottato gli "interest_index" dei Top 200 geni della lista mesh nutrizionali ordinati in ordine decrescente
# In verde i valori per gli stessi 200 geni però usando la lista di termini neurologici.
# In grigio i valori ottenuti usando 17 liste mesh random.
# 
# GRPMX interest index plot:
# directory color code
# grpm_survey_nutri  red
# grpm_survey_neuro  green
# 
# e in scala logaritimica
# 
# esempio di confronto: i dot verdi (neuro) che spiccano con valori alti stanno a significare che geni sono più interessanti dal punto di vista neurologico che nutrizionale

# ## misc

# #### Alternative "interest index" trial:
# 
#     ##GRPMX_report_new_int
#     if os.path.isfile('grpm_survey/GRPMX_report_new_int.csv'):
#         GRPMX_report_new_int = pd.read_csv('grpm_survey/GRPMX_report_new_int.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
#     else:
#         max_match_pmids = \
#             int(GRPMX_report['matching_pmids'].max()) #Pm(t)
#         max_pmids = \
#             int(GRPMX_report[GRPMX_report['matching_pmids']==max_match_pmids]['starting_pmid']) #P(t)
#         max_gene = \
#             (GRPMX_report[GRPMX_report['matching_pmids']==max_match_pmids]['gene'])
#         max_match_pmids, max_pmids, max_gene
# 
#     # duplicate df:
#     GRPMX_report_new_int = GRPMX_report
# 
#     upper = GRPMX_report_new_int['matching_pmids_ratio']#GRPMX_report_new_int['matching_pmids']/GRPMX_report_new_int['starting_pmid']
#     down = max_match_pmids/max_pmids
# 
#     GRPMX_report_new_int['new_interest_index'] = \
#         round((upper/down),3)
#     GRPMX_report_new_int.set_index('gene').T.to_csv('grpm_survey/GRPMX_report_new_int.csv')
# 
#     GRPMX_report_new_int[['reference_mesh', 'starting_pmidmesh', 'starting_pmid','starting_mesh','starting_rsid', 'matching_pmidmesh', 'matching_pmids', 'matching_mesh','matching_rsid', 'dropped_rsid', 'matching_rsid_pmid10','matching_rsid_pmid100']] = GRPMX_report_new_int[['reference_mesh', 'starting_pmidmesh', 'starting_pmid','starting_mesh','starting_rsid', 'matching_pmidmesh', 'matching_pmids', 'matching_mesh','matching_rsid', 'dropped_rsid','matching_rsid_pmid10','matching_rsid_pmid100']].astype(int)
#     GRPMX_report_new_int[['matching_mesh_ratio', 'matching_pmids_ratio','matching_pmidmesh_ratio', 'matching_rsid_ratio','matching_pmids_index','interest_index','new_interest_index']] = GRPMX_report_new_int[['matching_mesh_ratio', 'matching_pmids_ratio','matching_pmidmesh_ratio','matching_rsid_ratio','matching_pmids_index','interest_index','new_interest_index']].astype(float)
# 
#     GRPMX_report_new_int_sort = GRPMX_report_new_int.sort_values(by='interest_index',ascending=False)
#     GRPMX_report_new_int_sort[['gene','matching_pmids_ratio','interest_index','new_interest_index']]

# #new interest index trial
# #Matching PMIDs in Database
# 
#     sort_by = 'interest_index'
#     GRPMX_report_new_int_sort = GRPMX_report_new_int.sort_values(by=sort_by,ascending=False)
# 
#     n=100
#     y_axis = 'matching_pmids'
#     y_axis2 = 'new_interest_index'
#     x = GRPMX_report_new_int_sort.gene.iloc[100*n:100*(n+1)]
#     y = GRPMX_report_new_int_sort[y_axis].iloc[100*n:100*(n+1)]
#     y2 = (GRPMX_report_new_int_sort[y_axis2].iloc[100*n:100*(n+1)])*40
#     plt.figure(figsize=(4, 30))
#     plt.title('Matching PMIDs in Database\n sorted by '+sort_by, loc='center',pad=10)
# 
#     plt.barh(x,y, color = '#04ade3',label='Series 1')
#     plt.barh(x,y2, color = 'green',label='Series 2', align='edge')
#     plt.gca().invert_yaxis()
#     plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#     #plt.xlabel('pmid count', position=(0.5, 1.08))
#     plt.ylabel('genes')
#     plt.xlabel(y_axis, position=(0.5, 1.08))
#     ax = plt.gca()
#     ax.xaxis.set_label_position('top')
#     plt.show
#     #plt.savefig(r'Path')

# #### DEBUGGING
# 
#     #check existence
#     GRPM_report = pd.read_csv(db_name+'/GRPM_report.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
#     mask = GRPM_report.gene.str.contains('SLC2A2')
#     GRPM_report[mask]
# 

# #Check genes in DB
# 
#     genes = protein_coding_genes_list
#     geneindb = []
#     genenotindb = []
#     for i in genes:
#         if os.path.isfile(db_name+'/'+i+'_litvar2_variants4gene.csv'):
#             geneindb.append(i)
#             pass
#         else:
#             #print(i+' not preset')
#             genenotindb.append(i)
#             pass
#     len(geneindb)
#     geneindb = pd.Series(geneindb)
#     genenotindb = pd.Series(genenotindb)
#     #genedupp.to_csv('genedup.csv')
#     print(len(geneindb),len(genenotindb))
#     genenotindb.to_csv('genes_notindb.csv')
#     geneindb.to_csv('genes_indb.csv')

#     #compare db and survery
#     len(GRPM_report.gene), len(GRPMX_report.gene)
#     mask = GRPM_report['gene'].isin(GRPMX_report['gene'])
#     missing_ser = GRPM_report.gene[-mask]
# 
#     missing_list = missing_ser.tolist()
#     pyperclip.copy(str(missing_list))
#     missing_list

# # Visualize GRPM DB report
#     print('GRPM_report table:')
#     display(GRPM_report.T)
#     GRPM_report.columns
#     #pyperclip.copy(str(GRPM_report.columns))

# ## Report Analysis

# In[ ]:


# Matching PMIDs in Database
sort_by = 'pubmed_pmid'
GRPM_report_sort = GRPM_report.sort_values(by= sort_by,ascending=False)
GRPM_report_sort[['gene','pubmed_pmid']]


# In[ ]:


print('Sorted by: ', sort_by)
y_axis = 'pubmed_pmid'
x = GRPM_report_sort.gene.iloc[:40]
y = GRPM_report_sort[y_axis].iloc[:40]
plt.figure(figsize=(4, (len(x)*0.25)))
plt.title('PMIDs in Database\n sorted by '+sort_by, loc='center',pad=10)

plt.barh(x,y, color='#0483ff')
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
plt.show()
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")


# In[ ]:





# In[ ]:


#LOG GET-----------------------------------------------------
genes_len = len(pd.read_csv('human_genes_repo/H_GENES_proteincoding_genes.csv'))
genes_found = len(GRPM_report)
print('\nLOG of base grpm database building:','\ngene query survey_l:',genes_len)
print('genes_found:',genes_found, round(genes_found/genes_len, 2))
print('no results on:')
litvar= 3481 #taken form Console LOG notes
litvar2= 346 #taken form Console LOG notes
litvar1= 3135 #taken form Console LOG notes
nbib= 236 #taken form Console LOG notes
print('litvar:',litvar, round(litvar/genes_len, 2),
      '\n   litvar2:',litvar2, round(litvar2/genes_len, 2),
      '\n   litvar1:',litvar1, round(litvar1/genes_len, 2),
      '\nnbib:',nbib, round(nbib/genes_len, 2))


# In[ ]:


master_df


# In[ ]:


# Visualize GRPMX_report
n = 0
grpmx_report_00 = repo_int_sort_list[n]
print(type(grpmx_report_00['matching_mesh_ratio'][0]))
print('GRPMX_report table:')
#pyperclip.copy(str(GRPMX_report.columns))
len(grpmx_report_00.gene)
grpmx_report_00


# In[ ]:


# Sort Short GRPMX_report
sort = 'interest_index'
print('GRPMX_report table ''sorted by '+sort+':')
grpmx_report_00_sort = grpmx_report_00.sort_values(by= sort,ascending=False).reset_index(drop=True)

grpmx_report_00[['gene', 'matching_pmidmesh', 'matching_pmids','matching_mesh', 'matching_rsid', 'dropped_rsid', 'matching_mesh_ratio','matching_pmids_ratio', 'matching_pmidmesh_ratio','matching_rsid_ratio', 'matching_rsid_pmid10', 'matching_rsid_pmid100','matching_pmids_index','interest_index','matching_top10mesh', 'matching_top10rsid']]#.to_clipboard() #to excel


# ## Defining Database most studied genes
# 

# In[ ]:


print('Full DB report')
GRPM_report[['gene','pubmed_pmid']].sort_values(by='pubmed_pmid', ascending = False).reset_index(drop=True)#.columns


# In[ ]:


GRPM_report_sort


# In[ ]:


# create Bar Diagram
sort_by = 'pubmed_pmid'
GRPM_report_sort = GRPM_report.sort_values(by= sort_by,ascending=False)

gene_f = 0
gene_l = len(GRPM_report_sort.gene)
if gene_l > 30:
    gene_l = 30
y_axis = 'pubmed_pmid'
x = GRPM_report_sort.gene.iloc[gene_f:gene_l]
y = GRPM_report_sort[y_axis].iloc[gene_f:gene_l]
plt.figure(figsize=(4, (gene_l-gene_f)*0.25))
plt.title('Gene associated PMIDs in Database\n sorted by '+sort_by, loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


# DENSITY PLOT
sort = 'pubmed_pmid'
data = GRPM_report_sort[sort].iloc[:100]

sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)
sns.rugplot(np.array(data), color='r')

#plt.yscale('log')
plt.title('pubmed pmid density')
#plt.yscale('log')
plt.show()


# ## Defining pmid number treshold (GRPM DB)

# In[ ]:


# Definining pmid number treshold (GRPM DB)-------------------------------
thr_pmid_1 = 7
#thr_pmid_2 = 150
upper = 8000

thr_pmid_3 = 1
upper_2 = thr_pmid_1-1

threshold_GRPM_report_sort_1 = GRPM_report_sort[(GRPM_report_sort['pubmed_pmid']>=thr_pmid_1) & (GRPM_report_sort['pubmed_pmid']<=upper)].reset_index(drop=True)

#threshold_GRPM_report_sort_2 = GRPM_report_sort[(GRPM_report_sort['pubmed_pmid']>=threshold_2) & (GRPM_report_sort['pubmed_pmid']<=upper)].reset_index(drop=True)

threshold_GRPM_report_sort_3 = GRPM_report_sort[(GRPM_report_sort['pubmed_pmid']>=thr_pmid_3) & (GRPM_report_sort['pubmed_pmid']<=upper_2)].reset_index(drop=True)


#pyperclip.copy(str(threshold)+' '+ str(len(threshold_GRPM_report_sort)))
print('threshold by pmid number:')
print(len(threshold_GRPM_report_sort_3),'genes has <',thr_pmid_1, 'pmid')
print(len(threshold_GRPM_report_sort_1),' genes has >',thr_pmid_1, 'pmid')
threshold_GRPM_report_sort_1


# In[ ]:


pyperclip.copy(str(threshold_GRPM_report_sort_1.gene.sample(200).to_list()))


# In[ ]:


# DENSITY PLOT THRESHOLD
data = threshold_GRPM_report_sort_1.pubmed_pmid.iloc[:100]

sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)
sns.rugplot(np.array(data), color='r')

#plt.yscale('log')
plt.title('pubmed pmid')
#plt.yscale('log')
plt.show()


# # Visualize GRPMX reports

# ## Defining interesting genes
# #### Comparison: ratio Vs Interest index

#     ['YBX1', 'FTO', 'GAL', 'SOX9', 'ADD1']
#     Index(['gene', 'reference_mesh', 'starting_pmidmesh', 'starting_pmid',
#            'starting_mesh', 'starting_rsid', 'matching_pmidmesh', 'matching_pmids',
#            'matching_mesh', 'matching_rsid', 'dropped_rsid', 'matching_mesh_%',
#            'matching_pmids_%', 'matching_pmidmesh_%', 'matching_rsid_%',
#            'matching_rsid_pmid10', 'matching_rsid_pmid100', 'matching_top10mesh',
#            'matching_top10rsid', 'total_runtime'],
#           dtype='object', name='Unnamed: 0')

# In[ ]:


master_df


# In[ ]:


# Comparison: matching_pmids_ratio Vs interest_index

#Set report to analyze
dir1 = 0
GRPMX_report_int = repo_int_sort_list[dir1]

# create Bar Diagram
sort_by = 'matching_pmids_ratio'
GRPMX_report_int_sort = GRPMX_report_int.sort_values(by= sort_by,ascending=False)
gene_f = 0
gene_l =   len(GRPMX_report_int_sort)
if gene_l  > 40:
    gene_l = 50
y_axis = 'interest_index'
x = GRPMX_report_int_sort.gene.iloc[:gene_l]
y = GRPMX_report_int_sort[y_axis].drop_duplicates().iloc[:gene_l]
plt.figure(figsize=(4, gene_l*0.25))
plt.title('Gene '+y_axis+ ':\n sorted by '+ sort_by, loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


# create Bar Diagram
sort_by = 'matching_rsid_ratio'
GRPMX_report_int_sort = GRPMX_report_int.sort_values(by= sort_by,ascending=False)
print('Sorted by: ', sort_by)
gene_f = 0
gene_l =   len(GRPMX_report_int_sort)
if gene_l  > 30:
    gene_l = 50
y_axis = 'matching_rsid'
x = GRPMX_report_int_sort.gene.iloc[:gene_l]
y = GRPMX_report_int_sort[y_axis].drop_duplicates().iloc[:gene_l]
plt.figure(figsize=(4, gene_l*0.25))
plt.title('Matching rsid in Database\n sorted by '+sort_by, loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
plt.show()
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")


# In[ ]:


# Matching PMIDs in Database
sort_by = 'matching_pmids_index'
GRPMX_report_int_sort = GRPMX_report_int.sort_values(by= sort_by,ascending=False)
print('Sorted by: ', sort_by)
print('in ', directory)

gene_f , gene_l = 0, 0
gene_l =   len(GRPMX_report_int_sort)
if gene_l  > 30:
    gene_l = 50
y_axis = 'matching_pmids_index'
x = GRPMX_report_int_sort.gene.iloc[:gene_l]
y = GRPMX_report_int_sort[y_axis].drop_duplicates().iloc[:gene_l]
plt.figure(figsize=(4, gene_l*0.25))
plt.title('Matching PMIDs in Database\n sorted by '+ sort_by, loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


# Matching PMIDs in Database
sort_by = 'interest_index'
GRPMX_report_int_sort = GRPMX_report_int.sort_values(by= sort_by,ascending=False)
print('Sorted by: ', sort_by)
print('in ', directory)
gene_f , gene_l, n = 0, 40, 1
gene_l =   len(GRPMX_report_int)
if gene_l  > 30:
    gene_l = 40
y_axis = 'interest_index'
x = GRPMX_report_int_sort.gene.iloc[gene_f*n:gene_l*n]
y = GRPMX_report_int_sort[y_axis].iloc[gene_f*n:gene_l*n]
plt.figure(figsize=(4, (gene_l-gene_f) * 0.25))
plt.title('Matching PMIDs in Database\n sorted by '+sort_by+'\n in'+directory, loc='center',pad=10)

plt.barh(x,y, color = '#04ade3')
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


# Matching PMIDs in Database

#Set report to analyze
n = 1
GRPMX_report_int_2 = repo_int_sort_list[n]

sort_by = 'interest_index'
GRPMX_report_int_2_sort = GRPMX_report_int_2.sort_values(by= sort_by,ascending=False)
print('Sorted by: ', sort_by)
print('in ', directories[n])
gene_f , gene_l, n = 0, 40, 1
gene_l =   len(GRPMX_report_int)
if gene_l  > 30:
    gene_l = 40
y_axis = 'interest_index'
x = GRPMX_report_int_sort.gene.iloc[gene_f*n:gene_l*n]
y = GRPMX_report_int_sort[y_axis].iloc[gene_f*n:gene_l*n]
plt.figure(figsize=(4, (gene_l-gene_f)*0.25))
plt.title('Matching PMIDs in Database\n sorted by '+sort_by+'\n in'+directories[n], loc='center',pad=10)

plt.barh(x,y, color = '#04ade3')
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# ## Matching PMIDs in Database

#     sort_by = 'new_matching_pmids_index'
#     GRPMX_report_new_int_sort = GRPMX_report_new_int.sort_values(by= sort_by,ascending=False)
#     print('Sorted by: ', sort_by)
# 
#     y_axis = 'new_matching_pmids_index'
#     x = GRPMX_report_new_int_sort.gene.iloc[:40]
#     y = GRPMX_report_new_int_sort[y_axis].iloc[:40]
#     plt.figure(figsize=(4, 8))
#     plt.title('Matching PMIDs in Database\n sorted by '+sort_by, loc='center',pad=10)
# 
#     plt.barh(x,y)
#     plt.gca().invert_yaxis()
#     plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#     #plt.xlabel('pmid count', position=(0.5, 1.08))
#     plt.ylabel('genes')
#     plt.xlabel(y_axis, position=(0.5, 1.08))
#     ax = plt.gca()
#     ax.xaxis.set_label_position('top')
#     #plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
#     plt.show()

# # Analyze GRPMX Surveys

# ## Cooccurance matrix module:

# In[ ]:


#for i in range(len(directories)):
#    master_df.at[i, 'grpm_mesh'] = grpm_out_list[i].mesh.nunique()
master_df
grpm_out_list[0].columns


# In[ ]:





# In[ ]:


# COOCCURANCE MATRIX MODULE------------
# second matrix build check

# choose value
value = 'gene'
# ['gene', 'rsid', 'pmids', 'mesh']

# Assuming grpm_out_list is a list of objects and i is the index
lists = []
for i in range(len(grpm_out_list)):
    lists.append(set(grpm_out_list[i][value].to_list()))

# set the lists
num_lists = len(lists)

# Initialize a 2D list of zeros with dimensions equal to the number of lists
cooccur_matrix = [[0] * num_lists for i in range(num_lists)]
type(cooccur_matrix)
print(len(set(lists[1]) & set(lists[1])))
print(len(set(lists[1]) & set(lists[3])))

# Loop over all pairs of lists and count the number of co-occurring elements
for i in range(num_lists):
    for j in range(num_lists):
        if i == j:
            cooccur_matrix[i][j] = len(lists[i])
        else:
            cooccur_matrix[i][j] = len(set(lists[i]) & set(lists[j]))

# Convert the 2D list to a Pandas DataFrame
cooccur_df = pd.DataFrame(cooccur_matrix,
                          columns=[i for i in directories],
                          index=[i for i in directories])
                          #columns=['grpm{}'.format(i+1) for i in range(num_lists)],
                          #index=['grpm{}'.format(i+1) for i in range(num_lists)])

# Print the resulting DataFrame
cooccur_df.to_clipboard()
cooccur_df.columns
print('Cooccurrance Matrix: '+value)
sub_cooccur_df = cooccur_df[['grpm_survey_pcg_nutri', 'grpm_survey_pcg_ob_bmi',
                             'grpm_survey_pcg_dmt2_ms', 'grpm_survey_pcg_cvd',
                             'grpm_survey_pcg_vitam', 'grpm_survey_pcg_eat_taste',
                             'grpm_survey_pcg_intol', 'grpm_survey_pcg_aller',
                             'grpm_survey_pcg_oxi_stress', 'grpm_survey_pcg_xeno']][:survey_l]
sub_cooccur_df


# In[ ]:


# Row max values------------------
plt.clf()
fig, ax = plt.subplots()

# Iterate over each row and plot individually
for i, row in enumerate(sub_cooccur_df[:survey_l].values):
    row_values = row.astype(float)  # Convert string values to floats
    max_value = np.max(row_values)  # Calculate the maximum value for the current row
    #heatmap = ax.imshow([row_values], cmap='YlGnBu', vmax=max_value, aspect='auto', extent=[-0.5, len(row)-0.5, i-0.5, i+0.5])
    heatmap = ax.imshow([row_values], cmap='YlGnBu', vmax=max_value, aspect='auto', extent=[+0.5, len(row)+0.5, i+0.5, i-0.5])
    #ax.text(len(row)-0.1, i, f"{max_value:.2f}", ha='left', va='center')
    ax.text(len(row)-0.1, i, f"{int(max_value)}", ha='left', va='center')

# Add colorbar
cbar = plt.colorbar(heatmap)

# Set the tick labels and axis labels
ax.set_xticks(np.arange(sub_cooccur_df[:survey_l].shape[1]))
ax.set_yticks(np.arange(sub_cooccur_df[:survey_l].shape[0]))
ax.set_xticklabels(list(master_df.directories[:survey_l]), rotation=45, ha='right')
ax.set_yticklabels(list(master_df.directories[:survey_l]))

# Set the axis labels
#ax.set_xlabel('Columns')
#ax.set_ylabel('Rows')

# Set the title
ax.set_title('Occurrence Matrix: ' + value)

# Show the plot
plt


# 
# but uhm uhm, at which point the scale must be updated accordingly
# 
# I mean, you scale the values and you end up with a range [0,1]
# 
# scaling all the rows to their maximum by having [0,1]values for each row would also simplify the code for the heatmap although we would leave behind the differences in the numeric data. I can do this
# 
# what is synoptic must be comparable.
# 

# In[ ]:


# modify coocur matrix with range [0,1]----
sub_cooccur_df.to_clipboard()
sub_cooccur_df_normalized = sub_cooccur_df.div(sub_cooccur_df.max(axis=1), axis=0) #div()= divide

sub_cooccur_df_normalized.to_clipboard()

#Overral values----------------------------
plt.clf()
#Create a figure and axis
fig, ax = plt.subplots()

#Create the heatmap
heatmap = ax.imshow(sub_cooccur_df_normalized[:survey_l], cmap='YlGnBu')

#Add colorbar
cbar = plt.colorbar(heatmap)

#Set the tick labels and axis labels
ax.set_xticks(np.arange(sub_cooccur_df_normalized[:survey_l].shape[1]))
ax.set_yticks(np.arange(sub_cooccur_df_normalized[:survey_l].shape[0]))
ax.set_xticklabels(list(master_df.directories[:survey_l]), rotation=45, ha='right')
ax.set_yticklabels(list(master_df.directories[:survey_l]))

#Set the axis labels
#ax.set_xlabel('Columns')
#ax.set_ylabel('Rows')

#Set the title
ax.set_title('Occurrence Matrix: '+value+' (normalized)')

#Show the plot
plt.show()


# In[ ]:


tag= 'oxi_stress'
ref = pd.read_csv('ref-mesh-archive/ref_mesh_'+tag+'.csv')
ref['Preferred Label'].drop_duplicates()


# In[ ]:


# Function to highlight cells
def highlight_cell(x):
    color = 'background-color: yellow'
    return color

# Apply styling to the DataFrame
styled_df = cooccur_df.style.applymap(highlight_cell)

# Display the styled DataFrame
styled_df


# ## Apply ‘interest index’ threshold

# In[ ]:


master_df


# #### (1) thresholding all surveys report

# In[ ]:


# def thresholding finction
threshold_f = 0.0055

def thresholding(dataframe):
    output = dataframe[(dataframe['interest_index']>=threshold_f) & (dataframe['interest_index']<=1)].reset_index(drop=True)
    return output

repo_thre_list = []
for repo, label in zip(repo_int_sort_list,master_df.directories):
    dataframe = thresholding(repo)
    #globals()[name] = dataframe
    repo_thre_list.append(dataframe)
    #print(label,'df thresholded')

repo_thre_list[0]

for i in range(len(comp_df_trial_all)):
    repo = repo_thre_list[i]
    comp_df_trial_all.loc[i, 'genes_thr0.0055'] = repo.gene.nunique()

#comp_df_trial_all


# grpm_out_list[0]

# #### (2) add interest index to grpm out list genes and sort it

# In[ ]:


grpm_out_int_sort = []
for grpm, i  in zip(grpm_out_list, range(len(grpm_out_list))):
    ## grpm [gene, interstr index] to merge wit grpmx output
    #print('add interest index to grpmx_output', directories[dir1])
    small_dummy = repo_int_sort_list[i][['gene','interest_index']]

    #grpmx_out_merge = pd.merge(grpmx_out_1, small_dummy, left_on='gene', right_on='gene')
    grpmx_out_merge_sort = pd.merge(grpm, small_dummy, left_on='gene', right_on='gene') \
        .sort_values(by=['interest_index','rsid'], ascending =False).reset_index(drop=True)
    grpm_out_int_sort.append(grpmx_out_merge_sort)

#grpm_out_int_sort[3]


# In[ ]:


grpm_out_int_sort[0][['gene','interest_index']].drop_duplicates()


# grpm_out_int_sort[3][['gene','interest_index']].drop_duplicates()

# #### (3) filtering grpmx_output_full with selected interesting genes

# In[ ]:


grpm_out_thre = []
for i in range(len(grpm_out_int_sort)):
    output = grpm_out_int_sort[i][grpm_out_int_sort[i].gene.isin(repo_thre_list[i].gene)]
    grpm_out_thre.append(output)

grpm_out_thre[2]['gene','interest_index'].drop_duplicates()


# #### (4) coocur matrix and visualise heatmap

# In[ ]:


# COOCCURANCE MATRIX MODULE------------
# second matrix build check

# choose value
value = 'gene'
# ['gene', 'rsid', 'pmids', 'mesh']

# Assuming grpm_out_list is a list of objects and i is the index
lists = []
for i in range(len(grpm_out_thre)):
    lists.append(set(grpm_out_thre[i][value].to_list()))

# set the lists
num_lists = len(lists)

# Initialize a 2D list of zeros with dimensions equal to the number of lists
cooccur_2_matrix = [[0] * num_lists for i in range(num_lists)]
type(cooccur_2_matrix)
print(len(set(lists[1]) & set(lists[1])))
print(len(set(lists[1]) & set(lists[3])))

# Loop over all pairs of lists and count the number of co-occurring elements
for i in range(num_lists):
    for j in range(num_lists):
        if i == j:
            cooccur_2_matrix[i][j] = len(lists[i])
        else:
            cooccur_2_matrix[i][j] = len(set(lists[i]) & set(lists[j]))

# Convert the 2D list to a Pandas DataFrame
cooccur_2_df = pd.DataFrame(cooccur_2_matrix,
                          columns=[i for i in directories],
                          index=[i for i in directories])
#columns=['grpm{}'.format(i+1) for i in range(num_lists)],
#index=['grpm{}'.format(i+1) for i in range(num_lists)])

# Print the resulting DataFrame
cooccur_2_df.to_clipboard()
cooccur_2_df.columns
print('Cooccurrance Matrix: '+ value)
sub_cooccur_2_df = cooccur_2_df[['grpm_survey_pcg_nutri', 'grpm_survey_pcg_ob_bmi',
                             'grpm_survey_pcg_dmt2_ms', 'grpm_survey_pcg_cvd',
                             'grpm_survey_pcg_vitam', 'grpm_survey_pcg_eat_taste',
                             'grpm_survey_pcg_intol', 'grpm_survey_pcg_aller',
                             'grpm_survey_pcg_oxi_stress', 'grpm_survey_pcg_xeno']][:survey_l]
sub_cooccur_2_df


# In[ ]:


# Row max values------------------
plt.clf()
fig, ax = plt.subplots()

# Iterate over each row and plot individually
for i, row in enumerate(sub_cooccur_2_df[:survey_l].values):
    row_values = row.astype(float)  # Convert string values to floats
    max_value = np.max(row_values)  # Calculate the maximum value for the current row
    #heatmap = ax.imshow([row_values], cmap='YlGnBu', vmax=max_value, aspect='auto', extent=[-0.5, len(row)-0.5, i-0.5, i+0.5])
    heatmap = ax.imshow([row_values], cmap='YlGnBu', vmax=max_value, aspect='auto', extent=[+0.5, len(row)+0.5, i+0.5, i-0.5])
    #ax.text(len(row)-0.1, i, f"{max_value:.2f}", ha='left', va='center')
    ax.text(len(row)-0.1, i, f"{int(max_value)}", ha='left', va='center')

# Add colorbar
cbar = plt.colorbar(heatmap)

# Set the tick labels and axis labels
ax.set_xticks(np.arange(sub_cooccur_2_df[:survey_l].shape[1]))
ax.set_yticks(np.arange(sub_cooccur_2_df[:survey_l].shape[0]))
ax.set_xticklabels(list(master_df.directories[:survey_l]), rotation=45, ha='right')
ax.set_yticklabels(list(master_df.directories[:survey_l]))

# Set the axis labels
#ax.set_xlabel('Columns')
#ax.set_ylabel('Rows')

# Set the title
ax.set_title('Occurrence Matrix: ' + value+'\ninterest_inxed threshold: 0.0055')

# Show the plot
plt


# In[ ]:


# modify coocur matrix with range [0,1]----
sub_cooccur_2_df.to_clipboard()
sub_cooccur_2_df_normalized = sub_cooccur_2_df.div(sub_cooccur_2_df.max(axis=1), axis=0) #div()= divide

sub_cooccur_2_df_normalized.to_clipboard()

#Overral values----------------------------
plt.clf()
#Create a figure and axis
fig, ax = plt.subplots()

#Create the heatmap
heatmap = ax.imshow(sub_cooccur_2_df_normalized[:survey_l], cmap='YlGnBu')

#Add colorbar
cbar = plt.colorbar(heatmap)

#Set the tick labels and axis labels
ax.set_xticks(np.arange(sub_cooccur_2_df_normalized[:survey_l].shape[1]))
ax.set_yticks(np.arange(sub_cooccur_2_df_normalized[:survey_l].shape[0]))
ax.set_xticklabels(list(master_df.directories[:survey_l]), rotation=45, ha='right')
ax.set_yticklabels(list(master_df.directories[:survey_l]))

#Set the axis labels
#ax.set_xlabel('Columns')
#ax.set_ylabel('Rows')

#Set the title
ax.set_title('Occurrence Matrix: '+value+' (normalized)\ninterest_inxed threshold: 0.0055')

#Show the plot
plt.show()


# ## Definining ‘interest index’ threshold process

# In[ ]:


# Definining ‘interest index’ threshold----------------------------------
threshold_f1 = 0.01
threshold_f2 = 0.0055

dir1 = 8 # choose dir
dir2 = 9 # choose dir
repo_int_sort_1 = repo_int_sort_list[dir1]
repo_int_sort_2 = repo_int_sort_list[dir2]

threshold_GRPMX_report_int_sort = repo_int_sort_1[(repo_int_sort_1['interest_index']>=threshold_f1) & (repo_int_sort_1['interest_index']<=1)].reset_index(drop=True)

print('thresholding',directories[dir1])
threshold_GRPMX_report_int_sort2 = repo_int_sort_1[(repo_int_sort_1['interest_index']>=threshold_f2) & (repo_int_sort_1['interest_index']<=1)].reset_index(drop=True)
print('thresholding',directories[dir2])
threshold_GRPMX_report_int_2_sort2 = repo_int_sort_2[(repo_int_sort_2['interest_index']>=threshold_f2) & (repo_int_sort_2['interest_index']<=1)].reset_index(drop=True)

print('threshold by \'interest index\':','\n  ', directories[dir1],'\n  ', directories[dir2])
#pyperclip.copy(str(threshold)+' '+ str(len(threshold_GRPMX_report_int_sort)))
threshold_GRPMX_report_int_sort2#.to_csv('threshold_GRPMX_report_int_sort.csv')#2#.to_clipboard()


# #### filtering grpmx_output_full with selected interesting genes

# In[ ]:


master_df


# In[ ]:


#grpmx_out_1 = pd.read_csv(directories[dir1]+'/grpmx_filtered_output.csv',index_col=0)
#grpmx_out_2 = pd.read_csv(directories[dir2]+'/grpmx_filtered_output.csv',index_col=0)

# the same as grpm_out_list
grpmx_out_1 = grpm_out_list[dir1]
grpmx_out_2 = grpm_out_list[dir2]


# In[ ]:


grpmx_out_1


# In[ ]:


## grpm [gene, interstr index] to merge wit grpmx output
print('add interest index to grpmx_output', directories[dir1])
small_dummy = repo_int_sort_list[dir1][['gene','interest_index']]

#grpmx_out_merge = pd.merge(grpmx_out_1, small_dummy, left_on='gene', right_on='gene')
grpmx_out_merge_sort = pd.merge(grpmx_out_1, small_dummy, left_on='gene', right_on='gene')\
    .sort_values(by=['interest_index','rsid'], ascending =False).reset_index(drop=True)
grpmx_out_merge_sort


# In[ ]:


## grpm [gene, interstr index] to merge wit grpmx output
print('add interest index to grpmx_output',directories[dir2])
small_dummy = repo_int_sort_list[dir2][['gene','interest_index']]
#grpmx_out_2_merge = pd.merge(grpmx_out_2, small_dummy, left_on='gene', right_on='gene')
grpmx_out_2_merge_sort = pd.merge(grpmx_out_2, small_dummy, left_on='gene', right_on='gene')\
    .sort_values(by=['interest_index','rsid'], ascending =False).reset_index(drop=True)
grpmx_out_2_merge_sort


# In[ ]:


#Appling for thershold genes-------------------------------------
print('Appling "interest index" threshold',threshold_f2,'on', directories[dir1])
grpmx_out_thr = grpmx_out_merge_sort[grpmx_out_merge_sort.gene.isin(threshold_GRPMX_report_int_sort2.gene)]
print('Appling "interest index" threshold',threshold_f2,'on', directories[dir2])
grpmx_out_thr_2 = grpmx_out_2_merge_sort[grpmx_out_2_merge_sort.gene.isin(threshold_GRPMX_report_int_2_sort2.gene)]

print('')
print(directories[dir1], 'threshold',threshold_f2,'applied')
print(grpmx_out_thr.gene.nunique(), 'genes')
print(grpmx_out_thr.rsid.nunique(), 'rsid with pmid with mesh')
print(grpmx_out_thr.mesh.nunique(), 'mesh')
print(round(grpmx_out_thr.rsid.nunique()/grpmx_out_thr.gene.nunique(), 0), 'risd for gene (mean)')
print ('Based on', repo_int_sort_list[dir1].reference_mesh[0], 'reference mesh terms')

print('')
print(directories[dir2], 'threshold',threshold_f2,'applied')
print(grpmx_out_thr_2.gene.nunique(), 'genes')
print(grpmx_out_thr_2.rsid.nunique(), 'rsid with pmid with mesh')
print(grpmx_out_thr_2.mesh.nunique(), 'mesh')
print(round(grpmx_out_thr_2.rsid.nunique()/grpmx_out_thr_2.gene.nunique(), 0), 'risd for gene (mean)')
print ('Based on', repo_int_sort_list[dir2].reference_mesh[0], 'reference mesh terms')


# In[ ]:


grpmx_out_thr_ngene = grpmx_out_thr.gene.nunique()
grpmx_out_thr_nrsid = grpmx_out_thr.rsid.nunique()
grpmx_out_thr_npmid = grpmx_out_thr.pmids.nunique()
grpmx_out_thr_nmesh = grpmx_out_thr.mesh.nunique()
print(directories[dir1])
print(grpmx_out_thr_ngene, 'genes')
print(grpmx_out_thr_nrsid, 'rsids')
print(grpmx_out_thr_npmid, 'pmids')
print(grpmx_out_thr_nmesh, 'meshs')
print('')
grpmx_out_thr_2_ngene = grpmx_out_thr_2.gene.nunique()
grpmx_out_thr_2_nrsid = grpmx_out_thr_2.rsid.nunique()
grpmx_out_thr_2_npmid = grpmx_out_thr_2.pmids.nunique()
grpmx_out_thr_2_nmesh = grpmx_out_thr_2.mesh.nunique()
print(directories[dir2])
print(grpmx_out_thr_2_ngene, 'genes')
print(grpmx_out_thr_2_nrsid, 'rsids')
print(grpmx_out_thr_2_npmid, 'pmids')
print(grpmx_out_thr_2_nmesh, 'meshs')


# In[ ]:


# Analyze intersections
directory_2 = directories[dir2]
directory   = directories[dir1]

grpmx_thr_genes = grpmx_out_thr.gene.drop_duplicates()
grpmx_thr_rsids = grpmx_out_thr.rsid.drop_duplicates()
grpmx_thr_pmids = grpmx_out_thr.pmids.drop_duplicates()
grpmx_thr_meshs = grpmx_out_thr.mesh.drop_duplicates()
grpmx_thr_2_genes = grpmx_out_thr_2.gene.drop_duplicates()
grpmx_thr_2_rsids = grpmx_out_thr_2.rsid.drop_duplicates()
grpmx_thr_2_pmids = grpmx_out_thr_2.pmids.drop_duplicates()
grpmx_thr_2_meshs = grpmx_out_thr_2.mesh.drop_duplicates()

mask_gene_1 = grpmx_out_thr_2.gene.isin( grpmx_thr_genes)
mask_rsid_1 = grpmx_out_thr_2.rsid.isin( grpmx_thr_rsids)
mask_pmid_1 = grpmx_out_thr_2.pmids.isin(grpmx_thr_pmids)
mask_mesh_1 = grpmx_out_thr_2.mesh.isin( grpmx_thr_meshs)
mask_gene_2 = grpmx_out_thr.gene.isin(grpmx_thr_2_genes)
mask_rsid_2 = grpmx_out_thr.rsid.isin(grpmx_thr_2_rsids)
mask_pmid_2 = grpmx_out_thr.pmids.isin(grpmx_thr_2_pmids)
mask_mesh_2 = grpmx_out_thr.mesh.isin(grpmx_thr_2_meshs)

grpmx_out_thr_maskgene = grpmx_out_thr[mask_gene_2]
grpmx_out_thr_maskrsid = grpmx_out_thr[mask_rsid_2]
grpmx_out_thr_maskpmid = grpmx_out_thr[mask_pmid_2]
grpmx_out_thr_maskmesh = grpmx_out_thr[mask_mesh_2]
grpmx_out_thr_2_maskgene = grpmx_out_thr_2[mask_gene_1]
grpmx_out_thr_2_maskrsid = grpmx_out_thr_2[mask_rsid_1]
grpmx_out_thr_2_maskpmid = grpmx_out_thr_2[mask_pmid_1]
grpmx_out_thr_2_maskmesh = grpmx_out_thr_2[mask_mesh_1]

print(directory,'vs',directory_2)
print(grpmx_out_thr_maskgene.gene.nunique() ,  'genes intersection on', grpmx_out_thr_ngene)
print(grpmx_out_thr_maskrsid.rsid.nunique() ,  'rsid intersection on',  grpmx_out_thr_nrsid)
print(grpmx_out_thr_maskpmid.pmids.nunique(),  'pmid intersection on',  grpmx_out_thr_npmid)
print(grpmx_out_thr_maskmesh.mesh.nunique() ,  'mesh intersection on',  grpmx_out_thr_nmesh)
print('')
print(directory_2,'vs',directory)
print(grpmx_out_thr_2_maskgene.gene.nunique(), 'genes intersection on', grpmx_out_thr_2_ngene)
print(grpmx_out_thr_2_maskrsid.rsid.nunique(), 'rsid intersection on',  grpmx_out_thr_2_nrsid)
print(grpmx_out_thr_2_maskpmid.pmids.nunique(),'pmid intersection on',  grpmx_out_thr_2_npmid)
print(grpmx_out_thr_2_maskmesh.mesh.nunique(), 'mesh intersection on',  grpmx_out_thr_2_nmesh)
print('')
print('Overlapping index:')
print(round(grpmx_out_thr_maskgene.gene.nunique() /grpmx_out_thr_ngene*100,2),'%','genes overlap')
print(round(grpmx_out_thr_maskrsid.rsid.nunique() /grpmx_out_thr_nrsid*100,2),'%', 'rsid overlap')
print(round(grpmx_out_thr_maskpmid.pmids.nunique()/grpmx_out_thr_npmid*100,2),'%', 'pmid overlap')
print(round(grpmx_out_thr_maskmesh.mesh.nunique() /grpmx_out_thr_nmesh*100,2),'%', 'mesh overlap')


# ## analyzing mesh on selected gene data

# In[ ]:


# just top10mesh from grpm report

type(threshold_GRPMX_report_int_sort2['matching_top10mesh'][0])
import ast
# Define a function to convert a string to a list
def str_to_list(s):
    return ast.literal_eval(s)

mesh_exp = []
for i in range(len(threshold_GRPMX_report_int_sort2)):
    #for mesh in threshold_GRPMX_report_int_sort2['matching_top10mesh'].apply(str_to_list)[i]:
    for mesh in threshold_GRPMX_report_int_sort2['matching_top10mesh'].apply(ast.literal_eval)[i]:
        f = threshold_GRPMX_report_int_sort2['gene'][i], mesh
        mesh_exp.append(f)

gene_mesh = pd.DataFrame(mesh_exp)


# In[ ]:


gene_mesh.columns = ['gene','mesh']
print('gene:', gene_mesh.gene.nunique())
print('mesh:', gene_mesh.mesh.nunique())


# In[ ]:


contingency_table = pd.crosstab(gene_mesh['gene'], gene_mesh['mesh'])
contingency_table


# In[ ]:


# groupby describe
gene_mesh_count = gene_mesh.groupby('mesh').describe().reset_index()
gene_mesh_count.columns = [('mesh_'),
                           ('gene_count'),
                           ('gene_unique'),
                           ('gene_top'),
                           ('gene_freq')]
gene_mesh_count.sort_values(by='gene_count', ascending =False).reset_index(drop=True)


# In[ ]:


#get pmid number range for threshold_GRPMX_report_int_sort2
subset = threshold_GRPMX_report_int_sort2[['gene','starting_pmid','matching_pmids']]
subset_sort = subset.sort_values(by='starting_pmid', ascending = False).reset_index(drop=True)
print(len(threshold_GRPMX_report_int_sort2), 'filtered genes with threshold',threshold_f2 )
print('    min matching pmid value: ', subset_sort.matching_pmids.min())
print('    min starting pmid value: ', subset_sort.starting_pmid.min())

pmid_num_db = GRPM_report_sort[GRPM_report_sort['pubmed_pmid'] >= subset_sort.starting_pmid.min()].count()['pubmed_pmid']
print(pmid_num_db, 'genes in GRPM db has >=', subset_sort.starting_pmid.min())
print( len(threshold_GRPMX_report_int_sort2),'genes filterd','over', pmid_num_db, 'equivalent genes in db')
subset_sort.gene.to_clipboard(index=False, header=False)


# In[ ]:


# DENSITY PLOT THRESHOLD
data = threshold_GRPMX_report_int_sort.interest_index#.iloc[:100]

sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)
sns.rugplot(np.array(data), color='r')

#plt.yscale('log')
plt.title('interest_index')
#plt.yscale('log')
plt.show()


# # Positive and negative control
# Truth table contructor
# 

# In[ ]:


#Truth table contructor:

genesthr01 = threshold_GRPMX_report_int_sort.gene
genesthr02 = threshold_GRPMX_report_int_sort2.gene
#genesthr03 = threshold_GRPMX_report_int_sort3.gene
#genesnewint = GRPMX_report_new_int.sort_values(by='new_interest_index',ascending=False).gene.iloc[:879]

genestesi = pd.Series(['ADIPOQ','ADRB1','ADRB3','APOA2','APOA5','DRD2','FTO','GHRL','GHSR','IL-1B','LEP','MC4R','NPY','PPARG','SLC6A4','ACE','ADRB2','AGT','AGTR1','APOA1','APOB','APOC3','APOE','CETP','CYBA','EDN1','F5','FGB','GJA4','HMGCR','LIPC','LPA','LPL','MMP3','MTHFR','NOS3','PON1','PPARA','PPARGC1A','PROCR','SELE','SREBF2','VEGF','FABP2','IRS1','LEPR','MTNR1B','TCF7L2','ALPL','BCO1','CBS','FADS1','FUT2','GC','MTR','MTRR','NBPF3','TCN2','VDR','ADD1','CD36','CYP11B2','SLC2A2','TAS2R38','ADH1B','ADH1C','ALDH2','ALDOB','AOC1','CYP1A2','HLA-DQA1','HLA-DQB1','HNMT','MCM6','EPHX1','G6PD','GPX1','GSTP1','GSTT1','SOD2','SOD3','SOUX','CRP','IL-10','IL-6','IL6R','TNF'])
genespnn = pd.Series(['ACE','ADH1B','ADIPOQ','ADRB2','ADRB3','AGT','ALDH2','APOA1','APOA2','APOB','APOC3','APOE','CBS','CD36','CETP','CRP','CYP11B2','CYP1A2','FTO','GC','GHSR','LEP','LEPR','LIPC','LPL','MTHFR','MTR','MTRR','NOS3','NPY','PPARA','PPARG','SLC2A2','SREBF2','TAS2R38','TCF7L2','TNF','VDR','ALDH1A2','FADS2','IL6','LTA','DRD2','MC4R','ABCA4','ABCB4','ABCC4','ABCG2','ABCG5','ABCG8','ACHE','ACSL5','ACTN3','ADCY3','ADH4','ADORA2A','ADRA2B','AGPAT2','AGRP','ALDH1A3','APOA4','AR','ARC','ATF4','ATP2A1','ATP7A','ATP7B','BCDIN3D','BDNF','BGLAP','BHMT','BIRC7','BMP6','BSCL2','BTG2','BUD13','C1S','C2','C2orf69','C3','C4orf46','C5','C7','C9','CALU','CAV2','CCDC50','CCK','CD14','CD4','CDH2','CDK14','CFH','CFI','CGAS','CHDH','CHKA','CHKB','CLOCK','COMT','CP','CPVL','CRY2','CS','CTH','CUBN','CUL4B','CYP24A1','CYP26A1','CYP26B1','CYP26C1','CYP27A1','CYP27B1','CYP2A6','CYP2C19','CYP2C9','CYP2E1','CYP2R1','CYP3A4','CYP3A5','CYP7A1','DBP','DHCR7','DHFR','DHRS3','DIS3L','DNAAF4','DUOX2','EBP','EEF2K','ELOVL2','ELOVL5','ENTPD6','ETV5','EZH2','F2','F7','F9','FAAH','FAAH2','FAIM2','FASN','FASTK','FGF21','FGF23','GAD2','GAS6','GCG','GGCX','GIP','GIPR','GLRX','GLT8D2','GNAS','GNAT3','GNB3','GNPDA2','GRK4','GRM8','GRP','H19','HAMP','HJV','HP','HR','HSD11B2','IDS','IGLJ3','INSIG2','IRX5','ISX','KCNH7','KCNH8','KCTD15','KIF5B','KSR2','LINC00474','LMNA','LPO','LRAT','LRP2','LRP5','MAF','MAX','MB','MC2R','MC3R','MC5R','MDS2','MGP','MIPEP','MLXIPL','MT2A','MTHFS','MYLIP','NADSYN1','NAE1','NAT2','NCBP2','NDN','NEDD4','NEDD4L','NHS','NLRP3','NLRP6','NMB','NR1H3','OCM','OGT','OR6N2','PAH','PANK2','PARP2','PARP3','PC','PEMT','PER2','PGC','PGR','PIP5K1B','PLTP','PNLIP','PNMT','PNPLA3','POMC','PPARD','PRRG2','PRRG3','PRRG4','PSMD3','PTH','PWRN2','PYY','RAB21','RAPGEF3','RBP4','RBPMS','RDH12','RDH5','RDH8','RGS18','RPE','RPE65','RXRA','SCAP','SCD','SEC16B','SGK1','SH2B1','SHMT2','SIRT2','SIRT3','SIRT4','SIRT5','SIRT6','SIRT7','SLC10A2','SLC13A3','SLC13A4','SLC16A1','SLC16A9','SLC17A1','SLC17A3','SLC19A2','SLC19A3','SLC22A1','SLC22A11','SLC22A12','SLC22A3','SLC22A6','SLC22A8','SLC23A1','SLC23A2','SLC24A3','SLC25A32','SLC26A4','SLC26A7','SLC2A5','SLC2A8','SLC2A9','SLC30A1','SLC30A2','SLC30A8','SLC39A11','SLC39A4','SLC40A1','SLC44A1','SLC4A5','SLC52A2','SLC52A3','SLC5A1','SLC5A12','SLC5A8','SLC8A1','SNORA48','SRR','SSB','STAT3','STRA6','SVIL','TAS1R2','TAS1R3','TCF4','TES','TF','TFAM','TFAP2B','TFEC','TFPI','TFR2','TG','THRA','TLR4','TM6SF2','TMEM18','TMPRSS6','TPO','TSLP','TTN','TTPA','TTR','TUB','TXNRD2','TXNRD3','TYMS','U3','UBE2L3','UCP2','UCP3','UMOD','UNG','UST','WNK4','XDH','XK','ZBTB7B','ZFHX3','ZFR2','ZNF169','ZNF664'])

pos_control = pd.Series(["ACSL1","ADAM17","ADRB2","ADRB3","ADD1","ADH1B","ADH1C","ADIPOQ","ADORA2A","ADRB1","AGT","ALDH2","ALDOB","ALPL","APOA1","APOA2","APOA5","APOB","APOC3","APOE","BCO1","C3","CBS","CD36","CETP","CHDH","CHRNA3","CHRNA5","COMT","CRP","CSK","CYP11B2","CYP1A2","CYP2R1","DNMT3B","DRD2","EDN1","EPHX1","F5","FABP2","FADS1","FAF1","FOXO3","FTH1","FTO","FUT2","G6PD","GC","GHRL","GIPR","GPX1","GSTP1","HFE","HJV","HLA-DQA1","IL10","IL1B","IL6","INSIG2","LEPR","LIPC","LPA","LPL","MC4R","MCM6","MTHFR","MTNR1B","MTR","MTRR","NADSYN1","NBPF3","NOS3","NPY","PCSK1","PEMT","PPARA","PPARG","PPARGC1A","SELE","SH2B1","SIRT1","SLC23A1","SLC2A2","SLC30A8","SLC40A1","SLC5A6","SOD2","SOD3","TAS2R38","TCF7L2","TCN2","TFR2","TNF","TTPA","UGT1A1","VDR","PNPLA3","AGTR1","AOC1","BCMO1","CDKN2B-AS1","CRY1","CYBA","FGB","GHSR","GSTT1","HAMP","HLA-DQ","HLA-DQB1","HMGCR","IL6R","IRS1","LEP","MMP3","PON1","PROCR","SLC6A4","SOUX","SREBF2"])

neg_control = pd.Series(["ABCB1","ABCB11","AGER","AGXT","AIRE","ANGPTL7","APOC1","ARMS2","ASXL1","BACH2","BAX","BMPR2","BRIP1","BTG4","CACNA1S","CACNB2","CAPN3","CARD8","CASP3","CD40","CDK4","CDKN2A","CFH","CHRNB4","CLEC16A","CLEC7A","CLPTM1L","CNNM2","COL11A1","COLEC10","CPS1","CR1","CRHR1","CXCL12","CYCS","DDC","DNAH8","DRD1","DRD4","EPM2AIP1","ERBB2","ERBB4","F9","FANCA","FCGR3A","FCGR3B","FGFR4","FKRP","GCH1","GFAP","GLB1","GRIN2B","HERC2","HMOX1","HOXB13","HSPB1","ICAM4","IDH1","IL1R1","IRF6","ITGB3","KIT","LIG4","LSP1","MC1R","MED12","MMP2","MMP9","MMUT","MPZ","MT-ND4","MT-ND4L","MT-ND5","MUC1","MYBPC3","MYH7","MYH9","NBN","NF1","NFKBIA","NLGN3","NOS1","NOTCH3","NR1I3","OBSL1","OCA2","OTC","PBRM1","PCSK9","PDGFRB","PIK3R1","PLCE1","PMS2","POLD1","PON2","PRPH2","PSRC1","PTCH1","PTEN","PTPN2","RETN","RHO","RNASEL","RORA","ROS1","SCN4A","SERPINC1","SF3B1","SLC22A12","SLC22A2","SLC26A4","SLC2A9","SMAD4","SMO","SMPD1","SORL1","SQSTM1","SRD5A2","STK11","SYN3","TARDBP","TERT","TIRAP","TMPRSS2","TNFSF15","TNPO3","TOMM40","TREM2","UBC","UCHL1","VWF","WRN","WT1","XPC","ZFHX3","ZNF804A"])
neg_control = neg_control.sample(119)

#threshold_GRPM_report_sort.gene.sample(200).to_clipboard()
neg_control2 = pd.Series(["ABCA7","ABCB11","ACADVL","ACTA2","ADA","ADAMTS13","AKT1","AOPEP","ARMS2","ATP2B1","AURKA","BAG6","BARD1","BCL2","BMPR2","BRCA1","BRIP1","BTK","C2","CACNA1C","CACNB2","CAPN3","CCR5","CD14","CD209","CDKAL1","CFB","CLU","CNNM2","CNR1","CNTNAP2","CPS1","CSMD1","CTNNA2","CX3CR1","CXCL8","CYP17A1","CYP1B1","CYP2B6","CYP2E1","DBH","DDC","DISC1","DLG4","DNMT1","DRD3","EFEMP1","ERAP1","ERCC4","F8","FAM171A2","FCGR3A","FCGR3B","FGF23","FLG","FOXP2","FSHR","GALNS","GALT","GATA4","GDF5","GHR","GJB2","GNAS","GRIN2A","GRIN2B","GZMM","HNF1A","HRAS","HTR1B","HTR2A","HTRA1","IL18","IL1A","IL1RL1","IL2","IL2RA","IL7R","INHA","JAK2","KCNH2","KCNJ5","KCNQ2","LHCGR","LRRC37A2","MED12","MED12L","MEFV","MLH1","MMP9","MPL","MT-ATP6","MT-ATP8","MT-CO3","MT-CYB","MT-ND2","MUC1","MYH14","NOS2","NOTCH4","NPHS2","NR3C2","NRAS","NT5C2","NTRK1","OCA2","OPA1","PAH","PDGFRB","POLG","PON2","POTEF","PRPH2","PSEN2","PTEN","PTPN2","PTPN22","RBFOX1","RET","RHD","RHO","ROS1","SCN1A","SCN9A","SERPINE1","SFTPC","SLC26A4","SLC2A9","SLC6A3","SOD1","SORL1","SYNE1","TGFBR1","TLR1","TNFRSF11B","TNFRSF1B","TNPO3","TOMM40","TRIM37","TSBP1","TSC1","TSHR","TSPO","TYR","VAPB","VEGFA","WRAP53","WRN","WT1","XRCC2","ZAR1L","ZBTB12","ZNF804A"])
neg_control2 = neg_control2.sample(119)

len(genesthr02),len(genesthr01), len(genestesi),len(genespnn)


# In[ ]:


# THE TRUTH
genelist = genesthr02

mask = genelist.isin(pos_control)
genesin = genelist[mask]

mask2 = pos_control.isin(genesin)
genesnotin = pos_control[-mask2]
print('positive control')
len(genesin), len(genesnotin)


# In[ ]:


mask = genelist.isin(neg_control)
genesin2 = genelist[mask]

mask2 = neg_control.isin(genesin2)
genesnotin2 = neg_control[-mask2]
print('negative control')
len(genesin2), len(genesnotin2)


# In[ ]:


mask = genelist.isin(neg_control2)
genesin2 = genelist[mask]

mask2 = neg_control2.isin(genesin2)
genesnotin2 = neg_control2[-mask2]
print('negative control 2')
len(genesin2), len(genesnotin2)


#                 geni tesi (87)		geni PNN (340)
#     GRPM-X	a^2/bc	geni (thr01) (459)	65 (+)	22 (-)	109 (+)	231 (-)
#             geni (thr005) (879)	72 (+)	15 (-)	150 (+)	190 (-)
#             ad/bc	geni (top 459)	20 (+)	67 (-)	34 (+)	306 (-)
#             geni (top 879)	28 (+)	57 (-)	73 (+)	267 (-)
# 

# 

# In[ ]:


# CUSTOM RANGE----------------------------
#Matching PMIDs in Database
sort_by = 'interest_index'
threshold_GRPMX_report_int_sort = threshold_GRPMX_report_int_sort.sort_values(by= sort_by,ascending=False)
print('Sorted by: ', sort_by)
m = 1
n = 10
seed = 1
base = 0
y_axis = 'matching_mesh'
x = threshold_GRPMX_report_int_sort.gene.iloc[base*m:seed*(n*2)]
y = threshold_GRPMX_report_int_sort[y_axis].iloc[base*m:seed*(n*2)]
plt.figure(figsize=(4, len(x)*0.25))
plt.title('Matching PMIDs in Database\n sorted by '+sort_by, loc='center',pad=10)

plt.barh(x,y, color = '#04ade3')
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel(y_axis, position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


#GET A SAMPLE ------------------------------------------------------------------
## get top 10 Genes:
top10_list = GRPMX_report_int_sort['gene'].iloc[:10].tolist()#.to_clipboard()
pyperclip.copy(str(top10_list))
top10_list


# 

# 
# 
# 

# # IMPORT DATA FROM DATABASE

# ##gene lookup threshold based

# In[ ]:


threshold_GRPMX_report_int_sort.iloc[:30]


# ### SINGLE GENE LOOKUP

# # IMPORT NBIB DATA FROM DATABESE
# time1 = datetime.now()
# #import gene-fullnbib
# dummy_nbib = pd.read_csv(db_name+'/complete_nbibtable.csv', index_col=0)
# dummy_nbib['pubmed_id'] = dummy_nbib['pubmed_id'].astype(str)
# time2 = datetime.now()
# print('time import nbib: ', time2-time1)
# print(dummy_nbib.memory_usage().sum() / 1024 / 1024, 'MB')

# #gene nbib lookup
# gene = 'APOA1'
# #gene = GRPM_report_sort.gene.iloc[1]
# 
# gene_nbib = dummy_nbib.loc[dummy_nbib['gene'] == gene]
# gene_nbib['descriptors']#.iloc[0]
# gene_nbib

# #mask full nbib with filtered pmids for gene query:
# mask = gene_nbib.pubmed_id.isin(gene_filtered_grpm.pmids)
# gene_nbib_filtered = gene_nbib[mask]
# gene_nbib_filtered.abstract#.to_clipboard()
# 
# #filter for rsid
# gene_filtered_grpm_rsid = gene_filtered_grpm[gene_filtered_grpm.rsid == 'rs2266788']
# 
# mask = gene_nbib.pubmed_id.isin(gene_filtered_grpm_rsid.pmids)
# gene_nbib_filtered_rsid = gene_nbib[mask]
# 
# gene_nbib_filtered.abstract#.to_clipboard()

# #FTO_nbib[['pubmed_id','publication_date','title','abstract']].to_csv('FTO_title-abstracts.csv')
# gene_pmids = list(gene_nbib['pubmed_id'].drop_duplicates())#.to_clipboard(index=False) # --> to VEP http://www.ensembl.org/Tools/VEP
# query = " OR ".join(gene_pmids)
# pyperclip.copy(query)

# In[ ]:


# IMPORT GRPM DATA FROM DATABESE-----------------------------------------------------
time2 = datetime.now()
#import gene-rsidpmidmesh
dummy_grpm = pd.read_csv(db_name+'/grpm_table_output.csv', index_col=0)
dummy_grpm['pmids'] = dummy_grpm['pmids'].astype(str) #convert pmid type in str
time3 = datetime.now()
print('time import grpm: ', time3-time2)
print(dummy_grpm.memory_usage().sum() / 1024 / 1024, 'MB')


# In[ ]:


#GRPM DB statistics
print('GRPM DB statistics:')
grpmdb_genes = dummy_grpm.gene.nunique()
grpmdb_rsids = dummy_grpm.rsid.nunique()
grpmdb_pmids = dummy_grpm.pmids.nunique()
grpmdb_meshs = dummy_grpm.mesh.nunique()
print(grpmdb_genes,'genes')
print(grpmdb_rsids,'rsid')
print(grpmdb_pmids,'pmid')
print(grpmdb_meshs,'mesh')
dummy_grpm


# In[ ]:


#gene grpm lookup
gene_grpm = dummy_grpm.loc[dummy_grpm['gene'] == gene]
gene_grpm = gene_grpm#[['gene', 'rsid', 'pmids', 'mesh', 'qualifier', 'major']]#.drop_duplicates().reset_index(drop=True)

gene_grpm#.head(len(gene_pmidrsidmesh))


#     # old merger
#     #import gene-rsidpmid
#     gene_litvar1_rsids2pmids = pd.read_csv(db_name+'/'+gene+'_litvar1_rsids2pmids.csv').drop(columns=['Unnamed: 0'])
#     gene_litvar1_rsids2pmids['pmids'] = gene_litvar1_rsids2pmids['pmids'].astype(str)
# 
#     # add a rsid-merger to gene_pubmed_pmidmesh
#     gene_rsidpmidmesh = pd.merge(gene_litvar1_rsids2pmids, gene_pubmed_pmidmesh, on='pmids')
#     gene_rsidpmidmesh

# # IMPORT GRPMX data from Survey

# #import grpmx_output
# dummy_grpm_filter = pd.read_csv(directory+'/grpmx_output.csv', index_col=0)
# dummy_grpm_filter['pmids'] = dummy_grpm_filter['pmids'].astype(str) #convert pmid type in str
# 
# #gene lookup
# gene = ''
# gene_grpm_filter = dummy_grpm_filter.loc[dummy_grpm_filter['gene'] == gene]
# gene_filtered_grpm = gene_grpm_filter.drop_duplicates().reset_index(drop=True)
# 
# len_pmids_01 = len(gene_filtered_grpm.pmids.drop_duplicates())
# print('interesting pmids for gene '+gene+':', len_pmids_01)
# 
# #to clipboard()
# gene_filtered_grpm.rsid.to_clipboard(index = False)
# gene_filtered_grpm.rsid.drop_duplicates().to_clipboard()
# 
# #module isin  on gene_rsidpmid
# mask = gene_grpm.mesh.isin(gene_filtered_grpm.mesh)
# gene_filtered_pmidrsidmesh_full = gene_grpm[mask].reset_index(drop=True)
# gene_filtered_pmidrsidmesh_full
# #gene_filtered_pmidrsidmesh

# #containing word Mesh-LOOKUP
# gene_grpm[gene_grpm.mesh.str.contains('poly', case=False)]

# In[ ]:


# IMPORT GRPMX DATA FROM SURVEY-----------------------------------------------------
time2 = datetime.now()
#import gene-rsidpmidmesh
dir1_grpmx = pd.read_csv(directory + '/grpmx_output.csv', index_col=0)
dir1_grpmx['pmids'] = dir1_grpmx['pmids'].astype(str) #convert pmid type in str
time3 = datetime.now()
print('time import grpmx: ', time3-time2, directory)


# In[ ]:


print('general grpmx statistics',directory)

dir1_grpmx_gene = dir1_grpmx.gene.nunique()
dir1_grpmx_rsid = dir1_grpmx.rsid.nunique()
dir1_grpmx_pmids = dir1_grpmx.pmids.nunique()
dir1_grpmx_mesh = dir1_grpmx.mesh.nunique()
print(dir1_grpmx_gene ,'gene match')
print(dir1_grpmx_rsid ,'rsid match')
print(dir1_grpmx_pmids,'pmid match')
print(dir1_grpmx_mesh ,'mesh match on',GRPMX_report.reference_mesh[0])

print('gene', round(dir1_grpmx_gene /grpmdb_genes,2), 'marching index on grpm db')
print('rsid', round(dir1_grpmx_rsid /grpmdb_rsids,2), 'marching index on grpm db')
print('pmid', round(dir1_grpmx_pmids/grpmdb_pmids,2), 'marching index on grpm db')
print('mesh', round(dir1_grpmx_mesh /grpmdb_meshs,2), 'marching index on grpm db')

dir1_grpmx


# In[ ]:


# IMPORT GRPMX DATA FROM SURVEY-----------------------------------------------------
time2 = datetime.now()
#import gene-rsidpmidmesh
dir2_grpmx = pd.read_csv(directory_2 + '/grpmx_output.csv', index_col=0)
dir2_grpmx['pmids'] = dir2_grpmx['pmids'].astype(str) #convert pmid type in str
time3 = datetime.now()
print('time import grpmx: ', time3-time2, directory)


# In[ ]:


print('general grpmx statistics',directory_2)

dir2_grpmx_gene = dir2_grpmx.gene.nunique()
dir2_grpmx_rsid = dir2_grpmx.rsid.nunique()
dir2_grpmx_pmids = dir2_grpmx.pmids.nunique()
dir2_grpmx_mesh = dir2_grpmx.mesh.nunique()
print(dir2_grpmx_gene ,'gene match')
print(dir2_grpmx_rsid ,'rsid match')
print(dir2_grpmx_pmids,'pmid match')
print(dir2_grpmx_mesh ,'mesh match on',GRPMX_report.reference_mesh[0])

print('gene', round(dir2_grpmx_gene /grpmdb_genes,2), 'marching index on grpm db')
print('rsid', round(dir2_grpmx_rsid /grpmdb_rsids,2), 'marching index on grpm db')
print('pmid', round(dir2_grpmx_pmids/grpmdb_pmids,2), 'marching index on grpm db')
print('mesh', round(dir2_grpmx_mesh /grpmdb_meshs,2), 'marching index on grpm db')

dir2_grpmx


# 

#     #### txt/csv MERGER Trial
#     # simple csv merger in column
#     import os
#     import csv
# 
#     # Define the directory containing the CSV files
#     directory = r'C:\Users\Public\Database_SSD\GRPM-X_proteincodinggens_db\grpm_db'
# 
#     genes =  ['ALG1', 'ALG12', 'ALG2', 'ALG3', 'ALG6', 'ALG8', 'ALG9', 'DDOST', 'DOLK', 'MPDU1', 'PGM1', 'PMM1', 'PMM2', 'RFT1']
# 
#     endswitch = '_pubmed_pmidmesh'
#     # Define the name of the merged output file
#     output_file = endswitch+".csv"
# 
# 
#     # Open the output file in write mode
#     with open(output_file, "w", newline="") as outfile:
#         for gene in genes:
#             # Create a CSV writer object
#             writer = csv.writer(outfile)
# 
#             # Iterate through each file in the directory
#             for filename in os.listdir(directory):
# 
#                 # Check if the file is a CSV file and ends with "_filtered_pmidrsidmesh.csv"
#                 if filename.endswith(gene+endswitch+".csv"):
#                     #if filename.endswith(".csv"):
# 
#                     # Open the file in read mode
#                     with open(os.path.join(directory, filename), "r", newline="") as infile:
# 
#                         # Create a CSV reader object
#                         reader = csv.reader(infile)
# 
#                         # Write the filename as the first row of the output file
#                         writer.writerow([filename])
# 
#                         # Write the contents of the file to the output file
#                         for row in reader:
#                             writer.writerow(row)
# 
#                         # Write a blank row between files
#                         writer.writerow([])
# 

# ## CSV merger (gene/rsid/pmid/mesh)
# #### converting old data storage system
# #goal
# [gene,rsid,pmids,mesh]

# In[ ]:


genes = GRPMX_report.gene
genes


# In[ ]:


genes = GRPMX_report.gene
genes.index[genes == gene][0]


#     #with append (deprecated)
#     complete_df = pd.DataFrame()
#     genes = ['FTO','INS']
#     #genes = GRPMX_report.gene
#     for gene in genes:
#         if os.path.isfile('grpm_survey/'+gene+'_filtered_pmidrsidmesh.csv'):
#             dfmatch_gene = pd.read_csv('grpm_survey/'+gene+'_filtered_pmidrsidmesh.csv', index_col=0)#).set_index('Unnamed: 0')
#             dfmatch_gene['gene'] = gene
#             dfmatch_gene = dfmatch_gene.reindex(columns=['gene','rsid', 'pmids', 'mesh'])
#             complete_df = complete_df.append(dfmatch_gene)
#             #print(genes.index[genes == gene][0])
# 

# #### Merge Survey

# In[ ]:


## Survey:

#with concat:
complete_df = pd.DataFrame()
#genes = ['FTO','INS']
genes = GRPMX_report.gene
for gene in genes[:10]:
    if os.path.isfile('grpm_survey/'+gene+'_filtered_pmidrsidmesh.csv'):
        dfmatch_gene = pd.read_csv('grpm_survey/'+gene+'_filtered_pmidrsidmesh.csv', index_col=0)
        dfmatch_gene['gene'] = gene

        complete_df = pd.concat([complete_df, dfmatch_gene])
        print(genes.index[genes == gene][0])

complete_df = complete_df.reindex(columns=['gene','rsid', 'pmids', 'mesh'])
complete_df.to_csv('complete_surveytable_genersidpmidmesh.csv')


# In[ ]:


complete_surveytable_genersidpmidmesh = pd.read_csv('complete_surveytable_genersidpmidmesh.csv', index_col=0)
complete_surveytable_genersidpmidmesh.gene.drop_duplicates()


# #### Merge grpm db: grpm table and nbib

# In[ ]:


## nbib db merger -----------------------------

#with concat:
complete_nbib = pd.DataFrame()
#genes = ['FTO','INS']
genes = GRPM_report.gene#.iloc[:5]

if os.path.isfile(db_name+'complete_nbibtable.csv'):
    complete_nbib = pd.read_csv(db_name+'complete_nbibtable.csv', index_col=0)
else:
    for gene in genes:
        if os.path.isfile(db_name+'/'+gene+'_nbib_table.csv'):
            #import gene_nbib
            gene_nbib = pd.read_csv(db_name+'/'+gene+'_nbib_table.csv', index_col=0)
            gene_nbib['pubmed_id'] = gene_nbib['pubmed_id'].astype(str)

            gene_nbib['gene'] = gene

            complete_nbib = pd.concat([complete_nbib, gene_nbib])
            print(genes.index[genes == gene][0])


# In[ ]:


pyperclip.copy(str(complete_nbib.columns.to_list()))

if os.path.isfile('complete_nbibtable.csv'):
    pass
else:
    complete_nbib = complete_nbib.reindex(columns=['gene','pubmed_id', 'citation_owner', 'nlm_status', 'last_revision_date', 'electronic_issn', 'print_issn', 'linking_issn', 'journal_volume', 'journal_issue', 'publication_date', 'title', 'pages', 'abstract', 'copyright', 'authors', 'language', 'grants', 'publication_types', 'electronic_publication_date', 'place_of_publication', 'journal_abbreviated', 'journal', 'nlm_journal_id', 'descriptors', 'pmcid', 'keywords', 'conflict_of_interest', 'received_time', 'revised_time', 'accepted_time', 'pubmed_time', 'medline_time', 'entrez_time', 'pii', 'doi', 'publication_status', 'corporate_author', 'secondary_source'])
    #complete_nbib.to_csv('grpm_survey/surveytable_genersidpmidmesh.csv')
    complete_nbib.to_csv('complete_nbibtable.csv')


# In[ ]:


print(len(complete_nbib.pubmed_id.drop_duplicates()))
complete_nbib


# In[ ]:


## grpm table db Merger--------------------

#with concat:
complete_df = pd.DataFrame()
#genes = ['FTO','INS']
genes = GRPM_report.gene

if os.path.isfile(db_name+'complete_dbtable_genersidpmidmesh.csv'):
    complete_df = pd.read_csv(db_name+'grpm_table_output.csv')
else:
    for gene in genes:
        if os.path.isfile(db_name+'/'+gene+'_pubmed_pmidmesh.csv'):
            #import gene-rsidpmidmesh
            gene_pubmed_pmidmesh = pd.read_csv(db_name+'/'+gene+'_pubmed_pmidmesh.csv', index_col=0)
            gene_pubmed_pmidmesh['pmids'] = gene_pubmed_pmidmesh['pmids'].astype(str)

            #import gene-rsidpmid
            gene_litvar1_rsids2pmids = pd.read_csv(db_name+'/'+gene+'_litvar1_rsids2pmids.csv', index_col=0)
            gene_litvar1_rsids2pmids['pmids'] = gene_litvar1_rsids2pmids['pmids'].astype(str)

            # add a rsid-merger to gene_pubmed_pmidmesh
            gene_rsidpmidmesh = pd.merge(gene_litvar1_rsids2pmids, gene_pubmed_pmidmesh, on='pmids')
            gene_rsidpmidmesh['gene'] = gene

            complete_df = pd.concat([complete_df, gene_rsidpmidmesh])
            print(genes.index[genes == gene][0])

if os.path.isfile('grpm_table_output.csv'):
    pass
else:
    complete_df = complete_df.reindex(columns=['gene','rsid', 'pmids', 'mesh', 'qualifier', 'major'])
    #complete_df.to_csv('grpm_survey/surveytable_genersidpmidmesh.csv')
    complete_df.to_csv('grpm_table_output.csv')
    #complete_df#.columns


# ## Complete GRPM DB TABLE
# 

# In[ ]:


## visualize complete dbtable:
complete_dbtable_grpm = pd.read_csv('grpm_table_output.csv', index_col=0)
complete_dbtable_grpm.gene.drop_duplicates()


# In[ ]:


gene = 'G6PD'
mask = complete_dbtable_grpm.gene == (gene)
LEP_lookup = complete_dbtable_grpm[mask]
LEP_lookup.mesh.drop_duplicates()


# In[ ]:


# mesh lookup for rsid
mask = complete_dbtable_grpm.mesh.str.contains('utrigen')
mesh_lookup = complete_dbtable_grpm[mask]

#complete_dbtable_grpm.rsid
mesh_lookup.gene.drop_duplicates().to_clipboard()
mesh_lookup


# In[ ]:


## Count with groupby-describe------------------

#complete_df_count =
type(complete_df.pmids.iloc[0])

#Analyze enrichment with groupby.describe method------------------------
complete_df['pmids'] = complete_df['pmids'].astype(str)

### groupby.describe analysis by rsid
complete_df_mesh_count = complete_df.groupby('mesh').describe().reset_index()
#pyperclip.copy(str(complete_df_mesh_count.columns.to_list()))

complete_df_mesh_count.columns = complete_df_mesh_count.columns.to_flat_index()
new_column_names = [('mesh'), ('gene_count'), ('gene_unique'), ('gene_top'), ('gene_freq'), ('rsid_count'), ('rsid_unique'), ('rsid_top'), ('rsid_freq'), ('pmids_count'), ('pmids_unique'), ('pmids_top'), ('pmids_freq'), ('qualifier_count'), ('qualifier_unique'), ('qualifier_top'), ('qualifier_freq'), ('major_count'), ('major_unique'), ('major_top'), ('major_freq')]
complete_df_mesh_count.columns = new_column_names
#complete_df_mesh_count.to_csv('grpm_survey/surveytable_genersidpmidmesh_countbymesh.csv')


# In[ ]:


complete_df_mesh_count


# In[ ]:


complete_df.pmids.drop_duplicates()


# In[ ]:


complete_df_mesh_count


# 

# ## Density Plots + Rug (grpm_report)

# In[ ]:


import numpy as np
import seaborn as sns
data = complete_df_mesh_count.pmids_count
sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)


# In[ ]:


import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

data = GRPMX_report_int.matching_pmids_ratio

sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)
sns.rugplot(np.array(data), color='r')

plt.title('Ig = (Pm/P)g')
plt.title('Interest index = a^2/bc')
plt.show()


# In[ ]:


GRPMX_report.matching_pmids_ratio


# ## Add study type with API

# In[ ]:


def get_study_type(pmids):

    Entrez.email = 'your_email@your_domain.com'

    # Retrieve the metadata for the articles
    handle = Entrez.esummary(db='pubmed', id=','.join(pmids), retmode='xml')
    records = Entrez.parse(handle)

    # Extract the article types from the metadata
    study_types = []
    for record in records:
        article_types = record['PubTypeList']
        #print(record['PubTypeList'])
        # Determine the study type based on the article types
        if 'Randomized Controlled Trial' in article_types:
            study_types.append('Randomized Controlled Trial')
        elif 'Controlled Clinical Trial' in article_types:
            study_types.append('Controlled Clinical Trial')
        elif 'Cohort Studies' in article_types:
            study_types.append('Cohort Study')
        elif 'Case-Control Studies' in article_types:
            study_types.append('Case-Control Study')
        elif 'Review' in article_types:
            study_types.append('Review')
        elif 'Clinical Trial' in article_types:
            study_types.append('Clinical Trial')
        elif 'Meta-Analysis' in article_types:
            study_types.append('Meta-Analysis')
        elif 'Multicenter Study' in article_types:
            study_types.append('Multicenter Study')
        # Add additional conditions to handle other study types as needed
        else:
            study_types.append('Unknown')

    return study_types


# In[ ]:


Entrez.email = 'your_email@your_domain.com'

# Retrieve the metadata for the articles
handle = Entrez.esummary(db='pubmed', id=','.join(['34556834', '26620191', '33006084']), retmode='xml')
records = Entrez.parse(handle)
dfa =pd.DataFrame()
for i in records:
    df = pd.json_normalize(i)
    dfa = pd.concat([dfa,df])

dfa.T


# In[ ]:


# Extract the article types from the metadata
study_types = []
for record in records:
        article_types = record['PubTypeList']
        print(record['PubTypeList'])
        # Determine the study type based on the article types
        if 'Randomized Controlled Trial' in article_types:
            study_types.append('Randomized Controlled Trial')
        elif 'Controlled Clinical Trial' in article_types:
            study_types.append('Controlled Clinical Trial')
        elif 'Cohort Studies' in article_types:
            study_types.append('Cohort Study')
        elif 'Case-Control Studies' in article_types:
            study_types.append('Case-Control Study')
        elif 'Review' in article_types:
            study_types.append('Review')
        # Add additional conditions to handle other study types as needed
        else:
            study_types.append('Unknown')


# In[ ]:


get_study_type(['34556834', '31533339', '33006084'])


# In[ ]:


#it takes time
ftopmids_str = list(map(str, ftopmids))
study_type = get_study_type(ftopmids_str)


# In[ ]:


pmids_studytype = pd.DataFrame(list(zip(ftopmids_str,study_type)),columns=[gene+'_PMID','study type'])
pmids_studytype_count = pmids_studytype.groupby('study type').describe().reset_index()
pmids_studytype_count.columns = pmids_studytype_count.columns.to_flat_index()
new_column_names = ['study_type', 'pmid-count', 'pmid-unique','pmid-top','pmid-freq']
pmids_studytype_count.columns = new_column_names
pmids_studytype_count


# In[ ]:


rand = pmids_studytype.loc[pmids_studytype['study type']=='Randomized Controlled Trial']
#rand = pmids_studytype.loc[pmids_studytype['study type']!='Unknown']
#rand = pmids_studytype.loc[pmids_studytype['study type']=='Unknown']
rand_gene_PMID = list(rand[gene+'_PMID'])
#conseqf = conseq.loc[conseq['SYMBOL'] == gene]
len(rand_gene_PMID)

#quert on pubmed:
pyperclip.copy(" OR ".join(list(map(str,rand_gene_PMID))))


# In[ ]:


#controllo
rand_FTO_PMID_str = list(map(str, rand.FTO_PMID))
study_type_control = get_study_type(rand_FTO_PMID_str)


# In[ ]:


#study_type_control


# #### Abstract analysis

# In[ ]:


import gensim
from gensim.utils import simple_preprocess
from gensim.parsing.preprocessing import STOPWORDS
from gensim import corpora
import nltk
from nltk.stem import PorterStemmer

# create a Porter stemmer object
p_stemmer = PorterStemmer()

# list of abstracts
abstracts = ['It has been suggested that Neel s "thrifty genotype" model may account for high body weights in some Oceanic populations, which presumably arose in modern times. In European populations, common variants (rs1421085-C, rs17817449-G, and rs9939609-A) in the fat mass and obesity (FTO associated) were recently found to be associated with body mass index (BMI) or obesity. In this study, we investigated the population frequencies of these variants in six Oceanic populations (Melanesians, Micronesians, and Polynesians) and tested for an association with BMI. Unlike European populations, the Oceanic populations displayed no significant association between the FTO polymorphisms and BMI. These variants were in strong linkage disequilibrium. The population frequencies ranged between 4.2 and 30.3% in the six Oceanic populations, and were similar to those in southeast and east Asian populations. Our study of the FTO polymorphisms has generated no evidence to support the thrifty genotype hypothesis for Oceanic populations.',
             'Variations in the fat-mass and obesity-associated gene (FTO) are associated with the obesity phenotype in many Caucasian populations. This association with the obesity phenotype is not clear in the Japanese. To investigate the relationship between the FTO gene and obesity in the Japanese, we genotyped single nucleotide polymorphisms (SNPs) in the FTO genes from severely obese subjects [n = 927, body mass index (BMI) > or = 30 kg/m2] and normal-weight control subjects (n = 1,527, BMI < 25 kg/m2). A case-control association analysis revealed that 15 SNPs, including rs9939609 and rs1121980, in a linkage disequilibrium (LD) block of approximately 50 kb demonstrated significant associations with obesity',
             'BACKGROUND: Recently, the role of the FTO (fat mass and obesity associated) gene in obesity development was described in Western European, but not in Oceanic, cohorts. The objective of this study was to test the hypothesis that the FTO single nucleotide polymorphism (SNP) is associated with body mass index (BMI) in the Slavic population and to analyze if there could be sex-specific effects of the SNP on BMI, waist-to-hip ratio (WHR), and lipid parameters. METHODS: We analyzed three large population-based samples comprising the post-MONICA study (1191 males, 1368 females) and the 3PMFs study (908 females). RESULTS: FTO rs17817449 SNP was related to BMI in males (p=0.014). In the females from both the post-MONICA and the 3PMFs study, FTO had no effect on BMI. Sub-analysis of females from the 3PMFs study demonstrated that FTO had an effect on BMI in postmenopausal females (p=0.035) but not in premenopausal females (follicle-stimulating hormone <40 U/L was used as marker of premenopausal status). WHR and lipid parameters were not associated with FTO in any of the analyzed groups. CONCLUSIONS: These results suggest that the effect of FTO SNP rs17817449 may be, in some populations at least, restricted to males and postmenopausal females.',
             'Participants analyzed actual and simulated longitudinal data from the Framingham Heart Study for various metabolic and cardiovascular traits. The genetic information incorporated into these investigations ranged from selected single-nucleotide polymorphisms to genome-wide association arrays. Genotypes were incorporated using a broad range of methodological approaches including conditional logistic regression, linear mixed models, generalized estimating equations, linear growth curve estimation, growth modeling, growth mixture modeling, population attributable risk fraction based on survival functions under the proportional hazards models, and multivariate adaptive splines for the analysis of longitudinal data. The specific scientific questions addressed by these different approaches also varied, ranging from a more precise definition of the phenotype, bias reduction in control selection, estimation of effect sizes and genotype associated risk, to direct incorporation of genetic data into longitudinal modeling approaches and the exploration of population heterogeneity with regard to longitudinal trajectories. The group reached several overall conclusions: (1) The additional information provided by longitudinal data may be useful in genetic analyses. (2) The precision of the phenotype definition as well as control selection in nested designs may be improved, especially if traits demonstrate a trend over time or have strong age-of-onset effects. (3) Analyzing genetic data stratified for high-risk subgroups defined by a unique development over time could be useful for the detection of rare mutations in common multifactorial diseases. (4) Estimation of the population impact of genomic risk variants could be more precise. The challenges and computational complexity demanded by genome-wide single-nucleotide polymorphism data were also discussed.',
             'This is the fifth abstract.']

# function to preprocess the abstracts
def preprocess_abstracts(abstracts):
    for abstract in abstracts:
        # tokenize the abstract into individual words
        tokens = simple_preprocess(abstract, deacc=True)
        # remove stop words
        stopped_tokens = [token for token in tokens if not token in STOPWORDS]
        # stem the tokens
        stemmed_tokens = [p_stemmer.stem(token) for token in stopped_tokens]
        yield stemmed_tokens

# preprocess the abstracts
processed_abstracts = list(preprocess_abstracts(abstracts))

# create a dictionary from the processed abstracts
dictionary = corpora.Dictionary(processed_abstracts)

# create a corpus from the dictionary and the processed abstracts
corpus = [dictionary.doc2bow(abstract) for abstract in processed_abstracts]

# train an LSI (Latent Semantic Indexing) model on the corpus
lsi_model = gensim.models.LsiModel(corpus=corpus, id2word=dictionary, num_topics=2)

# print the top 5 topics
print(lsi_model.print_topics(num_topics=8, num_words=20))

