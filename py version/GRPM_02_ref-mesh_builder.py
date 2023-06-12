#!/usr/bin/env python
# coding: utf-8

# # Ref MeSH List Builder

#  This notebook allows you to create Mesh lists for the GRPMX system starting with a topic of choice.
#  The system is based on the complete MESH datasheet.
#  The system uses ChatGTP to provide the initial "meshes."

# # Import Modules

# In[ ]:


import csv
#Import Modules
import os
import glob
import sys
import json
import openai
import IPython
#import nbib
#import requests as rq
from datetime import datetime
import time
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pyperclip

#from bs4 import BeautifulSoup
#from Bio import Entrez
#Entrez.email = "your_email@example.com"


# # Load Mesh dataset

# In[ ]:


#load MESH dataframe:
if os.path.exists('ref-mesh-archive/MESH.csv'):
    df = pd.read_csv('ref-mesh-archive/MESH.csv')
else:
    print("no MESH.csv avalable jump to 'Import/create MESH-STY-LITVAR' section")

#info on Proprieties: https://bioportal.bioontology.org/ontologies/MESH?p=properties
#Mesh Browser: https://meshb.nlm.nih.gov/search?searchMethod=SubString&searchInField=termDescriptor&sort=&size=20&searchType=exactMatch&from=0&q=Diabetes

# AQL: Allowable Qualifiers
AQL = ["blood (BL)","cerebrospinal fluid (CF)","chemically induced (CI)","classification (CL)","complications (CO)","congenital (CN)","diagnosis (DI)","diagnostic imaging (DG)","diet therapy (DH)","drug therapy (DT)","economics (EC)","embryology (EM)","enzymology (EN)","epidemiology (EP)","ethnology (EH)","etiology (ET)","genetics (GE)","history (HI)","immunology (IM)","metabolism (ME)","microbiology (MI)","mortality (MO)","nursing (NU)","parasitology (PS)","pathology (PA)","physiopathology (PP)","prevention & control (PC)","psychology (PX)","radiotherapy (RT)","rehabilitation (RH)","surgery (SU)","therapy (TH)","urine (UR)","veterinary (VE)","virology (VI)"]
aql=pd.DataFrame(AQL)
aql


# In[ ]:


print(df.columns)
df.iloc[:3]
df.T


# In[ ]:


#How many Mesh has mapping qualifier (subheadings)
all = df['Preferred Label'].drop_duplicates()
mapping_qf = df[['Mapping qualifier of','Preferred Label']].drop_duplicates().dropna()
has = df[['Has mapping qualifier','Preferred Label']].drop_duplicates().dropna()
has_len = has['Preferred Label'].nunique()
print('Total number of mesh terms:\n', df['Preferred Label'].nunique())
print('How many Mesh has mapping qualifier (subheadings):\n',has_len)
print(' %', round((has_len/len(all)*100),3))
#for i in df.columns:
#    print(df[i].drop_duplicates())


# In[ ]:


# search for single exact value
mesh = 'Disease or Syndrome'
print('search for single exact value:')
df[df['Preferred Label']==mesh]


# In[ ]:


#Extract subset of Mesh
str_1 = 'Fabry'
str_2 = 'Lysosomal Sto'
str_3 = ''
str_full = str_2 or str_1 or str_3

ref_search = df[df['Preferred Label'].str.contains(str_full).fillna(False)]
#sem.to_csv('ref-mesh-archive/ClassID-STY_SemanticTypes.csv')
print('look for mesh containing "',str_1,',',str_2,',',str_3,  '":')
ref_search
#ref_search.to_csv('ref-mesh-archive/ref_mesh_LSD.csv')


# ## Semantic Types Analysis

# In[ ]:


#Extract subset of all semantic types
df_sem = df[df['Class ID'].str.contains('STY').fillna(False)]
df_sem
#sem.to_csv('ref-mesh-archive/ClassID-STY_SemanticTypes.csv')
#df['Class ID'].to_csv('ref-mesh-archive/ClassID.csv')


# #### Characterize a Semantic Type

# In[ ]:


#how many mesh for semantic type?
semantictype = 'T192'
sam_label = df[df['Class ID'].str.contains(semantictype).fillna(False)].reset_index().at[0,'Preferred Label']
# Give me Semantic Type name for STY ID
print('semantic type:', sam_label)

mask2 = df['Semantic Types'].str.contains(semantictype).fillna(False)
mesh_sem = df[mask2][['Class ID','Preferred Label','Definitions','Semantic Types']]
print('number of mesh:', len(mesh_sem))
print('\nMesh for "',semantictype,sam_label,'" semantic type:')

mesh_sem


# #### Characterize single mesh

# In[ ]:


# search for single Mesh
mesh_classid ='D004048'
maskmesh = df['Class ID'].str.contains(mesh_classid).fillna(False)

sty_classid = df[maskmesh]['Semantic Types'].reset_index().iat[0, 1][-4:]
#print('note: this method catches only meshes that has sam_label at the in the felad "sematic Type" for multy sem type meshes')

masksty = df['Class ID'].str.contains(sty_classid).fillna(False)

masked = df[maskmesh][['Semantic Types','Preferred Label','Definitions']].reset_index()

#get mesh definition:
print('mesh:', masked.at[0,'Preferred Label'],
      '\nsemantic type:',df[masksty].reset_index().at[0,'Preferred Label'],sty_classid,
      '\ndescrition:',masked.at[0,'Definitions'])


# In[ ]:


# Row aggregation
data_mesh = []
data_sty = ['T021','T116']
dfd = pd.DataFrame()
for i in data_sty:
    #mask = df['Class ID'].str.contains(i).fillna(False) # aggregate sty
    mask = df['Semantic Types'].str.contains(i).fillna(False) # aggregate mesh for sty
    #dfd = dfd.append(df[mask])
    dfd = pd.concat([dfd, df[mask]])
dfd


# # Import/create MESH-STY and MESH-STY-LITVAR

# In[ ]:


#set options:
import_mesh_sty = True

# concat the exploded mesh table or import prearranged csv
if os.path.isfile('ref-mesh-archive/MESH_STY.csv') and import_mesh_sty == True:
    mesh_large_df_sty = pd.read_csv('ref-mesh-archive/MESH_STY.csv', index_col=0)
    print('MESH_STY sty imported from csv')
else:
    if import_mesh_sty == True:
        if os.path.exists('ref-mesh-archive/MESH.csv'):
            # correct list format
            df['STY_ID'] = df['Semantic Types'].str.replace(r'http://purl.bioontology.org/ontology/STY/','', regex = False)
            start_word = '[\"'
            end_word = '\"]'
            df['STY_ID'] = f'{start_word} ' + df['STY_ID'] + f' {end_word}'
            df['STY_ID'] = df['STY_ID'].str.replace(' ','', regex = False)
            df['STY_ID'] = df['STY_ID'].str.replace('|','\",\"', regex = False)
            print(df['STY_ID'].isna().sum())
            df.dropna(subset=['STY_ID'], inplace=True)
            df.reset_index(drop=True,inplace=True)
            type(df.STY_ID[0])

            from ast import literal_eval
            df['STY_ID'] = df['STY_ID'].apply(literal_eval)

            # Import exhcange table sty-code
            sty = pd.read_csv('ref-mesh-archive/MeshSTY-code.csv',sep=';')
            sty = sty.rename(columns={'ID':'Semantic Types'})
            sty = sty.rename(columns={'ID':'Semantic Types'})
            print(sty['Semantic Types'].nunique())

            #mesh_large = pd.DataFrame()
            mesh_large = []
            for i in range(len(df)):
                for sem in df['STY_ID'][i]: #dfrspost = mother table
                    out = df['Preferred Label'][i],df['Class ID'][i],sem,df['Synonyms'][i],df['Parents'][i],df['CUI'][i],df['AQL'][i],df['TERMUI'][i]
                    mesh_large.append(out)
                    #df_out = pd.DataFrame(out)
                    #pd.concat([mesh_large, df_out])

            mesh_large_df = pd.DataFrame(mesh_large)
            new_col_names = ['Preferred Label','Class ID','Semantic Types','Synonyms','Parents','CUI','AQL','TERMUI']
            mesh_large_df.columns = new_col_names
            mesh_large_df
            ## Add STY Labels
            mesh_large_df_sty = pd.merge(mesh_large_df, sty, on='Semantic Types', how='inner').reset_index(drop=True)
            #Add rsid coulmn con merge
            mesh_large_df_sty = mesh_large_df_sty.rename(columns={'Preferred Label_y':'Semantic Types Label','Preferred Label_x':'Preferred Label'})

            mesh_large_df_sty = mesh_large_df_sty[['Preferred Label', 'Semantic Types Label', 'Class ID', 'Semantic Types', 'Synonyms', 'Parents', 'CUI', 'AQL', 'TERMUI']]
            mesh_large_df_sty.to_csv('ref-mesh-archive/MESH_STY.csv')
            print('MESH_STY created')
        else:
            print('MESH.csv not avalable')


#### Create MESH-STY-LITVAR subset from MESH-STY.csv
if os.path.exists('ref-mesh-archive/MESH_STY_LITVAR1.csv'):
    mesh_litvar_sty = pd.read_csv('ref-mesh-archive/MESH_STY_LITVAR1.csv',index_col=0)
    print('MESH_STY_LITVAR1 imported from csv')
else:
    grpm_mesh = pd.read_csv('ref-mesh-archive/grpm_db_mesh.csv', index_col=0)
    mask = mesh_large_df_sty['Preferred Label'].isin(grpm_mesh.mesh)
    #mesh_large_df_sty['Preferred Label'].nunique()
    mesh_litvar_sty = mesh_large_df_sty[mask]
    mesh_litvar_sty.to_csv('ref-mesh-archive/MESH_STY_LITVAR1.csv')
    print('MESH_STY_LITVAR1 created')


# In[ ]:


if 'mesh_large_df_sty' in globals():
    print('MESH_STY', len(mesh_large_df_sty['Preferred Label']), ' rows')
    print('MESH_STY', mesh_large_df_sty['Preferred Label'].nunique(), 'mesh')

print('MESH_STY_LITVAR1', len(mesh_litvar_sty['Preferred Label']), ' rows')
print('MESH_STY_LITVAR1', mesh_litvar_sty['Preferred Label'].nunique(), 'mesh')
mesh_litvar_sty[['Preferred Label', 'Semantic Types Label', 'Class ID', 'mesh_id',
                 'Semantic Types']]#.columns


# In[ ]:


## ADD mesh_id col
if 'mesh_id' not in mesh_litvar_sty.columns:
    mesh_litvar_sty['mesh_id'] = mesh_litvar_sty['Class ID'].str.replace('http://purl.bioontology.org/ontology/MESH/', '')
    mesh_litvar_sty['mesh_id']
    mesh_litvar_sty.columns
    new_order = ['Preferred Label', 'Semantic Types Label', 'Class ID', 'mesh_id', 'Semantic Types',
                 'Synonyms', 'Parents', 'CUI', 'AQL', 'TERMUI']
    mesh_litvar_sty = mesh_litvar_sty[new_order]
    mesh_litvar_sty.to_csv('ref-mesh-archive/MESH_STY_LITVAR1.csv')
else:
    print('mesh_id aready added to csv')


# ## add mesh_id to grpm db

# In[ ]:


grpm_b_df = pd.read_csv('grpm_db_pcg/grpm_table_output.csv',index_col=0)
grpm_b_df


# In[ ]:


mesh_id_ref_df = mesh_litvar_sty[['Preferred Label','mesh_id']]
mesh_id_ref_df
grpm_db_merge_id = pd.merge(grpm_b_df, mesh_id_ref_df, left_on='mesh', right_on='Preferred Label')


# In[ ]:


#grpm_db_merge_id[['mesh','Preferred Label','mesh_id']].drop_duplicates()
grpm_db_merge_id = grpm_db_merge_id.drop('Preferred Label', axis=1)


# In[ ]:


grpm_db_merge_id.columns
new_order = ['gene', 'rsid', 'pmids', 'mesh_id']
grpm_db_merge_id = grpm_db_merge_id[new_order].drop_duplicates()
grpm_db_merge_id


# In[ ]:


grpm_db_merge_id.to_csv('grpm_db_pcg/grpm_table_output_id.csv')


# In[ ]:


check = pd.read_csv('grpm_db_pcg/grpm_table_output_id.csv', index_col=0)
check


# ### Sorting Mesh-STY

# In[ ]:


if 'mesh_large_df_sty' in globals():
    #modulo groupby and bar
    mesh_large_df_sty_less = mesh_large_df_sty[['Preferred Label','Semantic Types Label']]

    ### groupby.describe analysis by rsid--------------------
    mesh_large_df_sty_less_count = mesh_large_df_sty_less.groupby('Semantic Types Label').describe().reset_index()
    mesh_large_df_sty_less_count.columns = mesh_large_df_sty_less_count.columns.to_flat_index()
    new_column_names = ['Semantic Types Label', 'mesh-count', 'mesh-unique','mesh-top','mesh-freq']
    mesh_large_df_sty_less_count.columns = new_column_names

    mesh_large_df_sty_less_count_sort = mesh_large_df_sty_less_count.sort_values(by='mesh-count',ascending=False).reset_index(drop=True)
    print('MESH_STY')
    print(len(mesh_large_df_sty_less_count_sort))

mesh_large_df_sty_less_count_sort[['mesh-count','Semantic Types Label']]#.to_clipboard()


# In[ ]:


#Graph Barh
num = 20
x = mesh_large_df_sty_less_count_sort['Semantic Types Label'].iloc[:num]
y = mesh_large_df_sty_less_count_sort['mesh-count'].iloc[:num]
plt.figure(figsize=(4, len(x) *0.25))
plt.title('Global Mesh- Semantic Type enrichment', loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('Semantic Types')
plt.xlabel('mesh', position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
# use log scale
plt.gca().set_xscale('log')
#plt.savefig('Reference Mesh- Semantic Type enrichment.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


#'Global "Mesh-Semantic Type" Zipf s law'
x = mesh_large_df_sty_less_count_sort['Semantic Types Label'].iloc[:]
y = mesh_large_df_sty_less_count_sort['mesh-count'].sort_values().iloc[:]
plt.figure(figsize=(6, 5))
plt.title('Global "Mesh-Semantic Type" Zipf s law', loc='center',pad=10)

plt.plot(x,y)
plt.gca().invert_yaxis()

plt.ylabel('mesh number')
plt.xlabel('semantic type rank', position=(0.5, 1.08))

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.show()


# In[ ]:


## Density Rugged Plot
data = mesh_large_df_sty_less_count_sort['mesh-count'].iloc[:]

sns.set_style('whitegrid')
sns.kdeplot(np.array(data), bw_method=0.5)
sns.rugplot(np.array(data), color='r')

#plt.yscale('log')
plt.title('Mesh abundance for Semantic Type')
#plt.yscale('log')
print('Mesh abundance for Semantic Type:')
plt.show()


# 

# ## Analyze MESH-STY-LITVAR (subset of MESH-STY.csv)

# In[ ]:


print('mesh_litvar_sty mesh:',mesh_litvar_sty['Preferred Label'].nunique())

memory = mesh_litvar_sty.memory_usage().sum()
print(f'The memory_usage of mesh_litvar_sty is {memory/ (1024 * 1024):.2f} MB.')

file_size = os.path.getsize('ref-mesh-archive/MESH_STY_LITVAR1.csv')
print(f'The size of MESH_STY_LITVAR1.csv is {file_size/ (1024 * 1024):.2f} MB.')


# #### Sort Mesh-Sty-Litvar

# In[ ]:


#modulo groupby and bar
mesh_litvar_sty_less = mesh_litvar_sty[['Preferred Label','Semantic Types Label']]

### groupby.describe analysis by rsid--------------------
mesh_litvar_sty_less_count = mesh_litvar_sty_less.groupby('Semantic Types Label').describe().reset_index()
mesh_litvar_sty_less_count.columns = mesh_litvar_sty_less_count.columns.to_flat_index()
new_column_names = ['Semantic Types Label', 'mesh-count', 'mesh-unique','mesh-top','mesh-freq']
mesh_litvar_sty_less_count.columns = new_column_names

mesh_litvar_sty_less_count_sort = mesh_litvar_sty_less_count.sort_values(by='mesh-count',ascending=False).reset_index(drop=True)
print('MESH_STY_LITVAR1')
mesh_litvar_sty_less_count_sort[['mesh-count','Semantic Types Label']]#.to_clipboard()


# In[ ]:


#Graph Barh
num = 20
x = mesh_litvar_sty_less_count_sort['Semantic Types Label'].iloc[:num]
y = mesh_litvar_sty_less_count_sort['mesh-count'].iloc[:num]
plt.figure(figsize=(4, len(x)*0.25))
plt.title('Litvar1 Mesh- Semantic Type enrichment', loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('Semantic Types')
plt.xlabel('mesh', position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
# use log scale
plt.gca().set_xscale('log')
#plt.savefig('Reference Mesh- Semantic Type enrichment.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


#'Global "Mesh-Semantic Type" Zipf s law'
x = mesh_litvar_sty_less_count_sort['Semantic Types Label'].iloc[:]
y = mesh_litvar_sty_less_count_sort['mesh-count'].sort_values().iloc[:]
plt.figure(figsize=(6, 5))
plt.title('LitVar "Mesh-Semantic Type" Zipf s law', loc='center',pad=10)

plt.plot(x,y)
plt.gca().invert_yaxis()

plt.ylabel('mesh number')
plt.xlabel('semantic type rank', position=(0.5, 1.08))

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.show()


# In[ ]:


# VISUALIZE DIFFERENCES

#reload sorted df
if 'mesh_large_df_sty_less_count_sort' in locals() and isinstance(mesh_large_df_sty_less_count_sort, pd.DataFrame):
    pass
else:
    mesh_large_df_sty_less = mesh_large_df_sty[['Preferred Label','Semantic Types Label']]
    mesh_large_df_sty_less_count = mesh_large_df_sty_less.groupby('Semantic Types Label').describe().reset_index()
    mesh_large_df_sty_less_count.columns = mesh_large_df_sty_less_count.columns.to_flat_index()
    new_column_names = ['Semantic Types Label', 'mesh-count', 'mesh-unique','mesh-top','mesh-freq']
    mesh_large_df_sty_less_count.columns = new_column_names
    mesh_large_df_sty_less_count_sort = mesh_large_df_sty_less_count.sort_values(by='mesh-count',ascending=False).reset_index(drop=True)
if 'mesh_litvar_sty_less_count_sort' in locals() and isinstance(mesh_litvar_sty_less_count_sort, pd.DataFrame):
    pass
else:
    mesh_litvar_sty_less = mesh_litvar_sty[['Preferred Label','Semantic Types Label']]
    mesh_litvar_sty_less_count = mesh_litvar_sty_less.groupby('Semantic Types Label').describe().reset_index()
    mesh_litvar_sty_less_count.columns = mesh_litvar_sty_less_count.columns.to_flat_index()
    new_column_names = ['Semantic Types Label', 'mesh-count', 'mesh-unique','mesh-top','mesh-freq']
    mesh_litvar_sty_less_count.columns = new_column_names
    mesh_litvar_sty_less_count_sort = mesh_litvar_sty_less_count.sort_values(by='mesh-count',ascending=False).reset_index(drop=True)

# create two sample dataframes
df1 = mesh_large_df_sty_less_count_sort[['mesh-count', 'Semantic Types Label']]
df2 = mesh_litvar_sty_less_count_sort[['mesh-count', 'Semantic Types Label']]

# add a 'source' column to each dataframe
df1['source'] = 'global'
df2['source'] = 'litvar'

# combine the dataframes
combined_df = pd.merge(df1, df2, on=['Semantic Types Label'], how='outer', suffixes=('_df1', '_df2'))

# sort the dataframe by column 'A'
#combined_df = combined_df.sort_values('mesh-count')

# reset the index
combined_df = combined_df.reset_index(drop=True)

# display the combined dataframe
#combined_df['mesh-count_df2'] = combined_df['mesh-count_df2'].replace(np.nan, 0)
combined_df = combined_df[-(combined_df['Semantic Types Label'] == 'Drug Delivery Device')]
combined_df
mesh_large_df_sty_less['Preferred Label'].nunique()
combined_df = combined_df.sort_values(by='mesh-count_df2', ascending=False).reset_index(drop=True)


# define a formatting function that generates a proportional bar
def format_bar(value):
    max_value = combined_df['mesh-count_df1'].max().max()  # get the maximum value in the dataframe
    bar_width = int(value / max_value * 100)  # calculate the width of the bar as a percentage
    return f'<div style="background-color: blue; width: {bar_width}%">{value}</div>'

def format_bar_2(value):
    max_value = combined_df['mesh-count_df2'].max().max()  # get the maximum value in the dataframe
    bar_width = int(value / max_value * 100)  # calculate the width of the bar as a percentage
    return f'<div style="background-color: blue; width: {bar_width}%">{value}</div>'


# apply the formatting function to the dataframe
df_formatted = combined_df.style.format({'mesh-count_df1': format_bar, 'mesh-count_df2': format_bar_2})

# save the formatted dataframe to an HTML file
with open('formatted_dataframe.html', 'w') as f:
    f.write(df_formatted.render())

# display the formatted dataframe
df_formatted


# 

# ####Define query variables

# df.iloc[6]
# #df.loc['diet therapy']
# #find specific variable
# var1 = 'diet therapy'
# var2 = 'Diet, Food, and Nutrition'
# var3 = 'D004035'
# var4 = 'D000066888' #'Diet, Food, and Nutrition'
# var5 = 'Beverage'

# x = df[df.eq(var1).any(1)]
# y = df[df.eq(var2).any(1)]
# z = df[df.eq(var3).any(1)]
# k = df[df.eq(var4).any(1)]
# w = df[df.eq(var5).any(1)]
# 
# list = [x,y,z]

# In[ ]:





# # Reference Mesh list build:

# Build a 370 coherent and omic list of mesh term related to different topics
# - neurodegenerative diseases
# - skin diseases
# - infective diseases
# - reproductive phisyology
# - cancer

# ## Check avalable refs:

# In[ ]:


#Check avalable refs:
folder_path = "ref-mesh-archive"  # Replace with the actual folder path

# Create a file path pattern to match CSV files
file_pattern = os.path.join(folder_path, "*.csv")
file_pattern

# Use glob to get a list of file paths matching the pattern
csv_files = glob.glob(file_pattern)
csv_files_name = []
# Print the list of CSV files
for file in csv_files:
    file_name = os.path.basename(file)
    csv_files_name.append(file_name)

print('Available reference mesh lists:')
csv_files_df = pd.Series(csv_files_name)
ref_files_df = csv_files_df[csv_files_df.str.contains('ref_mesh_ng_')].reset_index(drop=True)
ref_files_df_small = csv_files_df[csv_files_df.str.contains('small_ref_mesh_ng_')].reset_index(drop=True)
ref_files_df_small


# In[ ]:


ref_names =  [ 'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_',
               'ref_' ]
ref_names = [string + str(i) for i, string in enumerate(ref_names)]

ref_list = []
for name, ref, save  in zip(ref_names, ref_files_df, ref_files_df):
    df = pd.read_csv('ref-mesh-archive/'+ref, index_col=0).reset_index(drop=True)
    df = df.drop('Semantic Types Label', axis=1).drop_duplicates()
    df = df.sort_values(by='mesh_id',ascending=True)
    globals()[name] = df
    df.to_csv('ref-mesh-archive/small_'+save)
    ref_list.append(df)


# In[ ]:


ref_0


# ### create JSON reference
# https://jsoneditoronline.org/#left=local.lejobi&right=local.lejobi
# https://www.convertjson.com/json-to-html-table.htm

# In[ ]:


import json
import csv

folder_path = "ref-mesh-archive"


num_elements = len(ref_files_df)  # Number of elements you want in the list
value_list = []
for i in range(1, num_elements + 1):
    element = "file_pat" + str(i)
    value_list.append(element)

jdata_list_n = ['jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_',
              'jdata_' ]
jdata_list_n = [string + str(i) for i, string in enumerate(jdata_list_n)]

jdata_list = []
for name in jdata_list_n:
    empty =[]
    globals()[name] = empty
    jdata_list.append(empty)


file_patterns = []
choose_ref = ref_files_df_small # or ref_files_df
for i in choose_ref:
    file_patterns.append(os.path.join(folder_path, i))

for i, j  in zip(range(len(ref_files_df)), jdata_list):

    with open(file_patterns[i], 'r') as file1:
        csv_reader = csv.DictReader(file1)
        # for row in csv_reader:
        #     j_data1.append(row)
        for row in csv_reader:
            modified_row = {key: value for key, value in row.items() if key != ''} #drop col
            j.append(modified_row)


# In[ ]:


#combined_data = {
#    'data1': j_data1,
#    'data2': j_data2 }

combined_data = {}

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
    "Xenobiotics Metabolism",
]

for tit, j_data in zip(titles, jdata_list):
    key = tit
    combined_data[key] = j_data

#for i, j_data in enumerate(jdata_list, start=1):
#    key = 'data' + str(i)
#    combined_data[key] = j_data

print(combined_data)
#pd.DataFrame(combined_data)
import pyperclip
pyperclip.copy(str(combined_data))

json_file = 'nutrigenetis_referencemesh_small.json'
with open(json_file, 'w') as outfile:
    json.dump(combined_data, outfile, indent=4)


# ### add mesh id to every mesh list

# In[ ]:


# add mesh id to every mesh list
label = 'ng_nutri'
dff = pd.read_csv('ref-mesh-archive/ref_mesh_'+label+'.csv',index_col=0).reset_index(drop=True)#.drop('Unnamed: 0', axis=1)

dff_id = pd.merge(dff, mesh_id_ref_df, left_on='Preferred Label', right_on='Preferred Label')
dff_id


# In[ ]:


dff_id.to_csv('ref-mesh-archive/ref_mesh_'+label+'.csv')
#grpm_db_merge_id = grpm_db_merge_id.drop('Preferred Label', axis=1)
check = pd.read_csv('ref-mesh-archive/ref_mesh_'+label+'.csv', index_col=0)
check


# In[ ]:


chack = check.drop('mesh_id_y',axis=1)
check = check.rename(columns={'mesh_id_y': 'mesh_id'})
check.to_csv('ref-mesh-archive/ref_mesh_'+label+'.csv')


# ### ref mesh cleaner

# In[ ]:


#Get rid from ref_mesh and run again
xeno = ['Polymorphism, Genetic', 'DNA Methylation','Liver Neoplasms','Energy Metabolism','Metabolism, Inborn Errors','X-Ray Absorption Spectroscopy','Blood-Testis Barrier']
oxi_stress = ['Food Deprivation','Food Analysis','Food, Formulated','Food Microbiology','Food Intolerance','Food Quality','Food Industry','Foods, Specialized','Food Chain','Foods, Specialized','Food Chain']
intol = ['Depression','Osteoporosis','Immunotherapy','Anxiety','Anti-Inflammatory Agents, Non-Steroidal','Immunotherapy, Adoptive','Desensitization, Immunologic','Peanut Hypersensitivity','Milk Hypersensitivity', 'Egg Hypersensitivity','Immunotherapy, Active','Neurologic Manifestations' ,'Infectious Anemia Virus, Equine','Nut and Peanut Hypersensitivity' ,'Diarrhea Virus 1, Bovine Viral','Eye Movement Desensitization Reprocessing','Diarrhea Virus 2, Bovine Viral']
eat_taste = ['Blood Glucose','Mouth Neoplasms','Glucose Transporter Type 1','Glucose Transporter Type 2','Auditory Perception','Citric Acid Cycle','Glucose Transporter Type 4','Loeys-Dietz Syndrome','United States Food and Drug Administration','Sodium-Glucose Transporter 2','Sodium-Glucose Transporter 2 Inhibitors','Glucose-6-Phosphate','Glucose Transporter Type 3','Pitch Perception','Depth Perception','Glucose-1-Phosphate Adenylyltransferase','Glucose Dehydrogenases','Glucosephosphates']
cvd = ['DNA, Mitochondrial','Cell Differentiation','Protein Transport','Mitochondrial Proteins','Toll-Like Receptors','Genome, Mitochondrial','Genes, Mitochondrial','Mitochondrial Precursor Protein Import Complex Proteins','Mitochondrial Precursor Protein Import Complex Proteins','Mitochondrial Proton-Translocating ATPases','Mitochondrial Dynamics','Mitochondrial Membranes','Leukemia, Plasma Cell','Heart Neoplasms']
dmt2_ms = [ 'MicroRNAs', 'Mitochondria', 'DNA, Mitochondrial', 'Mitochondrial Proteins', 'Pancreatic Neoplasms','Autophagy','Mitochondrial Diseases','Genome, Mitochondrial','Genes, Mitochondrial','Mitochondrial Precursor Protein Import Complex Proteins','Mitochondrial Precursor Protein Import Complex Proteins','Mitochondrial Dynamics','Mitochondrial Membranes','Mitochondrial Myopathies','RNA, Mitochondrial', 'Mitochondrial Encephalomyopathies','AIDS-Associated Nephropathy','Familial Primary Pulmonary Hypertension','Mitochondria, Heart','Autophagy-Related Protein-1 Homolog','Ocular Hypertension','Autophagy-Related Protein 7','Mitochondrial Trifunctional Protein','Mitochondrial Permeability Transition Pore','Mitochondrial Uncoupling Proteins','Mitochondrial ADP, ATP Translocases','Autophagy-Related Protein 8 Family','Mitochondrial Ribosomes','Mitochondrial Trifunctional Protein, alpha Subunit','Autophagy-Related Protein 12','Mitochondrial Turnover','Lysergic Acid Diethylamide','Mitochondrial Trifunctional Protein, beta Subunit','Chaperone-Mediated Autophagy','Mitochondrial Swelling','Mitochondrial Transmembrane Permeability-Driven Necrosis','Mitochondrial Size']


# In[ ]:


#check mesh ref

label = 'xeno'
dff = pd.read_csv('ref-mesh-archive/ref_mesh_'+label+'.csv', index_col=0)
dff['Preferred Label'].drop_duplicates()

# get rid exact
save_clean = True

#clean mesh ref
get_rid_list = xeno

dff_rid = dff.copy()
mask = dff['Preferred Label'].isin(get_rid_list)
dff_rid = dff_rid[-mask]
print(dff.describe())
print(dff_rid.describe())


# In[ ]:


if save_clean == True:
    dff_rid.to_csv('ref-mesh-archive/ref_mesh_'+label+'.csv')

dff_rid['Preferred Label'].drop_duplicates()


# 

# dff['Preferred Label'].drop_duplicates()
# check = 'Fistula'
# dff_check = dff[dff['Preferred Label'].str.contains(check)]
# dff_check
# dff['Preferred Label'].drop_duplicates()

# In[ ]:


label = 'ob_bmi_less'
dff = pd.read_csv('ref-mesh-archive/ref_mesh_'+label+'.csv', index_col=0)
dff['Preferred Label'].drop_duplicates()


# get rid contains
save_clean = True
tag = ''

#clean mesh ref
get_rid_list = ['Infant','Child','Adolescent','Elder','Maternal','Youth','Man','Woman','National','Neoplasm''Epidemiolo','Reproductive','Sexual','Genome','Animal','Doping','Social','Urban''Health C','Health E','Health F','Health I','Health P','Health S','Health T','ty Health','le Health','ic Health'
]
get_rid_list = ['Polymor','Genetic','Risk F']##['Fistula', 'Neoplasm','Labor', 'Chick']

dff_rid = dff.copy()
for get_rid in get_rid_list:
    mask = dff['Preferred Label'].str.contains(get_rid)
    dff_rid = dff_rid[-mask]

if save_clean == True:
    dff_rid.to_csv('ref-mesh-archive/ref_mesh_'+label+'_'+tag+'.csv')

dff_rid['Preferred Label'].drop_duplicates()


# In[ ]:


dff['Preferred Label'].drop_duplicates()


# ## Build new ref_mesh lists

# In[ ]:


mesh_litvar_sty


# ### Split nutritional topic into 8 different interest categories

# In[ ]:


ngdb_df = pd.read_csv('NGDB1.csv', encoding='utf-8', sep=';')
print('ngbd statistics:')
print('caterories', ngdb_df.Category.nunique())
print('genes', ngdb_df[' Gene'].nunique())
print('SNPs', ngdb_df.SNPId.nunique())
print('consequence type', ngdb_df['Gene_Consequence '].nunique())
ngdb_df.Category.drop_duplicates()


# 

# ### General Ref-Mesh from ChatGPT input
# Per lista generata da ChatGPT: cerca in Preferred Labels and Synonyms corrisponding mesh entry
#     Two general classes of mesh of interest:
#     - 'topic' physiology
#     - 'topic' pathology
#     To be concatenated or searched separately?

# In[ ]:


### Get ChatGPT - API

# https://pub.towardsai.net/how-to-use-chatgpt-api-for-direct-interaction-from-colab-or-databricks-39969a0ead5f
# https://platform.openai.com/account/api-keys
# https://platform.openai.com/docs/api-reference/introduction


openai.api_key = "your api key"
# please-paste-your-API-key-here

def responseGPT(prompt):
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "user", "content": prompt}])
    return print(completion)

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

def createImageWithGPT(prompt):
    completion = openai.Image.create(
        prompt=prompt,
        n=1,
        size="512x512"
    )
    return IPython.display.HTML("<img src =" + completion.data[0].url + ">")


# In[ ]:


nutritional_topic = ['diseases and disorders realted to nutrition and diet ', 'diet, food consuption, eating behaviour and nutrition']
infective_topic = ['infective agents, bacteria, virus and protozoan','infective diseases']
reproductive_topic = ['reproductive system physiology','reproductive system pathology', 'Assisted reproductive technology']

nutritional_topics = [
    ['Obesity, overweight and body weight control', 'compulsive eating behavior'],
    ['cardiovascular diseases','physiological processes realted to cardiovascular diseases','lipid metabolism in the context of cardiovascular diseases'],
    ['Diabetes Melitus Type II and metabolic syndrome'],
    ['Vitamin metabolism and Vitamins recommended intake levels','Micronutrients metabolism and Micronutrient recommended intake levels', 'disease related to vitamins and micronutrients deficiency'],
    ['eating behaviour and taste sensation'],
    ['food intolerances'],
    ['food allergies'],
    ['diet-induced oxidative stress'],
    ['metabolism of xenobiotics'],
]
pd.Series(nutritional_topics)
pd.Series(nutritional_topics).to_clipboard()


# In[ ]:


# GPT prompts

# parameters--------------------------------------
import pyperclip
num_mesh = 120
topics = nutritional_topics
topic_id = 7
#-----------------------------------------------
topic_01  = topics[topic_id][0]
topic_02  = topics[topic_id][1] if len(topics[topic_id])>=2 else None
topic_03 = topics[topic_id][2] if len(topics[topic_id])>=3 else None

prompt_01 = "give me a comprehensive pyhton list of "+str(num_mesh)+" real Pubmed Mesh terms related to "+ topic_01 +":\n gpt_01 = ['A','B','C',...]\n"
prompt_02 = "give me a comprehensive pyhton list of "+str(num_mesh)+" real Pubmed Mesh terms related to "+ topic_02 +":\n gpt_02 = ['A','B','C',...]\n" if len(topics[topic_id])>=2 else None
prompt_03 = "give me a comprehensive pyhton list of "+str(num_mesh)+" real Pubmed Mesh terms related to "+ topic_03 +":\n gpt_03 = ['A','B','C',...]" if len(topics[topic_id])>=3 else None

pyperclip.copy(prompt_01)
pyperclip.copy(prompt_01+prompt_02) if len(topics[topic_id])>=2 else None
pyperclip.copy(prompt_01+prompt_02+prompt_03) if len(topics[topic_id])>=3 else None


# In[ ]:


def getGPTterms(topic):
    completion = openai.ChatCompletion.create(
        model="gpt-3.5-turbo",
        messages=[
            {"role": "user", "content": "give me a comprehensive list of "+str(num_mesh)+" real Pubmed Mesh terms related to "+ topic +" in pathon list format (list = ['A','B','C',...]):\n gpt_01 = ['A','B','C',...]\n"}
        ]
    )
    return completion#.choices[0].message.content)

response = getGPTterms(topic_01)


# In[ ]:


print(response.choices[0].message.content)


# ### use GPT API
# persona = "Socrates"
# instruction = "**instruction: you are "+ persona +", answer like you are "+ persona +"** \n"
# 
# question = "Tell me about the four elements?"
# prompt = instruction +  question
# chatWithGPT(prompt)

# In[ ]:





# ### get gpt terms

# In[ ]:


folder_path = "ref-mesh-archive/gpt_terms"
file_pattern = os.path.join(folder_path, "*.csv")
csv_files = glob.glob(file_pattern)
csv_files_name = []
for file in csv_files:
    file_name = os.path.basename(file)
    csv_files_name.append(file_name)

print('Available gpt terms lists:')
csv_files_df = pd.Series(csv_files_name)
csv_files_df[csv_files_df.str.contains('gpt_')]


# In[ ]:


# Set parameters:--------------------------------
topic_label = 'xeno_TRIAL' # <== CHANGE THIS! Overwrite risk!
mother_df = mesh_litvar_sty

#load existing gpt terms:
overwrite_gpt = True
load_saved_gpt = False
topic_label = 'cvd' if load_saved_gpt == True else topic_label
#-----------------------------------------------

#ChatGPT list:
gpt_01 = ['Xenobiotics',
 'Metabolism',
 'Biotransformation',
 'Cytochrome P-450 Enzyme System',
 'Phase I Metabolic Detoxication',
 'Phase II Metabolic Detoxication',
 'Oxidation-Reduction',
 'Hydroxylation',
 'Dealkylation',
 'Deamination',
 'Glucuronidation',
          'Conjugation',
          'Amination',
          'Glycosylation',
 'Sulfation',
 'Methylation',
 'Acetylation',
 'Glutathione Transferase',
 'Conjugation, Genetic',
 'Pharmacokinetics',
 'Absorption',
 'Distribution',
 'Excretion',
 'Half-Life',
 'Clearance',
 'Bioavailability',
 'Blood-Brain Barrier',
 'Blood-Air Barrier',
 'Blood-Retinal Barrier',
 'Blood-Testis Barrier',
 'Placental Barrier',
 'Biological Transport',
 'Active Transport, Cell Nucleus',
 'Facilitated Diffusion',
 'Endocytosis',
 'Exocytosis',
 'Transcytosis',
          "Xenobiotics",
          "Drug Metabolism",
          "Cytochrome P-450 Enzyme System",
          "Phase I Metabolic Detoxication",
          "Phase II Metabolic Detoxication",
          "Toxicokinetics",
          "Biotransformation",
          "Glutathione Transferase",
          "Oxidoreductases",
          "Drug-Related Side Effects and Adverse Reactions",
          "Enzyme Induction",
          "Pharmacokinetics",
          "Enzyme Inhibitors",
          "Chemical and Drug Induced Liver Injury",
          "Drug Interactions",
          "Liver",
          "Toxicity Tests",
          "Pharmaceutical Preparations",
          "Metabolomics",
          "Drug Evaluation, Preclinical",
          "Liver Microsomes",
          "Liver Neoplasms, Experimental",
          "Drug Discovery",
          "Polymorphism, Genetic",
          "Drug Resistance",
          "Toxicity",
          "Hepatocytes",
          "Carboxylesterase",
          "Aryl Hydrocarbon Hydroxylases",
          "Mixed Function Oxygenases",
          "Uridine Diphosphate Glucuronic Acid",
          "Pharmacogenetics",
          "Drug-Induced Liver Injury",
          "Pharmaceutical Preparations",
          "Glutathione",
          "Metabolome",
          "Toxicokinetics",
          "Receptors, Aryl Hydrocarbon",
          "P-Glycoprotein",
          "Hepatocytes",
          "Drug Resistance",
          "Drug Interactions",
          "Cytochrome P-450",
          "CYP1A2",
          "Drug Metabolism, Phase I",
          "Drug Metabolism, Phase II",
          "Oxidative Stress",
          "Pharmacokinetics",
          "Uridine Diphosphate Glucuronic Acid",
          "Liver",
          "Enzyme Induction",
          "Drug-Related Side Effects and Adverse Reactions",
          "Pharmacogenetics",
          "Liver Diseases",
          "Pharmacogenomics",
          "Receptors, Drug",
          "Liver Neoplasms, Experimental",
          "Liver Microsomes",
          "Drug Evaluation, Preclinical",
          "Drug Resistance",
          "Liver Cirrhosis, Experimental",
          "Liver Failure",
          "Liver Function Tests",
          "Liver Regeneration",
          "Toxicology",
          "Toxicity Tests",
          "Drug Design",
          "Drug Screening Assays, Antitumor",
          "Liver Diseases, Alcoholic",
          "Liver Transplantation",
          "Toxicity",
          "Toxicology",
          "Pharmacy Services",
          "Polypharmacy",
          "Toxicology",
          "Toxicogenetics",
          "Toxicology",
          "Toxicology, Environmental",
          "Toxicology, Forensic",
          "Toxicology, Legal",
          "Toxicology, Pharmaceutical",
          "Toxicology, Veterinary",
          "Toxicology Tests",
          "Toxicology, Clinical",
          "Toxicity Tests",
          "Toxicity Tests, Acute",
          "Toxicity Tests, Chronic",
          "Toxicity Tests, Subacute",
          "Toxicity Tests, Subchronic",
          "Toxicity Tests, Sublethal",
          "Toxicity Tests, Mutagenicity",
          "Toxicity Tests, Reproduction",
          "Toxicity Tests, Teratogenicity",
          "Xenobiotics",
          "Environmental Pollutants",
          "Water Pollutants, Chemical",
          "Soil Pollutants",
          "Air Pollutants",
          "Pesticides",
          "Pharmaceutical Preparations",
          "Cosmetics",
          "Food Additives",
          "Radiation, Ionizing",
          "Radiation, Nonionizing",
          "Toxicity",
          "Toxicology",
          "Toxicology, Environmental",
          "Toxicology, Forensic",
          "Toxicology, Legal",
          "Toxicology, Pharmaceutical",
          "Toxicology, Veterinary",
          "Pharmacology",
          "Pharmacology, Clinical",
          "Pharmacology, Experimental",
          "Pharmacy",
          "Pharmacy Service, Hospital",
          "Pharmacy Services",
          "Polypharmacy",
          "Toxicogenetics",
          "Drug-Related Side Effects and Adverse Reactions",
          "Chemical and Drug Induced Liver Injury",
          "Liver Diseases",
          "Liver Cirrhosis, Experimental",
          "Liver Failure",
          "Liver Function Tests",
          "Liver Neoplasms, Experimental",
          "Liver Regeneration",
          "Liver Transplantation",
          "Carboxylesterase",
          "Aryl Hydrocarbon Hydroxylases",
          "Mixed Function Oxygenases",
          "Glutathione Transferase",
          "Receptors, Aryl Hydrocarbon",
          "P-Glycoprotein",
          "Cytochrome P-450 CYP1A2",
          "Oxidative Stress",
          "Drug Metabolism, Phase I",
          "Drug Metabolism, Phase II",
          "Pharmacokinetics",
          "Pharmacogenetics",
          "Pharmacogenomics"
          ]
#gpt_01 = []
gpt_02 =  []
#gpt_02 = []
gpt_03  = []# Paste here
#n.b. remember to remove  general terms like "Proteins" to avoid bias (get_rid section)

gpt_all_ser = pd.DataFrame()
if load_saved_gpt == True:
    if os.path.exists('ref-mesh-archive/gpt_terms/gpt_terms_'+topic_label+'.csv'):
        gpt_all_ser_df = pd.read_csv('ref-mesh-archive/gpt_terms_'+topic_label+'.csv')
        gpt_all_ser = gpt_all_ser_df['gpt_terms'].drop_duplicates()
        print('gpt terms loaded')
    else:
        print('no gpt terms to load')
else:
    gpt_01_df = pd.DataFrame({'gpt_terms': gpt_01})
    gpt_01_df['series'] = 'gpt_01'
    gpt_all_ser_df = gpt_01_df
    if len(gpt_02) > 0 :
        print('gpt_02 included')
        gpt_02_df = pd.DataFrame({'gpt_terms': gpt_02})
        gpt_02_df['series'] = 'gpt_02'
        gpt_all_ser_df = pd.concat([gpt_all_ser_df, gpt_02_df], axis=0)
    if len(gpt_03) > 0 :
        print('gpt_03 included')
        gpt_03_df = pd.DataFrame({'gpt_terms': gpt_03})
        gpt_03_df['series'] = 'gpt_03'
        gpt_all_ser_df = pd.concat([gpt_all_ser_df, gpt_03_df], axis=0)

#get rid of unwanted gpt terms: Exact match
print('raw', len(gpt_all_ser_df))
get_rid = ['Proteins','Protein','DNA','RNA',
           'Polymor','Genetic',
           'Infant',
           'Child',
           'Adolescent',
           'Elder',
           'Maternal',
           'Youth',
           'Man',
           'Woman',
           'National',
           'Neoplasm',
           'Reproductive',
           'Sexual',
           'Genome',
           'Animal',
           'Doping',
           'Social',
           'Urban',
           ] # generic terms creating bias
# attetion to 'Mitochondrial'
mask = gpt_all_ser_df['gpt_terms'].isin(get_rid)
gpt_all_ser_df = gpt_all_ser_df[-mask]
print('clean I:',len(gpt_all_ser_df))

#get rid of unwanted gpt terms: Contained str
get_rid_gpt = True
get_rid_spec = ['Pharma',
                'Cosme',
                'Drug',
                'Distribution'
                ]
if get_rid_gpt == True:
    for rid in get_rid_spec:
        mask = gpt_all_ser_df['gpt_terms'].str.contains(rid)#.isin(get_rid_spec)
        gpt_all_ser_df = gpt_all_ser_df[-mask]
    print('clean II:',len(gpt_all_ser_df))

print('drop dup:',gpt_all_ser_df.gpt_terms.nunique())

if overwrite_gpt == True:
    gpt_all_ser_df.to_csv('ref-mesh-archive/gpt_terms/gpt_terms_'+topic_label+'.csv')
    print('gpt terms saved')

gpt_all_ser = gpt_all_ser_df['gpt_terms'].drop_duplicates()
pyperclip.copy(str(gpt_all_ser.to_list()))
gpt_all_ser_df
#pyperclip.copy(str(gpt_all_ser.to_list()))


# Get rid of gpt terms LOG
# on nutri: [ 'Immune System',    'Nail Concentrations',    'Natural Killer Cells',    'Nerve Function',    'PYY3-36',    'Phenotype',    'Pregnancy',    'Single Nucleotide Polymorphisms',    'Transplantation',    'T-Lymphocytes',    'Dendritic Cells',    'Macrophages',    'Graft',    'Phagocytosis',    'Calcium',    'Iron',    'Phosphorus',    'Copper',    'Development',    'Developmental',    'Gene Expression']
# on ob_bmi: ['Transcriptome',
# 'RNA, Ribosomal, 16S','Proteomics','Prevalence','Epidemiologic Studies','Epidemiology','Epidemics','Sick Role']
# oxi_stress. ['Chemokines','Cytokines','Diet','Interleukin']
# xeno: ['Pharma','Cosme','Drug','Distribution']

# ### cross gpt terms with litvar_mesh

# In[ ]:


# Use gpt terms to filter mother_df and Get Mesh:
overwrite_ref = True

split_01_02 = False
if split_01_02 == True:
    print(topics[topic_id][0],'topic')
    #get mesh
    chat_mesh_01 = pd.DataFrame(columns=mother_df.columns) #Empty DataFrame, same structure
    for i in gpt_01:
        meshy = mother_df[mother_df['Preferred Label'].str.contains(i).fillna(False)]
        chat_mesh_01 = pd.concat([meshy,chat_mesh_01])
    print(chat_mesh_01['Preferred Label'].nunique())

    print(topics[topic_id][1],'topic')
    #get mesh
    chat_mesh_02 = pd.DataFrame(columns=mother_df.columns) #Empty DataFrame, same structure
    for i in gpt_02:
        meshy = mother_df[mother_df['Preferred Label'].str.contains(i).fillna(False)]
        chat_mesh_02 = pd.concat([meshy,chat_mesh_02])
    print(chat_mesh_02['Preferred Label'].nunique())
    if 1 >2 :
        #analyze semantics
        chat_mesh_01_grup = chat_mesh_01.groupby('Semantic Types Label').describe()
        chat_mesh_01_sem = chat_mesh_01_grup.index.to_list()
        chat_mesh_02_grup = chat_mesh_02.groupby('Semantic Types Label').describe()
        chat_mesh_02_sem = chat_mesh_02_grup.index.to_list()
        #get sum of semantic types
        chat_mesh_sem_sum = chat_mesh_02_sem + chat_mesh_01_sem
        chat_mesh_sem_sum_ser = pd.Series(chat_mesh_sem_sum)
        print('chat_mesh_sem_sum_ser',chat_mesh_sem_sum_ser.nunique())

#get mesh
chat_mesh_all = pd.DataFrame(columns=mother_df.columns) #Empty DataFrame, same structure
time1 =datetime.now()
for i in gpt_all_ser:
    #meshy = mother_df[mother_df['Preferred Label'].str.contains(i).fillna(False)]
    meshy = mother_df[mother_df['Preferred Label'].str.extract(f'({i})').notna().any(axis=1)]
    chat_mesh_all = pd.concat([meshy,chat_mesh_all])
time2 =datetime.now()
print(time2-time1)

#analyze semantics
chat_mesh_all_grup = chat_mesh_all.groupby('Semantic Types Label').describe()
chat_mesh_all_grup_sem = chat_mesh_all_grup.index.to_list()

new_name = 'chat_mesh_'+topic_label
exec(f'{new_name} = chat_mesh_all')
print(new_name, 'created')
print('   gpt terms:', gpt_all_ser.nunique())
print('   mesh found', chat_mesh_all['Preferred Label'].nunique())
print('   semantic groups',chat_mesh_all['Semantic Types Label'].nunique())

# GET ref mesh list
if os.path.exists('ref-mesh-archive/ref_mesh_'+topic_label+'.csv'):
    print('ref_mesh_'+topic_label+'.csv alreday exists')
    if overwrite_ref == True:
        globals()[new_name][['Preferred Label', 'Semantic Types Label']].to_csv('ref-mesh-archive/ref_mesh_'+topic_label+'.csv')
        print('overwritten')
else:
    globals()[new_name][['Preferred Label', 'Semantic Types Label']].to_csv('ref-mesh-archive/ref_mesh_'+topic_label+'.csv')
    print('ref_mesh_'+topic_label+'.csv saved')

chat_mesh_all[['Preferred Label','mesh_id','Semantic Types Label']].drop_duplicates()


# #compare
# print(len(pd.read_csv('ref-mesh-archive/ref_mesh_cvd.csv')['Preferred Label'].drop_duplicates()))
# print(len(pd.read_csv('ref-mesh-archive/ref_mesh_cvd_lipo.csv')['Preferred Label'].drop_duplicates()))
# gpt_all_ser

# In[ ]:


#check nunique:
check = 'nutri (no db)'
check2 = 'nutri_redux'
print(pd.read_csv('ref-mesh-archive/ref_mesh_'+check+'.csv')['Preferred Label'].nunique())
print(pd.read_csv('ref-mesh-archive/ref_mesh_'+check+'.csv')['Preferred Label'].nunique())


# ### clean output

# 

# In[ ]:


#check mesh ref
label = topic_label
dff = pd.read_csv('ref-mesh-archive/ref_mesh_'+label+'.csv')
dff['Preferred Label'].drop_duplicates()
check = 'Apop'
dff_check = dff[dff['Preferred Label'].str.contains(check)]
dff['Preferred Label'].drop_duplicates()
dff_check['Preferred Label'].drop_duplicates()


# In[ ]:


save_clean = True

#clean mesh ref
get_rid_list = ['oma','Apop','Cell','gen Ty','Cilia','Cytokinesis','DNA','Diabet','Fluorescein','Intercellular Signaling','Leukem','Lysergic Acid','Loeys-Dietz','Mitochondrial','Peptide Hormones','Toll-Like','Tumor']
dff_rid = dff.copy()
for get_rid in get_rid_list:
    mask = dff['Preferred Label'].str.contains(get_rid)
    dff_rid = dff_rid[-mask]

if save_clean == True:
    dff_rid.to_csv('ref-mesh-archive/ref_mesh_'+label+'.csv')

dff_rid['Preferred Label'].drop_duplicates()


# get rid log
# ob_bmi mesh ['Asthma','Cardio','Child H','Chromosome','Delivery','Diethyl','Drug','Employer','Electronic','Epigenomics','Ocular','Pulmon','Malignant','Organizat','Person','Occupat','Literacy','Level Se','h Prop','h Reso','es Res','for the Age','Indigen','h Syst','Deter','Antibod','Insulin-Like Growth Factor','Express','Patient','Pregnancy Complications','Rural','Sleep','Snow','Veter','Volunt','Water']
# 
# cvd_list = ['oma','Apop','Cell','gen Ty','Cilia','Cytokinesis','DNA','Diabet','Fluorescein','Intercellular Signaling','Leukem','Lysergic Acid','Loeys-Dietz','Mitochondrial','Peptide Hormones','Toll-Like','Tumor']

# ### ref_mesh semantic type ripper

# In[ ]:


#Sty groups
chat_mesh_all_grup


# In[ ]:


pd.DataFrame(gpt_nutri)


# In[ ]:


## Rip form complete ref mesh
save = False
rip_list = ['Cell',
            'Receptor',
            'Hormone',
            'Tissue',
            'Body Part, Organ, or Organ Component',
            'Congenital Abnormality',
            'Anatomical Abnormality',
            'Indicator, Reagent, or Diagnostic Aid']
rip_list += ['Neoplastic Process']
rip_list

#remove_value = chat_mesh_all[chat_mesh_all['Semantic Types Label'].isin(rip_list)]['Semantic Types Label'].reset_index(drop=True).to_list()
#count analyzer
#remove_value = ref_sty_less_count_sort[ref_sty_less_count_sort['mesh-count']> 56]['Semantic Types Label'].reset_index(drop=True).to_list()

mask = chat_mesh_all['Semantic Types Label'].isin(rip_list)
chat_mesh_all_redux = chat_mesh_all[-mask]

if save == True:
    chat_mesh_all_redux[['Preferred Label', 'Semantic Types Label']].to_csv('ref-mesh-archive/ref_mesh_'+topic_label+'_redux.csv')
chat_mesh_all_redux['Preferred Label'].drop_duplicates()


# #STY Ripper LOG:
# rip_list = ['Cell',
#             'Receptor',
#             'Hormone',
#             'Tissue',
#             'Body Part, Organ, or Organ Component',
#             'Congenital Abnormality',
#             'Anatomical Abnormality',
#             'Indicator, Reagent, or Diagnostic Aid']
# rip_list += ['Neoplastic Process']

# ### build ref mesh sterp by step
# #### neuro trial

# #ChatGPT list:
# #give me a comprehensive pyhton list (list= ['A','B','C',...]) of 100 real Mesh terms related to nervous system physiology
# #Now continue the  list with another 100 different Mesh tems in python list format.
# gpt_neuro_phys = [
#     'Action Potentials',
#     'Adrenoleukodystrophy',
#     'Alzheimer Disease',
#     'Amyloid beta-Peptides',
#     'Amyloid beta-Protein Precursor',
#     'Amyloidosis',
#     'Amyotrophic Lateral Sclerosis',
#     'Ataxia',
#     'Axons',
#     'Basal Ganglia',
#     'Blood-Brain Barrier',
#     'Brain',
#     'Brain Injuries',
#     'Brain-Derived Neurotrophic Factor',
#     'Calcium Signaling',
#     'Cerebellum',
#     'Cerebral Amyloid Angiopathy',
#     'Cerebral Cortex',
#     'Charcot-Marie-Tooth Disease',
#     'Cholinergic Neurons',
#     'Cognition',
#     'Copper',
#     'Cortical Spreading Depression',
#     'Creutzfeldt-Jakob Syndrome',
#     'Dementia',
#     'Dementia with Lewy Bodies',
#     'Dentate Gyrus',
#     'Depression',
#     'Dopamine',
#     'Dopaminergic Neurons',
#     'Electrophysiology',
#     'Encephalitis',
#     'Epilepsy',
#     'Excitatory Amino Acid Antagonists',
#     'Excitatory Postsynaptic Potentials',
#     'GABAergic Neurons',
#     'Gangliosidoses',
#     'Glial Fibrillary Acidic Protein',
#     'Gliosis',
#     'Glycine Plasma Membrane Transport Proteins',
#     'Hippocampus',
#     'Huntington Disease',
#     'Intracellular Signaling Peptides and Proteins',
#     'Ion Channels',
#     'Ischemia',
#     'Limbic System',
#     'Lipid Metabolism',
#     'Long-Term Potentiation',
#     'Manganese Poisoning',
#     'Membrane Potentials',
#     'Membrane Proteins',
#     'Memory',
#     'Mental Disorders',
#     'Microglia',
#     'Mitochondria',
#     'Molecular Chaperones',
#     'Molecular Motor Proteins',
#     'Motor Neurons',
#     'Multiple Sclerosis',
#     'Muscular Atrophy, Spinal',
#     'Myelin Sheath',
#     'Nerve Degeneration',
#     'Nerve Growth Factors',
#     'Nerve Regeneration',
#     'Nerve Tissue Proteins',
#     'Nervous System Diseases',
#     'Neural Inhibition',
#     'Neural Pathways',
#     'Neurofibrillary Tangles',
#     'Neurogenesis',
#     'Neuroglia',
#     'Neuronal Plasticity',
#     'Neurons',
#     'Neuroprotective Agents',
#     'Neurosecretion',
#     'Neurotransmitter Agents',
#     'Nicotinic Receptors',
#     'Nitric Oxide',
#     'Nuclear Proteins',
#     'Oligodendroglia',
#     'Oxidative Stress',
#     'Parkinson Disease',
#     'Phosphatidylinositol 3-Kinases',
#     'Phosphorylation',
#     'Pituitary-Adrenal System',
#     'Presenilin-1',
#     'Presynaptic Terminals',
#     'Protein Aggregates',
#     'Protein Folding',
#     'Protein Transport',
#     'Protein-Serine-Threonine Kinases',
#     'Protein-Tyrosine Kinases',
#     'Receptor, Adenosine A2A',
#     'Receptor, Muscarinic M1',
#     'Receptor, Nerve Growth Factor',
#     'Receptors, AMPA',
#     'Receptors, Dopamine',
#     'Receptors, Glutamate',
#     'Receptors, GABA-A',
#     'Receptors, GABA-B',
#     'Receptors, Metabotropic Glutamate',
#     'Receptors, N-Methyl-D-Aspartate',
#     'Receptors, Nicotinic',
#     'Receptors, Serotonin',
#     'Hippocampus',
#     'Amygdala',
#     'Cerebellum',
#     'Brain Stem',
#     'Substantia Nigra',
#     'Basal Ganglia',
#     'Corpus Striatum',
#     'Thalamus',
#     'Hypothalamus',
#     'Limbic System',
#     'Synapse',
#     'Neurotransmitter',
#     'Excitatory Postsynaptic Potential',
#     'Inhibitory Postsynaptic Potential',
#     'Action Potential',
#     'Ion Channels',
#     'Glutamate',
#     'GABA',
#     'Dopamine',
#     'Serotonin',
#     'Acetylcholine',
#     'Norepinephrine',
#     'Epinephrine',
#     'Histamine',
#     'Neuropeptides',
#     'Neuropeptide Y',
#     'Substance P',
#     'Cholecystokinin',
#     'Enkephalins',
#     'Beta-endorphins',
#     'Neuronal Plasticity',
#     'Long-Term Potentiation',
#     'Long-Term Depression',
#     'Neurogenesis',
#     'Neuronal Migration',
#     'Neuronal Differentiation',
#     'Axon Guidance',
#     'Neuronal Survival',
#     'Neural Circuits',
#     'Neural Networks',
#     'Neural Coding',
#     'Central Pattern Generators',
#     'Neuronal Oscillations',
#     'Brain Waves',
#     'Gamma Waves',
#     'Theta Waves',
#     'Delta Waves',
#     'Alpha Waves',
#     'Beta Waves',
#     'Sleep',
#     'REM Sleep',
#     'Non-REM Sleep',
#     'Sleep Disorders',
#     'Insomnia',
#     'Sleep Apnea',
#     'Narcolepsy',
#     'Cataplexy',
#     'Cataplexy',
#     'Restless Legs Syndrome',
#     'Parkinsonism',
#     'Parkinsonian Tremor',
#     'Parkinsonian Rigidity',
#     'Parkinsonian Bradykinesia',
#     'Parkinsonian Akinesia',
#     'Parkinsonian Gait',
#     'Parkinsonian Posture',
#     'Parkinsonism, Secondary',
#     'Huntington Disease',
#     'Choreiform Disorders',
#     'Chorea',
#     'Athetosis',
#     'Dystonia',
#     'Myoclonus',
#     'Tics',
#     'Tourette Syndrome',
#     'Friedreich Ataxia',
#     'Spinocerebellar Ataxias',
#     'Cerebellar Ataxia',
#     'Ataxia Telangiectasia',
#     'Restless Legs Syndrome',
#     'Periodic Limb Movements',
#     'Peripheral Neuropathies',
#     'Autonomic Neuropathy',
#     'Diabetic Neuropathies',
#     'Charcot-Marie-Tooth Disease',
#     'Guillain-Barre Syndrome',
#     'Chronic Inflammatory Demyelinating Polyneuropathy',
#     'Multiple Mononeuropathy',
#     'Mononeuritis Multiplex',
#     'Carpal Tunnel Syndrome',
#     'Brachial Plexus Neuropathies',
#     'Amyotrophic Lateral Sclerosis',
#     'Spinal Muscular Atrophy',
#     'Progressive Bulbar Palsy',
#     'Progressive Supranuclear Palsy',
#     'Corticobasal Degeneration',
#     'Multiple System Atrophy',
#     'Lewy Body Disease',
#     'Frontotemporal Dementia',
#     'Alzheimer Disease',
#     'Vascular Dementia']
# 
# #give me a comprehensive pyhton list (list= ['A','B','C',...]) of 100 real Mesh terms related to nervous neudegenarative diseases and neurologic pathology.
# #Now continue the  list with another 100 different Mesh tems in python list format .
# 
# gpt_neuro_path = ['Alzheimer Disease',
#                   'Amyotrophic Lateral Sclerosis',
#                   'Ataxia',
#                   'Basal Ganglia Diseases',
#                   'Brain Diseases',
#                   'Brain Infarction',
#                   'Brain Injury, Chronic',
#                   'Cauda Equina Syndrome',
#                   'Central Nervous System Diseases',
#                   'Cerebellar Diseases',
#                   'Cerebral Amyloid Angiopathy',
#                   'Cerebral Hemorrhage',
#                   'Cerebral Palsy',
#                   'Cerebrovascular Disorders',
#                   'Charcot-Marie-Tooth Disease',
#                   'Chorea',
#                   'Cortical Blindness',
#                   'Creutzfeldt-Jakob Syndrome',
#                   'Dementia',
#                   'Diffuse Axonal Injury',
#                   'Dystonia',
#                   'Encephalitis',
#                   'Encephalomyelitis',
#                   'Epilepsy',
#                   'Essential Tremor',
#                   'Friedreich Ataxia',
#                   'Gait Disorders, Neurologic',
#                   'Gerstmann-Straussler-Scheinker Disease',
#                   'Hallervorden-Spatz Syndrome',
#                   'Hemiplegia',
#                   'Hereditary Sensory and Autonomic Neuropathies',
#                   'Huntington Disease',
#                   'Hydrocephalus',
#                   'Intracranial Aneurysm',
#                   'Intracranial Arteriovenous Malformations',
#                   'Intracranial Embolism and Thrombosis',
#                   'Intracranial Hemorrhages',
#                   'Intracranial Hypertension',
#                   'Kuru',
#                   'Leigh Disease',
#                   'Leukoencephalopathies',
#                   'Lewy Body Disease',
#                   'Locked-In Syndrome',
#                   'Meningitis',
#                   'Meningoencephalitis',
#                   'Migraine Disorders',
#                   'Mild Cognitive Impairment',
#                   'Motor Neuron Disease',
#                   'Multiple Sclerosis',
#                   'Muscle Spasticity',
#                   'Muscular Atrophy',
#                   'Muscular Dystrophies',
#                   'Myasthenia Gravis',
#                   'Myoclonus',
#                   'Neuroaxonal Dystrophies',
#                   'Neurodegenerative Diseases',
#                   'Neurogenic Inflammation',
#                   'Neuromuscular Diseases',
#                   'Neuronal Ceroid-Lipofuscinoses',
#                   'Neuropathic Pain',
#                   'Neuropil Disorders',
#                   'Neurosyphilis',
#                   'Normal Pressure Hydrocephalus',
#                   'Parkinson Disease',
#                   'Peripheral Nervous System Diseases',
#                   'Phantom Limb',
#                   'Poliomyelitis',
#                   'Post-Concussion Syndrome',
#                   'Prion Diseases',
#                   'Progressive Bulbar Palsy',
#                   'Progressive Supranuclear Palsy',
#                   'Reflex Sympathetic Dystrophy',
#                   'Restless Legs Syndrome',
#                   'Retinal Diseases',
#                   'Retinal Ganglion Cell Loss',
#                   'Schilder Disease',
#                   'Shy-Drager Syndrome',
#                   'Spinal Cord Diseases',
#                   'Spinal Cord Injuries',
#                   'Spinal Muscular Atrophies of Childhood',
#                   'Stiff-Person Syndrome',
#                   'Stroke',
#                   'Subarachnoid Hemorrhage',
#                   'Subdural Hematoma',
#                   'Syringomyelia',
#                   'Tay-Sachs Disease',
#                   'Tourette Syndrome' ]
# 
# gpt_neuro = gpt_neuro_phys + gpt_neuro_path
# gpt_neuro_ser = pd.Series(gpt_neuro)
# 
# gpt_neuro_ser = gpt_neuro_ser.drop_duplicates()
# print(gpt_neuro_ser.nunique())
# gpt_neuro_ser
# 
# print('physiologic topic')
# #mother_df = mesh_large_df_sty
# mother_df = mesh_litvar_sty
# #get mesh
# chat_mesh_neuro_phys = pd.DataFrame(columns=mother_df.columns) #Empty DataFrame, same structure
# time1 = datetime.now()
# for i in gpt_neuro_phys:
#     meshy = mother_df[mother_df['Preferred Label'].str.contains(i).fillna(False)]
#     #chat_mesh_neuro = chat_mesh_neuro.append([meshy])
#     chat_mesh_neuro_phys = pd.concat([meshy,chat_mesh_neuro_phys])
#     #chatmesh_var = pd.concat([meshy, chatmesh_var.loc[:]])
# time2 = datetime.now()
# print(time2-time1)
# 

# In[ ]:


print('physiologic topic')
#mother_df = mesh_large_df_sty
mother_df = mesh_litvar_sty
#get mesh
chat_mesh_neuro_phys = pd.DataFrame(columns=mother_df.columns) #Empty DataFrame, same structure
time1 = datetime.now()
for i in gpt_neuro_phys:
    meshy = mother_df[mother_df['Preferred Label'].str.contains(i).fillna(False)]
    #chat_mesh_neuro = chat_mesh_neuro.append([meshy])
    chat_mesh_neuro_phys = pd.concat([meshy,chat_mesh_neuro_phys])
    #chatmesh_var = pd.concat([meshy, chatmesh_var.loc[:]])
time2 = datetime.now()
print(time2-time1)


# In[ ]:


print('mesh number:',chat_mesh_neuro_phys['Preferred Label'].nunique())
chat_mesh_neuro_phys


# neuro_phys on global, mesh: 824
# neuro_phys on litvar, mesh: 438

# In[ ]:


chat_mesh_neuro_phys_grup = chat_mesh_neuro_phys.groupby('Semantic Types Label').describe()
chat_mesh_neuro_phys_grup.columns = chat_mesh_neuro_phys_grup.columns.to_flat_index()
chat_mesh_neuro_phys_sem = chat_mesh_neuro_phys_grup.index.to_list()
chat_mesh_neuro_phys_grup


# In[ ]:


print('pahologic topic')
#get mesh
chat_mesh_neuro_path = pd.DataFrame(columns=mesh_large_df_sty.columns) #Empty DataFrame, same structure
time1 =datetime.now()
for i in gpt_neuro_path:
    meshy = mesh_large_df_sty[mesh_large_df_sty['Preferred Label'].str.contains(i).fillna(False)]
    #chat_mesh_neuro = chat_mesh_neuro.append([meshy])
    chat_mesh_neuro_path = pd.concat([meshy,chat_mesh_neuro_path])
    #chatmesh_var = pd.concat([meshy, chatmesh_var.loc[:]])
time2 =datetime.now()
print(time2-time1)


# In[ ]:


print(chat_mesh_neuro_path['Preferred Label'].nunique())
chat_mesh_neuro_path


# In[ ]:


#analyze semantics
chat_mesh_neuro_path_grup = chat_mesh_neuro_path.groupby('Semantic Types Label').describe()
chat_mesh_neuro_path_sem = chat_mesh_neuro_path_grup.index.to_list()
chat_mesh_neuro_path_grup


# In[ ]:


#get sum of semantic types
chat_mesh_sem_sum = chat_mesh_neuro_path_sem + chat_mesh_neuro_phys_sem
chat_mesh_sem_sum_ser = pd.Series(chat_mesh_sem_sum)
print(chat_mesh_sem_sum_ser.nunique())
chat_mesh_sem_sum_ser


# In[ ]:


#get mesh
chat_mesh_neuro = pd.DataFrame(columns=mesh_large_df_sty.columns) #Empty DataFrame, same structure
for i in gpt_neuro_ser:
    meshy = mesh_large_df_sty[mesh_large_df_sty['Preferred Label'].str.contains(i).fillna(False)]
    #chat_mesh_neuro = chat_mesh_neuro.append([meshy])
    chat_mesh_neuro = pd.concat([meshy,chat_mesh_neuro])
    #chatmesh_var = pd.concat([meshy, chatmesh_var.loc[:]])


# In[ ]:


print('mesh', chat_mesh_neuro['Preferred Label'].nunique())
chat_mesh_neuro_sem = chat_mesh_neuro['Semantic Types Label'].drop_duplicates()
print('sty',chat_mesh_neuro['Semantic Types Label'].nunique())
chat_mesh_neuro


# In[ ]:


#analyze semantics
chat_mesh_neuro.groupby('Semantic Types Label').describe()


# In[ ]:


chat_mesh_neuro['Preferred Label'].nunique()#.columns


# In[ ]:


# GET ref mesh list
if os.path.exists('ref-mesh-archive/ref_mesh_neuro.csv'):
    print('pass')
    pass
else:
    chat_mesh_neuro[['Preferred Label', 'Semantic Types Label']].to_csv('ref-mesh-archive/ref_mesh_neuro.csv')


# In[ ]:


from ast import literal_eval
# grupuby describe:
chat_mesh_neuro_count = chat_mesh_neuro.groupby('Semantic Types Label').describe()
old_col = str(chat_mesh_neuro_count.columns.to_list())
new_col = old_col.replace("\', \'", "_")
new_col = new_col.replace("(", "")
new_col = new_col.replace(")", "")
chat_mesh_neuro_count.columns = literal_eval(new_col)
chat_mesh_neuro_count.sort_values('Preferred Label_count', inplace=True, ascending=False)
chat_mesh_neuro_count


# In[ ]:


## Rip form complete ref mesh
rip_list = ['Cell', 'Receptor', 'Hormone', 'Tissue', 'Body Part, Organ, or Organ Component', 'Congenital Abnormality', 'Anatomical Abnormality', 'Indicator, Reagent, or Diagnostic Aid']

#new_ref_sty = ref_sty[ref_sty['Semantic Types Label'].isin(remove_value)]
mask = chat_mesh_neuro['Semantic Types Label'].isin(rip_list)
new_ref_sty_neuro = chat_mesh_neuro[-mask]
new_ref_sty_neuro_sty = new_ref_sty_neuro['Semantic Types Label'].drop_duplicates()

if os.path.exists('ref-mesh-archive/new_ref_sty_neuro.csv'):
    pass
else:
    new_ref_sty_neuro.to_csv('ref-mesh-archive/new_ref_sty_neuro.csv')


# In[ ]:


print(new_ref_sty_neuro['Preferred Label'].nunique())
chat_mesh_neuro_count = new_ref_sty_neuro.groupby('Semantic Types Label').describe()
old_col = str(chat_mesh_neuro_count.columns.to_list())
new_col = old_col.replace("\', \'", "_")
new_col = new_col.replace("(", "")
new_col = new_col.replace(")", "")
chat_mesh_neuro_count.columns = literal_eval(new_col)
chat_mesh_neuro_count.sort_values('Preferred Label_count', inplace=True, ascending=False)
chat_mesh_neuro_count


# In[ ]:


new_ref_sty_neuro_sty


# 

# #### nutri trial

# In[ ]:


#ChatGPT list
gpt_nutri = ["Nutrition",
"Diet",
"Energy Metabolism",
"Macronutrients",
"Carbohydrates",
"Lipids",
"Micronutrients",
"Vitamins",
"Minerals",
"Trace Elements",
"Antioxidants",
"Nutrition Disorders",
"Malnutrition",
"Obesity",
"Overweight",
"Energy Balance",
"Body Composition",
"Body Mass Index",
"Adiposity",
"Nutrient Absorption",
"Nutrient Transport",
"Nutrient Metabolism",
"Nutrient Requirements",
"Nutritional Status",
"Nutritional Assessment",
"Nutritional Surveillance",
"Nutritional Epidemiology",
"Nutritional Physiology",
"Nutrient Interactions",
"Nutrient Toxicity",
"Nutrient Deficiency",
"Food Intake",
"Appetite",
"Food Preferences",
"Feeding Behavior",
"Food Habits",
"Food Choice",
"Eating Disorders",
"Nutritional Genomics",
"Nutrigenetics",
"Nutrigenomics",
"Epigenetics",
"Genetic Variation",
"Genotype",
"Phenotype",
"Gene Expression",
"Single Nucleotide Polymorphisms",
"MicroRNAs",
"Enzymes",
"Hormones",
"Insulin",
"Glucagon",
"Leptin",
"Adiponectin",
"Ghrelin",
"Cortisol",
"Thyroid Hormones",
"Growth Hormones",
"Inflammation",
"Oxidative Stress",
"Immune System",
"Aging",
"Development",
"Reproduction",
"Maternal Nutrition",
"Lactation",
"Infant Nutrition",
"Childhood Nutrition",
"Adolescent Nutrition",
"Elderly Nutrition",
"Sports Nutrition",
"Physical Activity",
"Exercise",
"Endurance",
"Strength",
"Bodybuilding",
"Resistance Training",
"Nutrition Supplementation",
"Nutrient Timing",
"Diet Supplements",
"Ergogenic Aids",
"Weight Management",
"Low Carbohydrate Diets",
"High Protein Diets",
"Low Fat Diets",
"Vegetarian Diets",
"Vegan Diets",
"Gluten-Free Diets",
"Dairy-Free Diets",
"Soy-Free Diets",
"Detox Diets",
"Intermittent Fasting",
"Time-Restricted Feeding",
"Meal Replacements",
"Energy Bars",
"Meal Supplements",
"Hydration",
"Water Requirements",
"Nutrient Density",
"Food Composition",
"Food Processing",
"Food Preservation",
"Food Safety",
"Food Additives",
"Food Contaminants",
"Food Fortification",
"Food Supplementation",
"Food Sources",
"Food Availability",
"Food Access",
"Food Insecurity",
"Food Assistance Programs",
"Food Labeling",
"Nutrient Recommendations",
"Dietary Reference Intakes",
"Recommended Dietary Allowances",
"Adequate Intake",
"Tolerable Upper Intake Level",
"Nutrient Recommendations for Special Populations",
"Pregnancy",
"Infancy",
"Childhood",
"Adolescence",
"Elderly",
"Athletes",
"Chronic Diseases",
"Cardiovascular Disease",
"Hypertension",
"Diabetes",
"Metabolic Syndrome",
"Cancer",
"Gastrointestinal Diseases",
"Inflammatory Bowel Disease",
"Non-alcoholic Fatty Liver Disease",
"Nutrition Support",
"Enteral Nutrition",
"Parenteral Nutrition",
"Artificial Feeding",
"Nutritional Requirements in Critical Illness",
"Nutritional Support in Surgery",
"Nutritional Management of Chronic Wounds",
"Nutritional Interventions",
"Nutrition Education",
"Food Skills Training",
"Cooking Classes",
"Meal Planning",
"Food Budgeting",
"Food Shopping",
"Food Store Environment",
"Food Marketing",
"Food Advertising",
"Food Industry",
"Food Systems",
"Agricultural Systems",
"Livestock Production",
"Seafood Production",
"Food Trade",
"Food Policy",
"Food Security",
"Food Sovereignty",
"Food Justice",
"Food System Sustainability",
"Food Waste",
"Sustainable Diets",
"Plant-based Diets",
"Mediterranean Diet",
"Paleolithic Diet",
"Traditional Diets",
"Fermented Foods",
"Probiotics",
"Prebiotics",
"Gut Microbiota",
"Gut Health",
"Gut-Brain Axis",
"Neural Control of Feeding",
"Brain-Derived Neurotrophic Factor",
"Cholecystokinin",
"Peptide YY",
"Amylin",
"Glucagon-Like Peptide 1",
"Oxyntomodulin",
"PYY3-36",
"Agouti-Related Peptide",
"Neuropeptide Y",
"Cocaine- and Amphetamine-Regulated Transcript",
"Melanin-Concentrating Hormone",
"Orexins",
"Pro-opiomelanocortin",
"Corticotropin-Releasing Hormone",
"Stress",
"Sleep",
"Circadian Rhythms",
"Light-Dark Cycles",
"Social Environment",
"Culture",
"Food Traditions",
"Food Beliefs and Practices",
"Carbohydrate Metabolism",
"Lipid Metabolism",
"Protein Metabolism",
"Nitrogen Balance",
"Macronutrient Metabolism",
"Micronutrient Metabolism",
"Vitamin Metabolism",
"Mineral Metabolism",
"Trace Element Metabolism",
"Antioxidant Status",
"Gut-Liver Axis",
"Adipose Tissue",
"Body Weight",
"Waist Circumference",
"Body Fat Distribution",
"Resistin",
"Insulin Resistance",
"Insulin Sensitivity",
"Insulin Secretion",
"Glucose Homeostasis",
"Hyperglycemia",
"Hypoglycemia",
"Type 2 Diabetes",
"Impaired Glucose Tolerance",
"Impaired Fasting Glucose",
"Dyslipidemia",
"Hypercholesterolemia",
"Hypertriglyceridemia",
"Low-Density Lipoprotein",
"High-Density Lipoprotein",
"Lipoprotein(a)",
"Atherosclerosis",
"Kidney Disease",
"Renal Failure",
"Renal Insufficiency",
"Chronic Kidney Disease",
"Endocrine Disorders",
"Parathyroid Hormones",
"Adrenal Hormones",
"Pancreatic Hormones",
"Reproductive Hormones",
"Menopause",
"Andropause",
"Bone Health",
"Osteoporosis",
"Fractures",
"Bone Mineral Density",
"Bone Turnover",
"Calcium Homeostasis",
"Vitamin D Metabolism",
"Muscle Mass",
"Sarcopenia",
"Muscle Strength",
"Physical Performance",
"Physical Function",
"Physical Disability",
"Endurance Training",
"High-Intensity Interval Training",
"Physical Inactivity",
"Sedentary Behavior",
"Geriatrics",
"Age-Related Changes",
"Cognitive Function",
"Dementia",
"Alzheimer's Disease",
"Parkinson's Disease",
"Huntington's Disease",
"Multiple Sclerosis",
"Neurodegeneration",
"Nerve Function",
"Neural Plasticity",
"Sleep Disorders",
"Mood",
"Anxiety",
"Depression",
"Psychological Well-Being",
"Quality of Life",
"Life Expectancy",
"Health Span",
"Hunger",
"Satiety",
"Body Image",
"Anorexia Nervosa",
"Bulimia Nervosa",
"Binge Eating Disorder",
"Dietary Patterns",
"Dietary Guidelines",
"Nutrient Adequacy",
"Nutrient Bioavailability",
"Nutrient Storage",
"Nutrient Utilization",
"Food Preference",
"Food Aversion",
"Food Sensitivities",
"Food Allergies",
"Food Intolerance",
"Food-Drug Interactions",
"Nutrient Supplements",
"Nutrient Fortification",
"Nutrient Deficiencies",
"Protein-Energy Malnutrition",
"Kwashiorkor",
"Marasmus",
"Vitamin Deficiencies",
"Vitamin A Deficiency",
"Vitamin C Deficiency",
"Vitamin D Deficiency",
"Vitamin E Deficiency",
"Vitamin K Deficiency",
"Thiamin Deficiency",
"Riboflavin Deficiency",
"Niacin Deficiency",
"Vitamin B6 Deficiency",
"Vitamin B12 Deficiency",
"Folate Deficiency",
"Mineral Deficiencies",
"Iron Deficiency",
"Zinc Deficiency",
"Calcium Deficiency",
"Magnesium Deficiency",
"Phosphorus Deficiency",
"Potassium Deficiency",
"Sodium Deficiency",
"Iodine Deficiency",
"Selenium Deficiency",
"Chromium Deficiency",
"Copper Deficiency",
"Manganese Deficiency",
"Molybdenum Deficiency",
"Diet-Gene Interactions",
"Genotype-Phenotype Interactions",
"Genome-Wide Association Studies",
"Personalized Nutrition",
"Precision Nutrition",
"Nutritional Status Assessment",
"Dietary Assessment",
"Biomarkers of Nutritional Status",
"Serum Concentrations",
"Plasma Concentrations",
"Red Blood Cell Concentrations",
"White Blood Cell Concentrations",
"Hair Concentrations",
"Nail Concentrations",
"Saliva Concentrations",
"Urine Concentrations",
"Stool Concentrations",
"Gastrointestinal Function",
"Gut Microbiome",
"Gut Permeability",
"Irritable Bowel Syndrome",
"Short-Chain Fatty Acids",
"Synbiotics",
"Fiber",
"Resistant Starch",
"Fermentable Carbohydrates",
"Fermentation",
"Digestion",
"Fat Metabolism",
"Nucleic Acid Metabolism",
"Cholesterol Metabolism",
"Triglycerides",
"Fatty Acids",
"Chylomicrons",
"Lipoproteins",
"Low-Density Lipoprotein Cholesterol",
"High-Density Lipoprotein Cholesterol",
"Very Low-Density Lipoprotein Cholesterol",
"Intermediate-Density Lipoprotein Cholesterol",
"Lipid Peroxidation",
"Lipid Oxidation",
"Free Radicals",
"Reactive Oxygen Species",
"Mitochondrial Function",
"Glucose",
"Glycogen",
"Gluconeogenesis",
"Type 2 Diabetes Mellitus",
"Endocrine Pancreas",
"Islets of Langerhans",
"Beta Cells",
"Alpha Cells",
"Gastrointestinal Hormones",
"Secretin",
"Gastrin",
"Somatostatin",
"Glucagon-Like Peptide-1",
"Interleukins",
"Cytokines",
"Tumor Necrosis Factor-Alpha",
"Adipokines",
"Hormone-Sensitive Lipase",
"Peroxisome Proliferator-Activated Receptor-Gamma",
"Peroxisome Proliferator-Activated Receptor-Alpha",
"Nuclear Receptor Subfamily 1, Group H, Member 3",
"Steroid Hormones",
"Estrogens",
"Androgens",
"Progestogens",
"Adrenal Cortex Hormones",
"Adrenal Medulla Hormones",
"Calcitonin",
"Vitamin D",
"Phosphorus",
"Calcium",
"Magnesium",
"Bone Metabolism",
"Bone Density",
"Bone Mass",
"Innate Immunity",
"Adaptive Immunity",
"Antibodies",
"Complement System",
"Phagocytosis",
"Inflammatory Response",
"Allergic Reactions",
"Hypersensitivity Reactions",
"Anaphylaxis",
"Immune Tolerance",
"Immune Suppression",
"Autoimmunity",
"Immune-Mediated Diseases",
"Transplantation",
"Graft Rejection",
"T-Lymphocytes",
"B-Lymphocytes",
"Natural Killer Cells",
"Monocytes",
"Macrophages",
"Dendritic Cells",
"Lymphocytes",
"Interferons",
"Cytotoxic T-Lymphocytes",
"Regulatory T-Lymphocytes"]
gpt_nutri_ser = pd.Series(gpt_nutri).drop_duplicates()

#get mesh
chat_mesh_nutri = pd.DataFrame(columns=mesh_large_df_sty.columns) #Empty DataFrame, same structure
for i in gpt_nutri_ser:
    meshy = mesh_large_df_sty[mesh_large_df_sty['Preferred Label'].str.contains(i).fillna(False)]
    #chat_mesh_nutri = chat_mesh_nutri.append([meshy])
    chat_mesh_nutri = pd.concat([meshy, chat_mesh_nutri])
    #chatmesh_var = pd.concat([meshy, chatmesh_var.loc[:]])


# Defined new_fer_mesh.csv for nutritional topic --> import it and model neuro un it

# new_nutri_refmesh = pd.read_csv('ref-mesh-archive/new_ref_mesh_large.csv', index_col=0)
# new_nutri_refmesh_count = new_nutri_refmesh.groupby('Semantic Types Label').describe()
# new_nutri_refmesh['Preferred Label'].drop_duplicates()

# In[ ]:


# grupuby describe:
new_nutri_refmesh = pd.read_csv('ref-mesh-archive/new_ref_mesh_large_corrected.csv', index_col=0)
new_nutri_refmesh_sty = new_nutri_refmesh['Semantic Types Label'].drop_duplicates()

new_nutri_refmesh_count = new_nutri_refmesh.groupby('Semantic Types Label').describe()
old_col = str(new_nutri_refmesh_count.columns.to_list())
new_col = old_col.replace("\', \'", "_")
new_col = new_col.replace("(", "")
new_col = new_col.replace(")", "")
new_nutri_refmesh_count.columns = literal_eval(new_col)
new_nutri_refmesh_count.sort_values('Preferred Label_count', inplace=True, ascending=False)
new_nutri_refmesh_count


# In[ ]:


new_nutri_refmesh['Preferred Label'].nunique()


# In[ ]:


#evalue STY overlapping
neuro_nutri_sty_overrlap  = new_ref_sty_neuro_sty[new_ref_sty_neuro_sty.isin(new_nutri_refmesh_sty)]
print( new_ref_sty_neuro_sty.nunique())
print( new_nutri_refmesh_sty.nunique())
print(neuro_nutri_sty_overrlap.nunique())


# 

# # Generate Random Mesh

# ## Random grpm mesh list from MESH-STY

#     list1 = [1, 2, 3, 4]
#     list2 = [3, 4, 5, 6]
#     list3 = [4, 5, 6, 7]
# 
#     lists = [list1, list2, list3]
#     num_lists = len(lists)
# 
#     # Initialize a 2D list of zeros with dimensions equal to the number of lists
#     cooccur_matrix = [[0] * num_lists for i in range(num_lists)]
# 
#     # Loop over all pairs of lists and count the number of co-occurring elements
#     for i in range(num_lists):
#         for j in range(num_lists):
#             if i == j:
#                 cooccur_matrix[i][j] = len(lists[i])
#             else:
#                 cooccur_matrix[i][j] = len(set(lists[i]) & set(lists[j]))
# 
#     # Print the resulting matrix with row and column headers
#     print(', '.join([''] + ['list{}'.format(i+1) for i in range(num_lists)]))
#     for i in range(num_lists):
#         row = ['list{}'.format(i+1)]
#         row.extend(cooccur_matrix[i])
#         print(', '.join(str(x) for x in row))
# 

# In[ ]:


mesh_large_df_sty


# In[ ]:


# Generate Random mesh list from complete df
print(mesh_large_df_sty['Preferred Label'].nunique())
mesh_large_df_sty_mesh = mesh_large_df_sty['Preferred Label'].drop_duplicates()
random_path = 'ref-mesh-archive/random_lists/'

size = 603
sample_01 = mesh_large_df_sty_mesh.sample(n=size, random_state = 78954)
sample_02 = mesh_large_df_sty_mesh.sample(n=size, random_state = 12245)
sample_03 = mesh_large_df_sty_mesh.sample(n=size, random_state = 87498)
sample_04 = mesh_large_df_sty_mesh.sample(n=size, random_state = 56798)
sample_05 = mesh_large_df_sty_mesh.sample(n=size, random_state = 34565)
sample_06 = mesh_large_df_sty_mesh.sample(n=size, random_state = 76523)
sample_07 = mesh_large_df_sty_mesh.sample(n=size, random_state = 78968)
sample_08 = mesh_large_df_sty_mesh.sample(n=size, random_state = 56845)
sample_09 = mesh_large_df_sty_mesh.sample(n=size, random_state = 76624)
sample_10 = mesh_large_df_sty_mesh.sample(n=size, random_state = 23845)

sample_01.to_csv(random_path+'random_01.csv')
sample_02.to_csv(random_path+'random_02.csv')
sample_03.to_csv(random_path+'random_03.csv')
sample_04.to_csv(random_path+'random_04.csv')
sample_05.to_csv(random_path+'random_05.csv')
sample_06.to_csv(random_path+'random_06.csv')
sample_07.to_csv(random_path+'random_07.csv')
sample_08.to_csv(random_path+'random_08.csv')
sample_09.to_csv(random_path+'random_09.csv')
sample_10.to_csv(random_path+'random_10.csv')


# ## Random grpm mesh list from MESH-STY-LITVAR
# Generate Random Mesh list from grpm df
# 
#     grpm_db_mesh = pd.read_csv('ref-mesh-archive/grpm_db_mesh.csv')
#     grpm_db_mesh

# In[ ]:


import pandas as pd
import os
if os.path.exists('ref-mesh-archive/grpm_db_mesh.csv'):
    grpm_db_mesh = pd.read_csv('ref-mesh-archive/grpm_db_mesh.csv',index_col=0)
    print('grpm mesh imported from archive')
else:
    grpmdb_table = pd.read_csv(r'H:\Il mio Drive\GRPM db\GRPM db - Vitam\grpm_db_pcg\grpm_table_output.csv')
    grpm_db_mesh = pd.DataFrame(grpmdb_table.mesh.drop_duplicates())
    type(grpm_db_mesh)


# rand = pd.read_csv('ref-mesh-archive/random_grpm_03.csv')
# mask= grpm_db_mesh.mesh.isin(rand['mesh'])
# grpm_db_mesh[mask]

# In[ ]:


# randomize 10 samples
size = 450
sample_grpm_01 = grpm_db_mesh.mesh.sample(n=size, random_state = 78954)
sample_grpm_02 = grpm_db_mesh.mesh.sample(n=size, random_state = 12245)
sample_grpm_03 = grpm_db_mesh.mesh.sample(n=size, random_state = 87498)
sample_grpm_04 = grpm_db_mesh.mesh.sample(n=size, random_state = 56798)
sample_grpm_05 = grpm_db_mesh.mesh.sample(n=size, random_state = 34565)
sample_grpm_06 = grpm_db_mesh.mesh.sample(n=size, random_state = 76523)
sample_grpm_07 = grpm_db_mesh.mesh.sample(n=size, random_state = 78968)
sample_grpm_08 = grpm_db_mesh.mesh.sample(n=size, random_state = 56845)
sample_grpm_09 = grpm_db_mesh.mesh.sample(n=size, random_state = 76624)
sample_grpm_00 = grpm_db_mesh.mesh.sample(n=size, random_state = 23845)


# In[ ]:


#check fidelity
mask= grpm_db_mesh.mesh.isin(sample_grpm_01)
grpm_db_mesh[mask]
#type(sample_grpm_10)


# In[ ]:


# save random samples
random_path = 'ref-mesh-archive/random_lists/'

if os.path.exists('ref-mesh-archive/random_grpm_00.csv'):
    sample_grpm_01 = pd.read_csv(random_path+'random_grpm_01.csv',index_col=0)
    sample_grpm_02 = pd.read_csv(random_path+'random_grpm_02.csv',index_col=0)
    sample_grpm_03 = pd.read_csv(random_path+'random_grpm_03.csv',index_col=0)
    sample_grpm_04 = pd.read_csv(random_path+'random_grpm_04.csv',index_col=0)
    sample_grpm_05 = pd.read_csv(random_path+'random_grpm_05.csv',index_col=0)
    sample_grpm_06 = pd.read_csv(random_path+'random_grpm_06.csv',index_col=0)
    sample_grpm_07 = pd.read_csv(random_path+'random_grpm_07.csv',index_col=0)
    sample_grpm_08 = pd.read_csv(random_path+'random_grpm_08.csv',index_col=0)
    sample_grpm_09 = pd.read_csv(random_path+'random_grpm_09.csv',index_col=0)
    sample_grpm_00 = pd.read_csv(random_path+'random_grpm_00.csv',index_col=0)

save_rand = False
if save_rand == True:
    sample_grpm_01.to_csv(random_path+'random_grpm_01.csv')
    sample_grpm_02.to_csv(random_path+'random_grpm_02.csv')
    sample_grpm_03.to_csv(random_path+'random_grpm_03.csv')
    sample_grpm_04.to_csv(random_path+'random_grpm_04.csv')
    sample_grpm_05.to_csv(random_path+'random_grpm_05.csv')
    sample_grpm_06.to_csv(random_path+'random_grpm_06.csv')
    sample_grpm_07.to_csv(random_path+'random_grpm_07.csv')
    sample_grpm_08.to_csv(random_path+'random_grpm_08.csv')
    sample_grpm_09.to_csv(random_path+'random_grpm_09.csv')
    sample_grpm_00.to_csv(random_path+'random_grpm_00.csv')


# In[ ]:


type(sample_grpm_01)


# In[ ]:


sample_grpm_01_list = sample_grpm_01.to_list()
sample_grpm_02_list = sample_grpm_02.to_list()
sample_grpm_03_list = sample_grpm_03.to_list()
sample_grpm_04_list = sample_grpm_04.to_list()
sample_grpm_05_list = sample_grpm_05.to_list()
sample_grpm_06_list = sample_grpm_06.to_list()
sample_grpm_07_list = sample_grpm_07.to_list()
sample_grpm_08_list = sample_grpm_08.to_list()
sample_grpm_09_list = sample_grpm_09.to_list()
sample_grpm_00_list = sample_grpm_00.to_list()

lists= [sample_grpm_01_list,
        sample_grpm_02_list,
        sample_grpm_03_list,
        sample_grpm_04_list,
        sample_grpm_05_list,
        sample_grpm_06_list,
        sample_grpm_07_list,
        sample_grpm_08_list,
        sample_grpm_09_list,
        sample_grpm_00_list]

num_lists = len(lists)


# ## Coocuurance method

# In[ ]:


#(1) Initialize a 2D list of zeros with dimensions equal to the number of lists
cooccur_matrix = [[0] * num_lists for i in range(num_lists)]
cooccur_matrix = np.zeros((len(lists), len(lists)), dtype=int)

type(cooccur_matrix)
print(len(set(lists[1]) & set(lists[1])))
print(len(set(lists[1]) & set(lists[3])))

# Loop over all pairs of lists and count the number of co-occurring elements
for i in range(num_lists):
    for j in range(num_lists):
        if i == j:
            cooccur_matrix[i][j] = len(set(lists[i]))
        else:
            cooccur_matrix[i][j] = len(set(lists[i]) & set(lists[j])) # use set() to extract unique elements

print(type(cooccur_matrix))

# Convert the 2D list to a Pandas DataFrame
cooccur_df = pd.DataFrame(cooccur_matrix,
                          columns=['random{}'.format(i+1) for i in range(num_lists)],
                          index=['random{}'.format(i+1) for i in range(num_lists)])

# Print the resulting DataFrame
cooccur_df


# ## Read Random grpm mesh list without rep

# In[ ]:


# COOCCURANCE MATRIX MODULE------------
# second matrix build check
sample_grpm_01_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep1.csv')
sample_grpm_02_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep2.csv')
sample_grpm_03_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep3.csv')
sample_grpm_04_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep4.csv')
sample_grpm_05_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep5.csv')
sample_grpm_06_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep6.csv')
sample_grpm_07_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep7.csv')
sample_grpm_08_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep8.csv')
sample_grpm_09_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep9.csv')
sample_grpm_00_norep = pd.read_csv(random_path+'random_grpm_mesh_norep/random_grpm_norep10.csv')

lists= [sample_grpm_01_norep.mesh.to_list(),
        sample_grpm_02_norep.mesh.to_list(),
        sample_grpm_03_norep.mesh.to_list(),
        sample_grpm_04_norep.mesh.to_list(),
        sample_grpm_05_norep.mesh.to_list(),
        sample_grpm_06_norep.mesh.to_list(),
        sample_grpm_07_norep.mesh.to_list(),
        sample_grpm_08_norep.mesh.to_list(),
        sample_grpm_09_norep.mesh.to_list(),
        sample_grpm_00_norep.mesh.to_list()]

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
cooccur_df = pd.DataFrame(cooccur_matrix, columns=['random{}'.format(i+1) for i in range(num_lists)], index=['random{}'.format(i+1) for i in range(num_lists)])

# Print the resulting DataFrame
cooccur_df


# In[ ]:


# Print the resulting matrix with row and column headers
print(', '.join([''] + ['list{}'.format(i+1) for i in range(num_lists)]))

for i in range(num_lists):
    row = ['list{}'.format(i+1)]
    row.extend(cooccur_matrix[i])
    print(', '.join(str(x) for x in row))


# crea 2 set da 5 random mesh senza ripetizioni

# ## Create Random grpm mesh list without rep from MESH-STY-LITVAR

# In[ ]:


# randomize
sample_grpm_full = grpm_db_mesh.mesh.sample(n=len(grpm_db_mesh), random_state = 782954)
sample_grpm_full


# In[ ]:


import numpy as np
# chunk size
size = 450
# Split the DataFrame into smaller DataFrames of chunk length 450
chunks = np.array_split(sample_grpm_full, len(sample_grpm_full) // 450 + 1)

# Print the number of chunks and the length of each chunk
for i, chunk in enumerate(chunks):
    print('Chunk {}: {} rows'.format(i+1, len(chunk)))
    pass


# In[ ]:


# Save each DataFrame chunk to a separate CSV file
for i, chunk in enumerate(chunks):
    chunk.to_csv('ref-mesh-archive/random_grpm_norep{}.csv'.format(i+1), index=False)


# ## Try random mesh list based on STY proportions

# In[ ]:


#define ref_sty
ref_sty = pd.concat([new_ref_sty_neuro_sty, new_nutri_refmesh_sty]).drop_duplicates()
strip  = ['Enzyme', 'Bacterium','Organism']
strip2  = ['Disease or Syndrome']
result = [item for item in new_nutri_refmesh_sty.to_list() if item not in strip]
result2 = [item for item in result if item not in strip2]
ref_sty_nutri = pd.Series(result)
ref_sty_nutri_rest = pd.Series(result2)

#mother df
print('mother df',mesh_large_df_sty['Preferred Label'].nunique(), len(mesh_large_df_sty['Preferred Label']))

#filter mother df with ref sty
mesh_large_df_sty_subset = mesh_large_df_sty[mesh_large_df_sty['Semantic Types Label'].isin(ref_sty_nutri)]
print('filtered df', mesh_large_df_sty_subset['Preferred Label'].nunique(), len(mesh_large_df_sty_subset['Preferred Label']))
#mesh_large_df_sty_subset['Semantic Types Label'].drop_duplicates()
mesh_large_df_sty_subset_rest = mesh_large_df_sty_subset[mesh_large_df_sty_subset['Semantic Types Label'].isin(ref_sty_nutri_rest)]

mesh_large_df_sty_subset_rest


# In[ ]:


mesh_large_df_sty_subset.groupby('Semantic Types Label').describe()


# In[ ]:


mesh_large_df_sty['Preferred Label'].nunique()


# In[ ]:


#get random mesh list from filtered df
sample_rest = 603 - 255
disease_sample = 255

mesh_large_df_sty_subset_disease = mesh_large_df_sty_subset.loc[mesh_large_df_sty_subset['Semantic Types Label'] == 'Disease or Syndrome']
random_mesh_disease_01 = mesh_large_df_sty_subset_disease.sample(n=disease_sample, random_state = 42)
random_mesh_disease_02 = mesh_large_df_sty_subset_disease.sample(n=disease_sample, random_state = 124)
random_mesh_01 = mesh_large_df_sty_subset.sample(n=sample_rest, random_state = 42)
random_mesh_02 = mesh_large_df_sty_subset.sample(n=sample_rest, random_state = 124)

random_mesh_01_con = pd.concat([random_mesh_disease_01, random_mesh_01])
random_mesh_02_con = pd.concat([random_mesh_disease_02, random_mesh_02])
random_mesh_02_con


# In[ ]:


if os.path.exists('ref-mesh-archive/random_mesh_01_con.csv'):
    pass
else:
    random_mesh_01_con.to_csv('ref-mesh-archive/random_mesh_01_con.csv')
    random_mesh_02_con.to_csv('ref-mesh-archive/random_mesh_02_con.csv')


# In[ ]:


random_mesh_02_con.groupby('Semantic Types Label').describe()
random_mesh_01_con.groupby('Semantic Types Label').describe()


# In[ ]:


disease_sample = 255
amino_sample = 92
biolog_sample = 65
new_nutri_refmesh_count


# # Analyze Reference Mesh list

# In[ ]:


import ast
ref = pd.read_csv('ref-mesh-archive/ref_mesh_nest.csv', sep=";")
print('nested ref mesh', ref['Class ID'].nunique())

ref['Semantic Types'] = ref['Semantic Types'].apply(ast.literal_eval)
#otherwise= dfg['col1'] = df['col1'].str.strip('[]').str.split(', ').apply(lambda x: [int(i) for i in x])
sty = pd.read_csv('ref-mesh-archive/MeshSTY-code.csv',sep=';')
sty = sty.rename(columns={'ID':'Semantic Types'})
ref


# In[ ]:


ref_large = []
for i in range(len(ref)):
    for sem in ref['Semantic Types'][i]: #dfrspost = mother table
        out = ref['Preferred Label'][i],ref['Class ID'][i],sem,ref['Synonyms'][i],ref['Parents'][i],ref['CUI'][i],ref['AQL'][i],ref['TERMUI'][i]
        ref_large.append(out)

ref_large_df = pd.DataFrame(ref_large)
new_col_names = ['Preferred Label','Class ID','Semantic Types','Synonyms','Parents','CUI','AQL','TERMUI']
ref_large_df.columns = new_col_names
ref_large_df
#ref_large_df.to_csv('ref-mesh-archive/ref_mesh_large.csv')


# ## Mesh semantics analyze on Ref list

# In[ ]:


## set df on ref mesh sty large

if os.path.exists('ref-mesh-archive/ref_mesh_sty_large.csv'):
    ref_sty = pd.read_csv('ref-mesh-archive/ref_mesh_sty_large.csv', index_col=0)
else:
    ref_sty = pd.merge(ref_large_df, sty, on='Semantic Types', how='inner').reset_index(drop=True)
    #Add rsid coulmn con merge
    #rspmidmesh_merge = pd.merge(pmidmesh, rsidpmid, on= 'pmids', how='inner').drop_duplicates().reindex(columns=['pmids', 'rsid', 'mesh'])
    ref_sty = ref_sty.rename(columns={'Preferred Label_y':'Semantic Types Label','Preferred Label_x':'Preferred Label'})
    #ref_sty = ref_sty.drop('Column1',axis=1)
    #ref_sty.to_csv('ref-mesh-archive/ref_mesh_sty_large.csv')

ref_sty = ref_sty[['Preferred Label', 'Semantic Types Label', 'Class ID', 'Semantic Types', 'Synonyms', 'Parents', 'CUI', 'AQL', 'TERMUI']]


# In[ ]:


print(ref_sty['Preferred Label'].nunique(), 'mesh total')
ref_sty


# In[ ]:


#modulo groupby and bar
ref_sty_less = ref_sty[['Preferred Label','Semantic Types Label']]

### groupby.describe analysis by rsid--------------------
ref_sty_less_count = ref_sty_less.groupby('Semantic Types Label').describe().reset_index()
ref_sty_less_count.columns = ref_sty_less_count.columns.to_flat_index()
new_column_names = ['Semantic Types Label', 'mesh-count', 'mesh-unique','mesh-top','mesh-freq']
ref_sty_less_count.columns = new_column_names
ref_sty_less_count_sort = ref_sty_less_count.sort_values(by='mesh-count',ascending=False).reset_index(drop=True)

x = ref_sty_less_count_sort['Semantic Types Label'].iloc[:]
y = ref_sty_less_count_sort['mesh-count'].iloc[:]
plt.figure(figsize=(4, len(ref_sty_less_count_sort)*0.25))
plt.title('Reference Mesh- Semantic Type enrichment', loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('Semantic Types')
plt.xlabel('mesh', position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
# use log scale
plt.gca().set_xscale('log')
#plt.savefig('Reference Mesh- Semantic Type enrichment.png',dpi=300, bbox_inches = "tight")
plt.show()


# In[ ]:


ref_sty_less_count_sort#['Semantic Types Label']


# In[ ]:


# analizing abundance
#count = [12, 21]
#remove_value = ref_sty_less_count_sort[ref_sty_less_count_sort['mesh-count'].isin(count)]['Semantic Types Label'].reset_index(drop=True).to_list()
remove_value = ref_sty_less_count_sort[ref_sty_less_count_sort['mesh-count']> 56]['Semantic Types Label'].reset_index(drop=True).to_list()
ref_sty[ref_sty['Semantic Types Label'].isin(remove_value)]['Semantic Types Label'].drop_duplicates()


# ## STY Ripper form complete ref mesh

# In[ ]:


## Rip form complete ref mesh
save = False
rip_list = ['Cell', 'Receptor', 'Hormone', 'Tissue', 'Body Part, Organ, or Organ Component', 'Congenital Abnormality', 'Anatomical Abnormality', 'Indicator, Reagent, or Diagnostic Aid']

remove_value = ref_sty_less_count_sort[ref_sty_less_count_sort['Semantic Types Label'].isin(rip_list)]['Semantic Types Label'].reset_index(drop=True).to_list()
remove_value = ref_sty_less_count_sort[ref_sty_less_count_sort['mesh-count']> 56]['Semantic Types Label'].reset_index(drop=True).to_list()

#new_ref_sty = ref_sty[ref_sty['Semantic Types Label'].isin(remove_value)]
mask = ref_sty['Semantic Types Label'].isin(rip_list)
new_ref_sty = ref_sty[-mask]

if save == True:
    new_ref_sty.to_csv('ref-mesh-archive/new_ref_mesh_large_corrected.csv')


# In[ ]:





# In[ ]:


new_ref_sty['Preferred Label'].drop_duplicates().to_csv('ref-mesh-archive/new_ref_mesh_corrected.csv')
new_ref_sty['Preferred Label'].drop_duplicates()


# # Global Mesh df break through:

# #### Search for all sons from a parent ID:
#     goal: merge all child mesh together

# In[ ]:


#Search for all brothers from a parent ID:
var4 = 'D000066888' #'Diet, Food, and Nutrition'
var5 = 'D044623' #'Nutrition Therapy'
var6 = 'D014808' #Nutritional and Metabolic Diseases
var7 = 'D002318'#Cardiovascular Diseases
var8 = 'D004066'#Digestive System Diseases
var9 = 'D004700'#Endocrine System Diseases
var10 = 'D006967'#Hypersensitivity

sub = df[df['Parents'].str.contains(var4).fillna(False)]
# use .dropna(inplace=True)   OR   df.fillna(0, inplace=True)
sub


# #### Concatenare i risultati
# 

# In[ ]:


dfgb = pd.DataFrame(columns=df.columns) #Empty DataFrame, same structure
#dfgg = pd.DataFrame()
dfgg = sub
# DEVO ricercare in 'Parents' l'ID della categoria superiore

#Parents Class ID list:
kids = sub['Class ID'].tolist() #search for nested childs

# Multiple search:
for i in kids:
    child = df[df['Parents'].str.contains(i).fillna(False)]
    dfgg = dfgg.append([child])
    dfgb = pd.concat([child, dfgb.loc[:]])
    kids = child['Class ID'].tolist()
    for i in kids:
        child = df[df['Parents'].str.contains(i).fillna(False)]
        dfgg = dfgg.append([child])
        #dfgb = pd.concat([child, dfgb.loc[:]])
        kids = child['Class ID'].tolist()
        for i in kids:
            child = df[df['Parents'].str.contains(i).fillna(False)]
            dfgg = dfgg.append([child])
            #dfgb = pd.concat([child, dfgb.loc[:]])
            kids = child['Class ID'].tolist()
            for i in kids:
                child = df[df['Parents'].str.contains(i).fillna(False)]
                dfgg = dfgg.append([child])
                #dfgb = pd.concat([child, dfgb.loc[:]])
                kids = child['Class ID'].tolist()
                for i in kids:
                    child = df[df['Parents'].str.contains(i).fillna(False)]
                    dfgg = dfgg.append([child])
                    #dfgb = pd.concat([child, dfgb.loc[:]])
                    kids = child['Class ID'].tolist()


# In[ ]:


dfgg.loc[:].drop_duplicates()


# In[ ]:


#dfgb.drop_duplicates()


# In[ ]:


chat_mesh_nutri#.drop_duplicates()


# In[ ]:


len(data)


# for i in data:
#     meshy = df[df['Synonyms'].str.contains(i).fillna(False)]
#     chatmesh = chatmesh.append([meshy])

# In[ ]:


#chatmesh.iloc[0]


# chatmeshsub = chatmesh[['Preferred Label','Class ID','Synonyms']].drop_duplicates()
# chatmeshsub.to_csv('ref-mesh-archive/df5_chatmeshsub_unique-2312mesh.csv')
# chatmeshsub

# In[ ]:


len(chatmeshsub2['Preferred Label'].drop_duplicates())


# In[ ]:





# dfgb = pd.DataFrame(columns=df.columns) #Empty DataFrame, same structure
# #fgb = pd.concat([z, dfgb.loc[:]])
# #DEVO  ricercare in 'Parents' l'ID della categoria superiore
# dfgb
# #Multiple search:
# for i in lisid:
#     dfgb = pd.concat([i, dfgb.loc[:]])
# dfgb.loc[:]
# 
# 
# 

# 
# varl1 = ['D004035', 'D003922']
# #serch for class ID
# df[df['Class ID'].str.contains(var3)]
# 
# dfgg = pd.DataFrame(columns=df.columns)
# for i in varl1:
#     #dfg.append(df[df['Class ID'].str.contains(i)], ignore_index=True)
#     row = df[df['Class ID'].str.contains(i)]
#     dbv= pd.concat([row, dfgg.loc[:]])
# dbv # dont' work

# In[ ]:


#Create new empty dataframe
dfg = pd.DataFrame(columns=df.columns)

ds1 = df[df['Class ID'].str.contains(var3)]
ds2 = df[df['Class ID'].str.contains(var4)]
#dfg.append(arr, ignore_index=True)
#dfg2 = dfg.append(ds1)
dfg2 = pd.concat([ds2, dfg2.loc[:]])
dfg2
#USE CONCAT!


# In[ ]:


# Importing pandas as pd
import pandas as pd

# Creating the first Dataframe using dictionary
df1 = df = pd.DataFrame({"a":[1, 2, 3, 4],
                         "b":[5, 6, 7, 8]})

# Creating the Second Dataframe using dictionary
df2 = pd.DataFrame({"a":[1, 2, 3],
                    "b":[5, 6, 7]})

# Print df1
print(df1, "\n")

# Print df2
print(df2)
# to append df2 at the end of df1 dataframe
df1.append(df2)


# In[ ]:


dft = pd.DataFrame()
# Create a new dataframe with the desired columns and values
new_data = pd.DataFrame({'Name': ['John', 'Jane'], 'Age': [25, 30]})

# Append the new dataframe to the existing one
dft = dft.append(new_data, ignore_index=True)

for i in varl1:
    dft.append(i, ignore_index=True)
dft


# In[ ]:





# In[ ]:


dd= df[['Preferred Label','MeSH Frequency']]
ddf= dd.dropna(subset='Preferred Label')
print(ddf.shape)
print(dd.shape)
import matplotlib.pyplot as plt

#plt.scatter(dd['Preferred Label'], dd['MeSH Frequency'])
ddf


# In[ ]:


dd


# In[ ]:


#df["Preferred Label"].__contains__("human")
df2=df[df["Preferred Label"].str.contains('Polymorphism')]


# In[ ]:


df3= df2[['Preferred Label','MeSH Frequency']]
df3


# In[ ]:


df2=df[df["Preferred Label"].str.contains('Diabetes')]
df3= df2[['Preferred Label','MeSH Frequency']]
df3


# In[ ]:


df2=df[df["Preferred Label"].str.contains('Fabry')]
df3= df2[['Preferred Label','MeSH Frequency']]
df3


# In[ ]:


#https://www.nlm.nih.gov/databases/download/data_distrib_main.html
#https://www.nlm.nih.gov/databases/download/mesh.html
#https://id.nlm.nih.gov/mesh/swagger/ui


# ## Confronto di specifiche MESH list con il MESH big data

# In[ ]:


dn = pd.read_csv('ref-mesh-archive/nbdb1-mesh.csv')

print(len(dn['mesh'].drop_duplicates()))
dga= df[df['Preferred Label'].isin(dn['mesh'])]

#dgb = dn[-df['Preferred Label'].isin(dn['mesh'])].drop_duplicates()

#see not contained in all-MESH
dgaa= dga['Preferred Label']
dgb = dn['mesh'].drop_duplicates()
dgg = dgb[dgb.isin(dgaa)==False]
dgg


# #### Interpolazione con mesh selected fprm semantic type

# In[ ]:


dn = pd.read_csv('ref-mesh-archive/MeshEX-meshsorted.csv',sep=';')
print(len(dn['Preferred Label'].drop_duplicates()))
dnmask = dn.Column2.isin(['green'])
dnfilter = dn[dnmask]
dnfilter
# moldulo 'Allineator'


# #### moldulo 'Allineator' for list value
# dnallign = []
# for i in range(len(dnfilter)):
#     for pmid in dnfilter['Semantic Types'][i]: #dfrspost = mother table
#         #out = (dfrspost['rsid'][i], pmid)
#         out = dnfilter['Preferred Label'][i], pmid
#         dnallign.append(out)
# print(type(dnallign))

# #### MODULO per table explosion:
# 
# data.csv = mesh; 'Semantic Types'; ID
# Bone Density,"T201,T081",D015519
# Ideal Body Weight,"T201,T074,T081",D056865ù
# 
# split the SEM column by comma and explode the values to create a new row for each SEM value
# df = df.assign(SEM=df['SEM'].str.split(',')).explode('SEM')

# In[ ]:


# MOLDULO per table explosion:

# split the SEM column by comma and explode the values to create a new row for each SEM value
dnfilter = dnfilter.assign(SEM=dnfilter['Semantic Types'].str.split(',')).explode('SEM')

# print the updated table
dnfilerstimple = dnfilter[['Preferred Label','Class ID','Semantic Types','SEM']]
dnfiltsupersimple = dnfilter[['Preferred Label','SEM']].drop_duplicates()
dnfiltsupersimple.to_clipboard()


# In[ ]:


#Modulo analisi groupby:
dnfilercount= dnfiltsupersimple.groupby('SEM').describe().reset_index()
dnfilercount.columns = dnfilercount.columns.to_flat_index()
new_column_names = ['SEM', 'Labes_count', 'Labes_count_unique','Labes_top','Labels_freq']
dnfilercount.columns = new_column_names
dnfilercountsort = dnfilercount.sort_values(by= 'Labes_count' ,ascending = False)
dnfilercountsortsmall= dnfilercountsort[['SEM','Labes_count']]


# In[ ]:


# Add semantic name column #Use Merge
semcode = pd.read_csv('ref-mesh-archive/MeshSTY-code.csv',sep=';')
new_column_names = ['Label', 'SEM']
semcode.columns = new_column_names
dnfilercountsortsmall2 = pd.merge(dnfilercountsortsmall, semcode, on='SEM', how='inner')
print(len(dnfilercountsortsmall.SEM))
dnfilercountsortsmall2


# import matplotlib.pyplot as plt
# x = dnfilercountsortsmall2.iloc[:,2]
# y = dnfilercountsortsmall2.iloc[:,1]
# #plt.figure(figsize=(15, 12))
# plt.figure(figsize=(12, 6))
# plt.scatter(x, y)
# 
# plt.title('Scatter Plot: "nutritional physiology" reference mesh plot')
# plt.xticks(rotation=90)
# 
# plt.show()
# 

# In[ ]:


import matplotlib.pyplot as plt
w = dnfilercountsortsmall2.Labes_count[dnfilercountsortsmall2['Labes_count'] >1]
u = dnfilercountsortsmall2.Label[dnfilercountsortsmall2['Labes_count'] >1]

#plt.figure(figsize=(15, 12))
plt.figure(figsize=(12, 6))
plt.scatter(u,w)

plt.title('Scatter Plot: "nutritional physiology" reference mesh plot')
plt.xticks(rotation=90)

plt.show()
print(len(dnfilercountsortsmall2.Label[dnfilercountsortsmall2['Labes_count'] >1]))


# In[ ]:


dnfilercountsortsmall2[dnfilercountsortsmall2['Labes_count'] >1].head(8)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# dga= df[df['Preferred Label'].isin(dnfilter['Preferred Label'])]
# 
# #dga['Preferred Label'].drop_duplicates()
# dga

# In[ ]:




