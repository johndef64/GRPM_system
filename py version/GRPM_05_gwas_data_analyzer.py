#!/usr/bin/env python
# coding: utf-8

# # Import modules

# In[ ]:


#Import Modules
import json
import sys
import nbib
import pandas as pd
import requests as rq
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pyperclip
from datetime import datetime
import os
from bs4 import BeautifulSoup

from Bio import Entrez
Entrez.email = "your_email@example.com"


# # Import dataset

# In[ ]:


#import dataset
df_gwas = pd.read_table('gwas_catalog_data/gwas_catalog_v1.0.2-associations_e109_r2023-03-27.tsv')
df_gwas[['PUBMEDID','SNP_ID_CURRENT']] = df_gwas[['PUBMEDID','SNP_ID_CURRENT']].astype(str)
print('genes: ', len(df_gwas['MAPPED_GENE'].drop_duplicates()))
print('studies: ', len(df_gwas['STUDY'].drop_duplicates()))
print('rsid: ', len(df_gwas['SNP_ID_CURRENT'].drop_duplicates()))
df_gwas.columns


# In[ ]:


#df_gwas


# In[ ]:


#display selected columns
df_gwas['STRONGEST SNP-RISK ALLELE'].drop_duplicates()
df_gwas[['MAPPED_GENE','DISEASE/TRAIT','MAPPED_TRAIT','SNP_ID_CURRENT','STRONGEST SNP-RISK ALLELE','RISK ALLELE FREQUENCY']].drop_duplicates()


# In[ ]:


print('SNPs:',len(df_gwas.SNPS.drop_duplicates()))
print('DISEASE/TRAIT:',len(df_gwas['DISEASE/TRAIT'].drop_duplicates()))
print('MAPPED_TRAIT:',len(df_gwas['MAPPED_TRAIT'].drop_duplicates()))
df_gwas['MAPPED_TRAIT'].value_counts()


# # Filter dataset

#  #Workflow------------------------
# 
# in STRONGEST SNP-RISK ALLEL drop "?"
# filter for your rsids
# --> I need the Merged data from GRPM-X:
#       - [gene, rsids, pmids, mesh(?)]
# merge tables on rsids
# comparare GENE-MESH vs
# get the STRONGEST SNP-RISK ALLELE
# 
# Column Selection Mode = Alt+Maiusc+Ins
# #https://www.jetbrains.com/help/rider/Multicursor.html#add-and-remove-carets

# In[ ]:


# Drop no risk allele:
mask = df_gwas['STRONGEST SNP-RISK ALLELE'].str.contains("\?")
df_gwas_drop = df_gwas[-mask].reset_index(drop=True)

# Drop complementary base allele (risk allele freq missing)
df_gwas_drop_nonan = df_gwas_drop.dropna(subset=['RISK ALLELE FREQUENCY'],axis=0).reset_index(drop=True)
print('Drop no risk allele:')
print('SNPs:',len(df_gwas_drop_nonan.SNPS.drop_duplicates()))
print('DISEASE/TRAIT:',len(df_gwas_drop_nonan['DISEASE/TRAIT'].drop_duplicates()))
print('MAPPED_TRAIT:',len(df_gwas_drop_nonan['MAPPED_TRAIT'].drop_duplicates()))


# Full dataset:
# SNPs: 267372
# DISEASE/TRAIT: 21399
# MAPPED_TRAIT: 7690

# ## LOOKUP FOR SINGLE RSID

# In[ ]:


# LOOKUP FOR SINGLE RSID
rsid_mask = df_gwas_drop_nonan['SNPS'].str.contains('rs1421085')
df_gwas_drop_nonan_rsid = df_gwas_drop_nonan[rsid_mask]
df_gwas_drop_nonan_rsid[['MAPPED_GENE','DISEASE/TRAIT','MAPPED_TRAIT','SNP_ID_CURRENT','STRONGEST SNP-RISK ALLELE','RISK ALLELE FREQUENCY']].drop_duplicates()
df_gwas_drop_nonan_rsid.value_counts('DISEASE/TRAIT')
df_gwas_drop_nonan_rsid.value_counts('STRONGEST SNP-RISK ALLELE')


# In[ ]:


# Display and to clipboard
df_gwas_drop_nonan[['MAPPED_GENE', 'DISEASE/TRAIT', 'SNPS','STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()#to_clipboard()


# In[ ]:


df_gwas_drop_nonan_rsid[['MAPPED_GENE', 'DISEASE/TRAIT', 'SNPS','STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()#.to_clipboard()


# # Merge GWAS and GRPMX data

# In[ ]:


# load my GRPMx Data
directory = 'grpm_survey_pcg_'+'nutri'
df_grpmx = pd.read_csv(directory+'/grpmx_filtered_output.csv', index_col=0)
df_grpmx_repo = pd.read_csv(directory+'/GRPMX_report_int.csv')

#filter for 0.95 quantile
df_grpmx_repo_95 = df_grpmx_repo[df_grpmx_repo.interest_index >= df_grpmx_repo.interest_index.quantile(0.95)]
df_grpmx_95 = df_grpmx[df_grpmx.gene.isin(df_grpmx_repo_95.gene)]
print('df_grpmx_95 genes:', df_grpmx_95.gene.nunique())
df_grpmx_95


# In[ ]:


#threshold
df_grpmx_t = df_grpmx_95

## ADD gene-interest index as common sorting handle
small_dummy = df_grpmx_repo[['gene','interest_index']]
df_grpmx_int = pd.merge(df_grpmx, small_dummy, left_on='gene', right_on='gene')
df_grpmx_95_int = pd.merge(df_grpmx_95, small_dummy, left_on='gene', right_on='gene')

print('GRPMX threshold Statistics:')
print('genes:', df_grpmx_t.gene.nunique())
print('rsid:', df_grpmx_t.rsid.nunique())
print('mesh:', df_grpmx_t.mesh.nunique())
df_grpmx_int


# In[ ]:


complete_merge = False
if complete_merge == True:
    #common handle sort
    df_grpmx_int = df_grpmx_int.sort_values(by=['interest_index','rsid','mesh'], ascending =False).reset_index(drop=True)
else:
    df_grpmx_95_int = df_grpmx_95_int.sort_values(by=['interest_index','rsid','mesh'], ascending =False).reset_index(drop=True)


# In[ ]:


if complete_merge == True:
    # Merge two df on rsid:
    df_merged = pd.merge(df_grpmx_int, df_gwas, left_on='rsid', right_on='SNPS')
    df_merged[['pmids','PUBMEDID']] = df_merged[['pmids','PUBMEDID']].astype(str)

    df_merged_drop = pd.merge(df_grpmx_int, df_gwas_drop_nonan, left_on='rsid', right_on='SNPS')
    df_merged_drop[['pmids','PUBMEDID']] = df_merged_drop[['pmids','PUBMEDID']].astype(str)


# In[ ]:


if complete_merge == True:
    print('complete grpmx merge')
    print(len(df_merged_drop), len(df_merged))
    df_merged_drop[['gene','rsid','mesh','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE']].drop_duplicates()
    print(df_merged_drop.rsid.nunique(), 'rsid')
    print(df_merged_drop.mesh.nunique(), 'mesh')
    print(df_merged_drop['DISEASE/TRAIT'].nunique(), 'DISEASE/TRAIT')
    print(df_merged_drop['STRONGEST SNP-RISK ALLELE'].nunique(), 'STRONGEST SNP-RISK ALLELE')

    df_merged_drop[['gene','rsid','mesh','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE']].drop_duplicates()


# In[ ]:


merge_nonan = True
merge_both = True

if complete_merge == False:
    # Merge two df on rsid:
    if merge_nonan == False:
        df_merged_95 = pd.merge(df_grpmx_95_int, df_gwas, left_on='rsid', right_on='SNPS')
        df_merged_95[['pmids','PUBMEDID']] = df_merged_95[['pmids','PUBMEDID']].astype(str)
    if merge_nonan == True:
        df_merged_95_drop = pd.merge(df_grpmx_95_int, df_gwas_drop_nonan, left_on='rsid', right_on='SNPS')
        df_merged_95_drop[['pmids','PUBMEDID']] = df_merged_95_drop[['pmids','PUBMEDID']].astype(str)
    if merge_both == True:
        df_merged_95 = pd.merge(df_grpmx_95_int, df_gwas, left_on='rsid', right_on='SNPS')
        df_merged_95[['pmids','PUBMEDID']] = df_merged_95[['pmids','PUBMEDID']].astype(str)
print('genes merged:', df_merged_95_drop.gene.nunique())
df_merged_95_drop


# In[ ]:


print('threshold (q0.95) grpmx merge')
print('full 95 merge pmids:',df_merged_95.pmids.nunique(), 'drop95 pmids:',df_merged_95_drop.pmids.nunique())
#df_merged_95_drop[['gene','rsid','mesh','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE']].drop_duplicates()
print('grpm merged with gwas stats:')
print(' ',df_merged_95_drop.gene.nunique(), 'gene')
print(' ',df_merged_95_drop.pmids.nunique(), 'pmid')
print(' ',df_merged_95_drop.mesh.nunique(), 'mesh')
print(' ',df_merged_95_drop.rsid.nunique(), 'rsid')
print(' ',df_merged_95_drop['DISEASE/TRAIT'].nunique(), 'DISEASE/TRAIT')
print(' ',df_merged_95_drop['STRONGEST SNP-RISK ALLELE'].nunique(), 'STRONGEST SNP-RISK ALLELE')


# In[ ]:


df_merged_less = df_merged.drop(columns=['DATE ADDED TO CATALOG','MERGED','SNP_ID_CURRENT'])
len(df_merged_less['STUDY'].drop_duplicates())
print('GWAS TRAIT count con entire grpmx')
df_merged_less.value_counts('DISEASE/TRAIT')


# In[ ]:


print('GWAS TRAIT count (nonan) con entire grpmx')
df_merged_drop_less = df_merged_drop.drop(columns=['DATE ADDED TO CATALOG','MERGED','SNP_ID_CURRENT'])
len(df_merged_less['STUDY'].drop_duplicates())
df_merged_drop_less.value_counts('DISEASE/TRAIT')


# In[ ]:


#rename columns:
df_merged_drop = df_merged_drop.rename(columns={'gene':'LITVAR GENE', 'rsid':'LIVAR RSID', 'pmids':'LITVAR PMID','mesh':'PUBMED MESH'})
df_merged_drop_nonan = df_merged_drop


# In[ ]:


#df_merged_drop_nonan.to_csv('df_merged_drop_nonan.csv') # diversi GB!!!!!!!!!!!!!!!!
df_usage = df_merged.memory_usage()/1048576
df_usage.sum()


# In[ ]:


df_merged_drop[['LIVAR RSID','LITVAR GENE']].drop_duplicates()


# #### Define Mother Dataframe:

# In[ ]:


df_merged_drop


# In[ ]:


# LOOKUP FOR SINGLE RSID
rsid_mask = df_merged_drop_nonan['LIVAR RSID'].str.contains('rs1421085')
df_merged_drop_nonan_rsid = df_merged_drop_nonan[rsid_mask]


# In[ ]:


# Display and to clipboard
df_merged_drop[['LITVAR GENE','LIVAR RSID','MAPPED_GENE','PUBMED MESH', 'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()
df_merged_drop_nonan_rsid[['LITVAR GENE','LIVAR RSID','MAPPED_GENE','PUBMED MESH', 'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()#.to_clipboard()


# In[ ]:


df_merged_drop_nonan_rsid['PUBMED MESH'].drop_duplicates()


# #### Get risk Allele list and use it to filter mother table

# In[ ]:


df_merged_drop_nonan.value_counts('SNPS')


# In[ ]:


#df_merged_drop_rsid_nonan[['STRONGEST SNP-RISK ALLELE','PUBMED MESH']].groupby('STRONGEST SNP-RISK ALLELE').describe().reset_index()
rsid_mask = df_merged_drop_nonan['LIVAR RSID'].str.contains('rs1421085')
df_merged_drop_nonan_rsid = df_merged_drop_nonan[rsid_mask]
df_merged_drop_nonan_rsid[['LITVAR GENE','LIVAR RSID','MAPPED_GENE','PUBMED MESH', 'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()
type(df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts())#.head(1))
risk_allele = df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts()#.index[0]
risk_allele


# In[ ]:


#---> creare una lista programmaticamente di tutti i 'risk allele' by count and use it to filter mother dataframe with isin module!

# get all rsid list
rsid_list = df_merged_drop_nonan['LIVAR RSID'].drop_duplicates().to_list()
len(rsid_list)


# In[ ]:


# risk allele pickup (part1)
time_start = datetime.now()
risk_allele_list = []
for i in rsid_list[:2000]:
    rsid_mask = df_merged_drop_nonan['LIVAR RSID'].str.contains(i)
    df_merged_drop_nonan_rsid = df_merged_drop_nonan[rsid_mask]
    risk_allele = df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts().index[0]
    risk_allele_list.append(risk_allele)
    #print(str(risk_allele))
finish_start = datetime.now()
pd.Series(risk_allele_list).to_csv('risk_allele_list_0-2000.csv')


# In[ ]:


# risk allele pickup (part2)
for i in rsid_list[2000:]:
    rsid_mask = df_merged_drop_nonan['LIVAR RSID'].str.contains(i)
    df_merged_drop_nonan_rsid = df_merged_drop_nonan[rsid_mask]
    risk_allele = df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts().index[0]
    risk_allele_list.append(risk_allele)
    #print(str(risk_allele))
finish_start = datetime.now()
print(finish_start - time_start)
pd.Series(risk_allele_list).to_csv('risk_allele_list_0-2000.csv')


# In[ ]:


# import back risk allele list

risk_allele_df = pd.read_csv('risk_allele_list_4376.csv', index_col=0)
risk_allele_list = risk_allele_df['0'].to_list()
risk_allele_list
# --> ora filtrare Mother Df per i risk allele and ...BAM!
risk_allele_df


# In[ ]:





# In[ ]:


# genera la merged  with GRPMX dropped!
#df_merged_drop_nonan = pd.read_csv('df_merged_drop_nonan.csv', index_col=0)

df_merged_drop_less_gwa = df_merged_drop_nonan[['MAPPED_GENE','SNPS', 'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()
#df_merged_drop_less_gwa.to_csv('df_merged_drop_less_gwa.csv')#.SNPS.drop_duplicates()
df_merged_drop_less_gwa.SNPS.drop_duplicates().sample(10)


# In[ ]:


type(df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts())#.head(1))
risk_allele = df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].value_counts()#.index[0]
risk_allele


# In[ ]:


#qualti sono monorischio??
rsid_list = df_merged_drop_less_gwa.SNPS.drop_duplicates().to_list()
time_start = datetime.now()
monorisk_list = []

for i in rsid_list:
    rsid_mask = df_merged_drop_less_gwa['SNPS'].str.contains(i)
    df_merged_drop_less_gwa_rsid = df_merged_drop_less_gwa[rsid_mask].drop_duplicates()
    risk_allele = df_merged_drop_less_gwa_rsid['STRONGEST SNP-RISK ALLELE'].value_counts()#.index[0]
    if len(risk_allele)==1:
            monorisk = risk_allele.index[0]
            monorisk_list.append(monorisk)
            print(risk_allele)


# In[ ]:


#print(risk_allele)
pd.Series(monorisk_list)#.to_csv('monorisk_list.csv')
#len(rsid_list)


# In[ ]:


# risk allele pickup: ricerca dei valori di conteggio ambigui------------------

rsid_list = df_merged_drop_less_gwa.SNPS.drop_duplicates().to_list()
time_start = datetime.now()
ambiguity_list = []

for i in rsid_list:
    rsid_mask = df_merged_drop_less_gwa['SNPS'].str.contains(i)
    df_merged_drop_less_gwa_rsid = df_merged_drop_less_gwa[rsid_mask].drop_duplicates()
    risk_allele = df_merged_drop_less_gwa_rsid['STRONGEST SNP-RISK ALLELE'].value_counts()#.index[0]
    if len(risk_allele)>1:
        if risk_allele[0] == risk_allele[1]:
            ambiguity = risk_allele.index[0]
            ambiguity_list.append(ambiguity)
            print(risk_allele)

print(risk_allele)
ambiguity_list


# In[ ]:


ambiguity_list = []
for i in rsid_list:
    rsid_mask = df_merged_drop_less_gwa['SNPS'].str.contains(i)
    df_merged_drop_less_gwa_rsid = df_merged_drop_less_gwa[rsid_mask].drop_duplicates()
    risk_allele = df_merged_drop_less_gwa_rsid['STRONGEST SNP-RISK ALLELE'].value_counts()#.index[0]
    if len(risk_allele)>1:
        if risk_allele[0] == risk_allele[1]:
            ambiguity = risk_allele.index[0], risk_allele[0]
            ambiguity_list.append(ambiguity)
            #print(risk_allele)
ambiguity_df = pd.DataFrame(ambiguity_list)


# In[ ]:


#ambiguity_df.to_csv('ambiguity_magg1_df.csv')
pd.read_csv('ambiguity_1_list.csv')


# In[ ]:


ambiguity_df


# In[ ]:


ambiguity_df.groupby(by=1).describe().to_clipboard()


# ### remove ambiguity

# In[ ]:


risk_allele_df#[0]
df_merged_drop_less_gwa
type(risk_allele_df.iloc[:,0])


# In[ ]:


# select just ambiugous count > 1
ambiguous_magg1 = pd.read_csv('ambiguity_magg1_df.csv', index_col=0)
# Only GWAS
mask = df_merged_drop_less_gwa['SNPS'].isin(ambiguous_magg1.rsid)
df_merged_drop_less_gwa_ambmagg1 = df_merged_drop_less_gwa[mask]
df_merged_drop_less_gwa_ambmagg1#.to_csv('df_merged_drop_less_gwa_ambmagg1.csv')


# In[ ]:


#GRPMX-GWAS
mask = df_merged_drop_nonan['SNPS'].isin(ambiguous_magg1.rsid)
df_merged_drop_nonan_ambmagg1 = df_merged_drop_nonan[mask]
df_merged_drop_nonan_ambmagg1[['LITVAR GENE','PUBMED MESH','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE','RISK ALLELE FREQUENCY','MAPPED_GENE']]#.to_csv('df_merged_drop_less_gwa_ambmagg1.csv')


# In[ ]:


df_merged_drop_nonan[['LITVAR GENE','PUBMED MESH','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE','RISK ALLELE FREQUENCY','MAPPED_GENE']]
#quanti sono monorisk allele???


# In[ ]:


# before filter mother df with "riskallelelist

# Only GWAS
mask = df_merged_drop_less_gwa['STRONGEST SNP-RISK ALLELE'].isin(risk_allele_df.iloc[:,0])
df_merged_drop_less_gwa_riskall = df_merged_drop_less_gwa[mask]
#df_merged_drop_less_gwa_riskall.to_csv(r'C:\Users\Public\Database_SSD\GWAS calalog data\df_merged_drop_less_gwa_riskall.csv')
df_merged_drop_less_gwa_riskall


# In[ ]:


# GRPMX-GWAS merged
mask2 = df_merged_drop_nonan['STRONGEST SNP-RISK ALLELE'].isin(risk_allele_df.iloc[:,0])
df_merged_drop_nonan_riskall = df_merged_drop_nonan[mask2]
df_merged_drop_nonan_riskall.to_csv(r'C:\Users\Public\Database_SSD\GWAS calalog data\df_merged_drop_nonan_riskall.csv')


# In[ ]:


# then remove ambiugous rsids

# GRPMX-GWAS merged
ambiguous_rsids = pd.read_csv('ambiguous_rsids.csv')
mask_amb = df_merged_drop_nonan_riskall.SNPS.isin(ambiguous_rsids['ambiguous_rsids'])
df_merged_drop_nonan_riskall_unamb = df_merged_drop_nonan_riskall[-mask_amb]


# In[ ]:


df_merged_drop_nonan_riskall_unamb.columns


# In[ ]:


df_merged_drop_nonan_riskall_unamb[['LITVAR GENE', 'LIVAR RSID', 'LITVAR PMID', 'PUBMED MESH', 'PUBMEDID', 'DISEASE/TRAIT','MAPPED_GENE','STRONGEST SNP-RISK ALLELE', 'SNPS','RISK ALLELE FREQUENCY', 'P-VALUE','OR or BETA']][300:325].to_clipboard(sep=',')


# In[ ]:


# Only GWAS
mask_amb = df_merged_drop_less_gwa_riskall.SNPS.isin(ambiguous_rsids['ambiguous_rsids'])
df_merged_drop_less_gwa_riskall_unamb = df_merged_drop_less_gwa_riskall[-mask_amb]
df_merged_drop_less_gwa_riskall_unamb.to_csv(r'C:\Users\Public\Database_SSD\GWAS calalog data\df_merged_drop_less_gwa_unambiguous.csv')


# In[ ]:


df_merged_drop_less_gwa_riskall_unamb[300:325].to_clipboard(sep=',')


# In[ ]:


len(df_merged_drop_less_gwa_riskall_unamb.SNPS.drop_duplicates()), len(df_merged_drop_nonan.SNPS.drop_duplicates())
unambiguous_rsids = df_merged_drop_less_gwa_riskall_unamb.SNPS.drop_duplicates().reset_index(drop=True)
#unambiguous_rsids.to_csv('unambiguous_rsids.csv')
#df_merged_drop_nonan_unamb.to_csv(r'C:\Users\Public\Database_SSD\GWAS calalog data\df_merged_drop_nonan_unambiguous.csv')


# In[ ]:


df_merged_drop_less_gwa_unamb


# In[ ]:


# LOOKUP FOR SINGLE RSID
rsid_mask = df_merged_drop_less_gwa['SNPS'].str.contains('rs1558902')
df_merged_drop_less_gwa_rsid = df_merged_drop_less_gwa[rsid_mask]
df_merged_drop_less_gwa_rsid.to_clipboard()


# In[ ]:


# LOOKUP FOR SINGLE ALLELE
allele_mask = df_merged_drop_nonan_rsid['STRONGEST SNP-RISK ALLELE'].str.contains('-C')
df_merged_drop_rsid_allele = df_merged_drop_nonan_rsid[allele_mask]
df_merged_drop_rsid_allele[['LITVAR GENE','LIVAR RSID','MAPPED_GENE','PUBMED MESH', 'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE', 'RISK ALLELE FREQUENCY']].drop_duplicates()


# In[ ]:


dfx = pd.read_csv('grpm_survey/'+gene2+'_filtered_pmidrsidmesh.csv').drop(columns=['Unnamed: 0'])
gene2_pmidmesh = pd.read_csv('grpm_db/'+gene2+'_pubmed_pmidmesh.csv').drop(columns=['Unnamed: 0'])
dfy = gene2_pmidmesh[gene2_pmidmesh.pmids.isin(dfx.pmids)]

dfxy = pd.merge(dfx, dfy, left_on= '')

