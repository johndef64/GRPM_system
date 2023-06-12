#!/usr/bin/env python
# coding: utf-8

# # Import Modules

# In[ ]:


import os
import glob
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import sys
import numpy as np
#import pyperclip


# In[ ]:


#Only for for Google Colab
import sys
if 'google.colab' in sys.modules:
    from google.colab import drive
    drive.mount('/content/drive')
    os.chdir('/content/drive/MyDrive/GRPM db/GRPM db - Vitam')


# # Define context
#     - gene list
#     - survey directory
#     - ref-mesh list

# In[ ]:


#Check avalable refs:
ref_path = "ref-mesh-archive/"  # Replace with the actual ref mesh path

#---------------------------------
# use random mesh list?
random_mesh = True
if random_mesh == True:
    ref_path = "ref-mesh-archive/random_lists/"
#---------------------------------

# Create a file path pattern to match CSV files
file_pattern = os.path.join(ref_path, "*.csv")

# Use glob to get a list of file paths matching the pattern
csv_files = glob.glob(file_pattern)
csv_files_name = []
# Print the list of CSV files
for file in csv_files:
    file_name = os.path.basename(file)
    csv_files_name.append(file_name)

pd.set_option('display.max_rows', 100)
print('Available reference mesh lists:')
csv_files_df = pd.Series(csv_files_name)

if random_mesh == False:
    print(csv_files_df[csv_files_df.str.contains('ref_mesh_')].str.replace('ref_mesh_','').reset_index(drop=True))
else:
    print(csv_files_df[csv_files_df.str.contains('random')].str.replace('ref_mesh_','').reset_index(drop=True))


# history:
# redo= ['xeno','oxi_stress','intol','eat_taste','cvd','dmt2_ms','ob_bmi_less']
# pd.DataFrame(redo).sort_values(by=0)

# In[ ]:


#------------------------------------------------------
# define directory folder:
survey_path = '' # keep black for root path
topic_tag = 'random_grpm_09' # choose ref_mesh.csv tab from available
add = ''    # additional survey tag
db_tag = 'pcg'       # choose db
    # 'pcg'    = protein coding genes = grpm_db
    # 'rna'    = rna genes            = grpm_db_rna
    # 'pseudo' = pseudogenes          = grpm_db_pseudo

#------------------------------------------------------
# saving options:
save_plot = False
checkpoint = 200 #save data each x genes

#------------------------------------------------------

# Set directory folder name
if random_mesh == True:
    survey_path = 'grpm_random/'
    directory = survey_path + 'grpm_'+ db_tag  + '_' + topic_tag + add
else:
    directory = survey_path + 'grpm_survey_' + db_tag + '_' + topic_tag + add

# Import Mesh-reference list
if random_mesh == False:
    ref = pd.read_csv(ref_path + "ref_mesh_" + topic_tag + ".csv",index_col=0)
else:
    ref = pd.read_csv(ref_path + topic_tag + ".csv",index_col=0)

if 'mesh' in ref.columns:
    pass
else:
    ref = ref.rename(columns={'Preferred Label': 'mesh'})
ref_mesh_n = ref.mesh.nunique()
ref_mesh_list = ref['mesh'].drop_duplicates()

# create directory and load/create master dataframes-----
if not os.path.exists(directory):
    os.makedirs(directory)

if os.path.isfile(directory+'/grpmx_filtered_output.csv'):
    complete_df = pd.read_csv(directory+'/grpmx_filtered_output.csv',index_col=0)
else:
    complete_df = pd.DataFrame()

if os.path.isfile(directory+'/GRPMX_report.csv'):
    df_report_complete = pd.read_csv(directory+'/GRPMX_report.csv',index_col=0)
    restart = True
else:
    df_report_complete = pd.DataFrame()
    restart = False
#------------------------------------------------------
ref_mesh_list


# In[ ]:


# check checkpoint report:
df_report_complete.T


# # Load GRPM data from database

# In[ ]:


#Load GRPM db Report-----------------------------------------
db_name = 'grpm_db_'+db_tag

#get gene list from grpm report
GRPM_report = pd.read_csv(db_name+'/GRPM_report.csv',index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
grpm_genes_list = GRPM_report.gene.to_list()

#Import grpm data back-------------------------------------------
time_load_1 = datetime.now()

columns = ['gene', 'rsid', 'pmids', 'mesh']
dummy = pd.read_csv(db_name+'/grpm_table_output.csv', usecols=columns)

dummy['pmids'] = dummy['pmids'].astype(str) #convert pmid type in str
time_load_2 = datetime.now()
print(time_load_2-time_load_1)


# # Survey module

# In[ ]:


time_start = datetime.now()

if restart == True:
    restart_from = len(df_report_complete.T)
    gene_start = restart_from
    print('search restarted from '+str(restart_from))
else:
    gene_start = 0

# define gene list
genes = grpm_genes_list[gene_start:]

for gene in genes:

    time_alpha = datetime.now()
    timestamp = time_alpha.strftime('%Y%m%d%H%M%S')

    if gene in dummy.gene.drop_duplicates().to_list():

        dummy_gene = dummy.loc[dummy['gene'] == gene]
        rsidpmid = dummy_gene[['rsid','pmids']].drop_duplicates().reset_index(drop=True)
        pmidmesh = dummy_gene[['pmids','mesh']].drop_duplicates()
        #dfmesh = dummy_gene[['pmids', 'mesh', 'qualifier', 'major']].drop_duplicates().reset_index(drop=True)

        #Filter pmid for rsid with pmid>1
     #   rsidpmid_count = rsidpmid.groupby('rsid').describe().reset_index()
     #   rsidpmid_count.columns = rsidpmid_count.columns.to_flat_index()
     #   new_column_names = ['rsid', 'pmid_count', 'pmid_unique','pmid_top','pmid_freq']
     #   rsidpmid_count.columns = new_column_names
     #   outless = rsidpmid_count[rsidpmid_count.pmid_unique>1]
     #   mask = rsidpmid['rsid'].isin(outless.rsid)
     #   rsidpmid_less = rsidpmid[mask]
     #   #repo parametes
        lit1_rsid = dummy_gene.rsid.nunique()
     #   lit1_pmid_f = rsidpmid_less.pmids.nunique()

        # Correlazione su "pmidmesh.mesh"---------------------------------------
        mask = pmidmesh['mesh'].isin(ref_mesh_list)
        dfmatch = pmidmesh[mask]
        mask_full = dummy_gene['mesh'].isin(ref_mesh_list)
        dfmatch_full = dummy_gene[mask_full]
        #Statistics:
        pmidmesh_before =  pmidmesh.nunique()
        pmidmesh_after =    dfmatch.nunique()
        interesting_pmid =  dfmatch.nunique()

        #pmidmask = rsidpmid['pmids'].isin(dfmatch.pmids) #mymask
        #rsidlast = rsidpmid[pmidmask]  # mask on rsidpmid
        #rsidlastlist = rsidlast.rsid.drop_duplicates()

        #report parameters
        matching_rsid = dfmatch_full['rsid'].nunique()
        dropped_rsid = lit1_rsid - dfmatch_full['rsid'].nunique()

        starting_pmid = pmidmesh['pmids'].nunique()
        starting_mesh = pmidmesh['mesh'].nunique()
        starting_pmidmesh = len(pmidmesh)
        matching_pmids = dfmatch.pmids.nunique()
        matching_mesh = dfmatch.mesh.nunique()
        matching_pmidmesh = len(dfmatch)

        #dfmatch = pmidmesh[mask]
      #  mask = pmidmesh['mesh'].isin(ref['mesh'])
      #  pmidmeshint = pmidmesh.assign(interesting=mask)

        #Analyze enrichment with groupby.describe method--------------------
        dfmatch_less_ = dfmatch_full[['pmids', 'rsid', 'mesh']].drop_duplicates()
      #  interesting_rsid = dfmatch_less_.rsid.nunique()

        #dfmatch_less_['pmids'] = dfmatch_less_['pmids'].astype(str)

        ### groupby.describe analysis by rsid
        dfmatch_less_rsid = dfmatch_less_.groupby('rsid').describe().reset_index()
        dfmatch_less_rsid.columns = dfmatch_less_rsid.columns.to_flat_index()
        new_column_names = ['rsid', 'pmid-count', 'pmid-unique','pmid-top','pmid-freq','mesh-count', 'mesh-unique','mesh-top','mesh-freq']
        dfmatch_less_rsid.columns = new_column_names

        #statistics:
        matching_rsid_pmid10 = len(dfmatch_less_rsid[dfmatch_less_rsid['pmid-unique']>10])
        matching_rsid_pmid100 = len(dfmatch_less_rsid[dfmatch_less_rsid['pmid-unique']>100])

        #sortng, top10
        dfmatch_less_rsidless = dfmatch_less_rsid[['rsid','pmid-unique','mesh-unique']]
        dfmatch_less_rsidlesssort = dfmatch_less_rsidless.sort_values(by='pmid-unique', ascending= False).reset_index(drop=True)
        top10rsid = dfmatch_less_rsidlesssort['rsid'][:10].tolist()

        ### groupby.describe analysis by mesh
        dfmatch_less_mesh = dfmatch_less_.groupby('mesh').describe().reset_index()
        dfmatch_less_mesh.columns = dfmatch_less_mesh.columns.to_flat_index()
        #to handle generate df.groupby.describe, convert Multicolumn to single column
        #https://datascientyst.com/flatten-multiindex-in-pandas/
        new_column_names = ['mesh', 'pmid-count', 'pmid-unique','pmid-top','pmid-freq','rsid-count', 'rsid-unique','rsid-top','rsid-freq']
        dfmatch_less_mesh.columns = new_column_names

        dfmatch_less_mesh_less = dfmatch_less_mesh[['mesh','pmid-unique','rsid-unique']]
        #dfmatch_less_mesh_lesssort = dfmatch_less_mesh_less.sort_values(by='pmid-unique',ascending=False).reset_index(drop=True)

        # add frequency
        samplepmid_count = len(dfmatch.pmids.drop_duplicates())
        dfmatch_less_mesh_less_frq = dfmatch_less_mesh_less.copy()
        mesh_frq = dfmatch_less_mesh_less_frq.loc[:,'pmid-unique'].astype(float)/samplepmid_count
        dfmatch_less_mesh_less_frq.loc[:,'mesh frequency'] = round(mesh_frq,3)#*100
        dfmatch_less_mesh_less_frqsort = dfmatch_less_mesh_less_frq.sort_values(by='pmid-unique',ascending=False).reset_index(drop=True)

        top10mesh = dfmatch_less_mesh_less_frqsort['mesh'][:10].tolist()
        #^display(dfmatch_less_mesh_less_frqsort.head(20))

        if save_plot == True:
            # create a scatter plot
            x = dfmatch_less_mesh_less_frqsort['mesh'].head(30)
            y = dfmatch_less_mesh_less_frqsort['pmid-unique'].head(30)
            plt.figure(figsize=(5, 8))
            plt.title('Scatter Plot: '+gene+' pmid-mesh (filtered)', loc='center',pad=10)
            plt.scatter(y, x)
            plt.gca().invert_yaxis()
            #plt.subplots_adjust(left=0.3, right=0.9, bottom=0.3, top=0.9)
            #plt.xticks(rotation=90)
            plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
            plt.xlabel('pmid count', position=(0.5, 1.08))
            ax = plt.gca()
            ax.xaxis.set_label_position('top')
            #plt.show()
            plt.savefig(directory+'/'+gene+'_mesh_plot_'+timestamp+'_filtered.png',dpi=120, bbox_inches = "tight")
            plt.close()
        else:
            pass

        #STORE DATA----------------------------------------------------------------------
        #timestamp = time2.strftime('%Y%m%d%H%M%S')

        #comparison results:
        dfmatch_less_['gene'] = gene
        complete_df = pd.concat([complete_df, dfmatch_less_])


        #REPORT------------------------------------------------------------------

        report = { 'reference_mesh': ref_mesh_n,
                   #'filtered_pmidmesh': pmidmesh_after,
                   #'interesting_pmid': interesting_pmid,
                   #'interesting_rsid': interesting_rsid,
                   'starting_pmidmesh': starting_pmidmesh,
                   'starting_pmid' : starting_pmid,
                   'starting_mesh': starting_mesh,
                   'starting_rsid': lit1_rsid,
                   'matching_pmidmesh': matching_pmidmesh,
                   'matching_pmids': matching_pmids,
                   'matching_mesh': matching_mesh,
                   'matching_rsid': matching_rsid,
                   'dropped_rsid': dropped_rsid,
                   'matching_mesh_ratio': round((matching_mesh/starting_mesh),3),
                   'matching_pmids_ratio': round((matching_pmids/starting_pmid),3),
                   'matching_pmidmesh_ratio': round((matching_pmidmesh/starting_pmidmesh),3),
                   'matching_rsid_ratio': round((matching_rsid/lit1_rsid),3),
                   'matching_rsid_pmid10': matching_rsid_pmid10,
                   'matching_rsid_pmid100': matching_rsid_pmid100,
                   'matching_top10mesh':str(top10mesh),
                   'matching_top10rsid':str(top10rsid),
                   }

        df_report = pd.DataFrame(report, index=[gene]).transpose()

        # SLOW STEP!----------------------
        # generate fist report.csv

        #if os.path.isfile(directory+'/GRPMX_report.csv'):
        #    df_report_complete = pd.read_csv(directory+'/GRPMX_report.csv', index_col=0)#).set_index('Unnamed: 0')
        #    df_report_complete = pd.concat([df_report_complete, df_report], axis=1)
        #    df_report_complete.to_csv(directory+'/GRPMX_report.csv')
        #else:
        #    df_report.to_csv(directory+'/GRPMX_report.csv') # solo la prima volta

        # FAST ALT------------
        df_report_complete = pd.concat([df_report_complete, df_report], axis=1)

        time_omega = datetime.now()
        full_runtime = time_omega - time_alpha
        print(gene+'_runtime:', full_runtime, ' Genes processed:', genes.index(gene), 'on ', len(genes))
        total_seconds = full_runtime.total_seconds()

        # save checkpoint----------------------
        if genes.index(gene) > 1 and genes.index(gene) % checkpoint == 0:
            complete_df = complete_df.reindex(columns=['gene','rsid', 'pmids', 'mesh'])
            complete_df.to_csv(directory+'/grpmx_filtered_output.csv')
            df_report_complete.to_csv(directory+'/GRPMX_report.csv')
            print("saved checkpoint")
        else:
            pass

    else:
        print(gene+' not present in DataBase')
        pass

#directory = r'H:\Il mio Drive\GRPM db\GRPM proteincodinggenes\grpm_random_04'
# Save complete csv
complete_df = complete_df.reindex(columns=['gene','rsid', 'pmids', 'mesh'])
complete_df.to_csv(directory+'/grpmx_filtered_output.csv')

df_report_complete.to_csv(directory+'/GRPMX_report.csv')

# #Update gene values (remove previous gene entry)
GRPMX_report = pd.read_csv(directory+'/GRPMX_report.csv', index_col=0)
time_load_1 = datetime.now()
for gene in grpm_genes_list:
    if gene+'.1' in GRPMX_report.columns:
        GRPMX_report = GRPMX_report.drop(columns = gene)
        GRPMX_report = GRPMX_report.rename(columns={gene+'.1': gene})
        print(genes.index(gene))
    else:
        pass
time_load_2 = datetime.now()
print(time_load_2-time_load_1)
GRPMX_report.to_csv(directory+'/GRPMX_report.csv')

time_finish = datetime.now()
time_batch = time_finish - time_start
print('time batch:',time_batch)
print('runtime/gene:', time_batch/len(genes))


# In[ ]:


# save checkpoint
complete_df = complete_df.reindex(columns=['gene','rsid', 'pmids', 'mesh'])
complete_df.to_csv(directory+'/grpmx_filtered_output.csv')

df_report_complete.to_csv(directory+'/GRPMX_report.csv')

restart_from = len(df_report_complete.T)
print('partial job, genes in survey '+topic_tag+' report:', restart_from)


# In[ ]:


df_report_complete.T


# In[ ]:


df_read = pd.read_csv(directory+'/grpmx_filtered_output.csv', index_col=0)#.gene.drop_duplicates()
#df_read.to_clipboard()
print('genes matching:', df_read.gene.nunique())
print('mesh matching:', df_read.mesh.nunique())
print('apply threshold in Analizer Module')
df_read


# # Check results

# In[ ]:


# Visualize GRPMX_report.csv
GRPMX_report = pd.read_csv(directory+'/GRPMX_report.csv', index_col=0).transpose().reset_index().rename(columns={'index':'gene'})
GRPMX_report.gene.drop_duplicates().to_clipboard()
print(len(GRPMX_report.gene.drop_duplicates()))
GRPMX_report


# In[ ]:


GRPMX_report[['reference_mesh', 'starting_pmidmesh', 'starting_pmid','starting_mesh','starting_rsid', 'matching_pmidmesh', 'matching_pmids', 'matching_mesh','matching_rsid', 'dropped_rsid', 'matching_rsid_pmid10','matching_rsid_pmid100']] = GRPMX_report[['reference_mesh', 'starting_pmidmesh', 'starting_pmid','starting_mesh','starting_rsid', 'matching_pmidmesh', 'matching_pmids', 'matching_mesh','matching_rsid', 'dropped_rsid','matching_rsid_pmid10','matching_rsid_pmid100']].astype(int)

GRPMX_report[['matching_mesh_ratio', 'matching_pmids_ratio','matching_pmidmesh_ratio', 'matching_rsid_ratio']] = GRPMX_report[['matching_mesh_ratio', 'matching_pmids_ratio','matching_pmidmesh_ratio','matching_rsid_ratio']].astype(float)

GRPMX_report_less = GRPMX_report[['matching_pmids','matching_pmids_ratio','matching_mesh','matching_rsid','matching_top10mesh','matching_top10rsid']]

GRPMX_report_sort = GRPMX_report.sort_values(by= 'matching_pmids',ascending=False)
GRPMX_report[['gene', 'matching_pmidmesh', 'matching_pmids',
'matching_mesh', 'matching_rsid', 'dropped_rsid', 'matching_mesh_ratio',
'matching_pmids_ratio', 'matching_pmidmesh_ratio',
'matching_rsid_ratio', 'matching_rsid_pmid10', 'matching_rsid_pmid100',
'matching_top10mesh', 'matching_top10rsid']].T


# In[ ]:


# Matching PMIDs in Database
GRPMX_report_sort = GRPMX_report.sort_values(by= 'matching_pmids',ascending=False)

x = GRPMX_report_sort.gene.iloc[:40]
y = GRPMX_report_sort['matching_pmids'].iloc[:40]
plt.figure(figsize=(5, 8))
plt.title('Matching PMIDs in Database', loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel('matching pmid', position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
#plt.clf()


# In[ ]:


# Add "interest index" to report:----------------------------------------------------------
max_match_pmids = int(GRPMX_report['matching_pmids'].max())
GRPMX_report_int = GRPMX_report
GRPMX_report_int['matching_pmids_index'] = round((GRPMX_report_int['matching_pmids']/max_match_pmids),3)

GRPMX_report_int['interest_index'] = round(GRPMX_report_int['matching_pmids_index'] * GRPMX_report_int['matching_pmids_ratio'],3)

GRPMX_report_int.set_index('gene').T.to_csv(directory+'/GRPMX_report_int.csv')


# In[ ]:


plt.clf()# Matching PMIDs in Database
GRPMX_report_sort = GRPMX_report.sort_values(by= 'matching_pmids_index',ascending=False)

x = GRPMX_report_sort.gene.iloc[:100]
y = GRPMX_report_sort['matching_pmids_index'].iloc[:100]
plt.figure(figsize=(5, 20))
plt.title('Matching PMIDs in Database', loc='center',pad=10)

plt.barh(x,y)
plt.gca().invert_yaxis()
plt.tick_params(axis='x', which='both', top=True, bottom=False, labeltop=True, labelbottom=False)
#plt.xlabel('pmid count', position=(0.5, 1.08))
plt.ylabel('genes')
plt.xlabel('matching pmid', position=(0.5, 1.08))
ax = plt.gca()
ax.xaxis.set_label_position('top')
#plt.savefig(r'C:\Users\Public\Database_SSD\GRPX_heavyload (RunOnly)\logs\Database_Gene('+str(len(GRPMX_report.gene))+')-PMID_filtered.png',dpi=300, bbox_inches = "tight")
#plt.clf()

