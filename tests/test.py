#!/usr/bin/env python
# coding: utf-8

# # Import modules

# In[1]:

import importlib
import subprocess

try:
    importlib.import_module('pygrpm')
except ImportError:
    subprocess.check_call(["pip", "install", "git+https://github.com/johndef64/GRPM_system.git"])

from pygrpm import *


# # GET Datasets

# In[3]:


### GET Datasets ###
if not os.path.exists('grpm_dataset/grpm_dataset.parquet'):
    get_and_extract('grpm_dataset', record_id='14052302')
    get_and_extract('nutrigenetic_dataset', record_id='14052302')


# # LOAD Datasets

# In[2]:


### LOAD Datasets ###
"""
1. GRPM Datset
2. Nutrigenetic Dataset
3. Nutrigenetic Dataset + GWAS
"""

pcg_grpm, rna_grpm, pseudo_grpm = grpm_importer()
grpm_nutrigen, grpm_nutrigen_int, grpm_nutrigen_int_gwas = nutrig_importer()


# # SHOW Stats

# In[4]:

# Full GRPM Dataset Build
pcg_grpm_stats = get_stats(pcg_grpm, group_by = 'gene')
display(grpm_nutrigen_int.head(10))
display(pcg_grpm_stats.head(20))


# In[7]:

# Nutrigenetic Dataset (10 Topics)
grpm_nutrigen_stats = get_stats(grpm_nutrigen, group_by = 'gene', gi_sort=True)
display(grpm_nutrigen.head(10))
display(grpm_nutrigen_stats.head(20))

# In[5]:

# Nutrigenetic Dataset + GWAS
grpm_nutrigen_int_gwas_stats = get_stats(grpm_nutrigen_int_gwas, group_by='GRPM_GENE', gi_sort=True)
display(grpm_nutrigen_int_gwas.head(10))
display(grpm_nutrigen_int_gwas_stats.head(20))


# # QUERY GRPM Dataset

# ## MeSH Query Example

# In[3]:


# LOAD MeSH

grpm_mesh= import_grpm_mesh()
grpm_mesh.head()


# In[5]:


# Random Query Example
mesh_query =  grpm_mesh['Preferred Label'].drop_duplicates().sample(10).to_list()

# Filter and get unique results
result = query_dataset(pcg_grpm, mesh_query, 'mesh')
display(result)


# ## Build MeSH Query
# [CUDA recommended] - In Colab:: load Runtime with GPU

# In[5]:


# LOAD Language Model
MODEL = 'dmis-lab/biobert-v1.1'
model = load_language_model(MODEL)
file_path = 'ref-mesh/GrpmMeshSynEmbeddings_biobert-v1.1.pkl'

# Get MeSH embeddings
test_cuda()
grpm_mesh_embeddings = get_mesh_embeddings(grpm_mesh, model, file_path)

grpm_meshes = grpm_mesh_embeddings['meshes']
mesh_embeddings = grpm_mesh_embeddings['embeddings']


# In[14]:


# User defined Topic Terms

user_query = "diet ketogenic, diet reducing, diet sodium-restricted, diet, dietary, dietetics, dyslipidemias, eating disorders, feeding and eating disorder, food hypersensitivity, foodborne diseases, gastrointestinal diseases, hypercholesterolemia, hyperglycemia, hyperlipidemias, hyperphagia, hypoglycemia, hypophagia, insulin resistance"  # comma separated list

topic_terms_list = user_query.split(',')
topic_terms = pd.Series(topic_terms_list)


# Extract MeSH Query
tab = create_corr_table(topic_terms, grpm_meshes, model, mesh_embeddings)

threshold = 0.84 # set similarity threshold
mesh_query = tab[tab.similarity >= threshold].list2.to_list()
print('\n\nMeSH Query:', mesh_query)


# MeSH Query: ['Diet', 'Diet', 'Dyslipidemias', 'Hypersensitivity', 'Gastrointestinal Diseases', 'Hypercholesterolemia', 'Hyperglycemia', 'Hyperlipidemias', 'Hyperphagia', 'Hypoglycemia']

# ## Execute MeSH Query

# In[ ]:


# Filter and get unique results
result = query_dataset(pcg_grpm, mesh_query, 'mesh')
display(result)


#

# # QUERY Nutrigenetic Dataset

# In[6]:


get_stats(grpm_nutrigen_int, 'topic')


# ## 1. Query by Nutritional Topic

# In[17]:


# Select Topic
topic = "Vitamin and Micronutrients Metabolism and Deficiency-Related Diseases"

# Filter and get unique results
topic_data = query_dataset(grpm_nutrigen_int, [topic], 'topic')
print(f'Displaying "{topic}" topic')
display(topic_data)


# In[ ]:


# Get Topic Data Stats
stats = get_stats(topic_data, "gene", gi_sort=True)
stats


# ### Select Topic on Nutrigenetic-GWAS dataset

# In[36]:


# Select Topic on Nutrigenetic-GWAS dataset
topic = "Vitamin and Micronutrients Metabolism and Deficiency-Related Diseases"

topic_data_gwas = query_dataset(grpm_nutrigen_int_gwas, [topic], 'GRPM_TOPIC')
display(topic_data_gwas.head())
# Get Topic Data Stats
stats = get_stats(topic_data_gwas, group_by = "GRPM_GENE")
stats


#

# ## 2. Advanced query (use case)
#
# Exploring the genetic determinants of nutritional status involves understanding how genetic variations influence the intake and utilization of micronutrients, impacting nutrient transport, metabolism, and cellular uptake.
#
# > *Micronutrients such as trace elements and vitamins are important as enzyme cofactors in the metabolism of all cells in the body and therefore key to determining nutritional status.*
#
#
# Build a composite query:
# - Nutritional Status, Mechanisms of Micronutrient Metabolism, and Micronutrient Measurement

# In[11]:


from pygrpm import *

# Download and import MeSH Embeddings
grpm_mesh_embeddings = import_mesh_embeddings()

grpm_meshes = grpm_mesh_embeddings['meshes']
mesh_embeddings = grpm_mesh_embeddings['embeddings']

# Define queries using natural language
QUERIES =[
    # 1. **Nutritional Status**:
    ["Measurement of nutritional status", 0.85],
    ["Assess essential micronutrients", 0.84],
    ["Focus on vitamins like vitamin A, D, and B-vitamins.", 0.84],
    ["trace minerals such as iron, zinc, and iodine.", 0.84],

    # 2"**Mechanisms of Micronutrient Metabolism**:
    ["cellular processes for micronutrient absorption.", 0.87],
    ["transport and transformation of nutrients.", 0.84],
    ["nutrients storage mechanisms.", 0.84],
    ["cofactors in nutrient utilization.", 0.84],
    ["homeostasis of nutrients.", 0.84],
]


mesh_query = []
for i in range(len(QUERIES)):
    query =  QUERIES[i]
    print("\n",query)
    meshes = get_mesh_query(query[0], grpm_meshes, model, mesh_embeddings=mesh_embeddings, threshold=query[1])
    mesh_query.extend(meshes)

mesh_query = list(set(mesh_query))


# Get Genes and Variants possibly related to the Query

# In[1]:


# Filter and get unique results
query_result = query_dataset(grpm_nutrigen_int, mesh_query, 'mesh')
query_result


# In[20]:


get_stats(query_result, "gene", gi_sort=True)


# In[26]:


## Select specific genes

# Gene Query
my_genes = "VDR, PNPLA3, PNPLA3"

# Filter and get unique results
gene_panel = query_dataset(query_result, my_genes.split(','), 'gene')
display(gene_panel)


#%%
