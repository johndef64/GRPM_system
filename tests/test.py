#%% md
# Import modules
import importlib
import subprocess
import torch
print("Torch version:",torch.__version__)
print("Is CUDA enabled?",torch.cuda.is_available())
if torch.cuda.is_available():
    print(torch.randn(1).cuda())

try:
    importlib.import_module('pygrpm')
except ImportError:
    subprocess.check_call(["pip", "install", "git+https://github.com/johndef64/GRPM_system.git"])

from pygrpm import *
#%% md
# GET Datasets
get_and_extract('grpm_dataset', record_id='14052302')
get_and_extract('nutrigenetic_dataset', record_id='14052302')
#%% md
# LOAD Datasets
pcg_grpm, rna_grpm, pseudo_grpm = grpm_importer()
grpm_nutrigen, grpm_nutrigen_int, grpm_nutrigen_int_gwas = nutrig_importer()

display(grpm_nutrigen_int)
#%% md
# SHOW Stats

pcg_grpm_stats = get_stats(pcg_grpm, group_by = 'gene')
display(pcg_grpm_stats)
#%%
grpm_nutrigen_stats = get_stats(grpm_nutrigen, group_by = 'gene')
display(grpm_nutrigen_stats)
#%% md
# QUERY GRPM Dataset

## MeSH Query Example

## GET MeSH Dataset ##
get_and_extract('ref-mesh', record_id='14052302')
get_topic_terms()

# LOAD MeSH
grpm_mesh = mesh_importer()
grpm_mesh.head()
#%%
# Random Query Example
mesh_query =  grpm_mesh['Preferred Label'].sample(10).to_list()

# Filter and get unique results
result = query_dataset(pcg_grpm, mesh_query, 'mesh')
display(result)
#%% md
## Build MeSH Query

# LOAD Language Model
sentence_transformer = load_language_model('dmis-lab/biobert-v1.1')

# Get MeSH embeddings
series2 = grpm_mesh['Preferred Label'].reset_index(drop=True)
mesh_embeddings = extract_embedding(series2.to_list(), sentence_transformer)
#%%
# User defined Topic Terms
topic_terms_sample = ["diet ketogenic",
                      "diet reducing",
                      "diet sodium-restricted",
                      "diet",
                      "dietary",
                      "dietetics",
                      "dyslipidemias",
                      "eating disorders",
                      "feeding and eating disorde",
                      "food hypersensitivity",
                      "foodborne diseases",
                      "gastrointestinal diseases",
                      "hypercholesterolemia",
                      "hyperglycemia",
                      "hyperlipidemias",
                      "hyperphagia",
                      "hypoglycemia",
                      "hypophagia",
                      "insulin resistance",
                      ]
series1 = pd.Series(topic_terms_sample)
#%%
# Extract MeSH Query
tab = create_corr_table(series1, series2, sentence_transformer, mesh_embeddings)
#%%
mesh_query = tab[tab.similarity >= 0.90].list2.to_list()
print('\n\nMeSH Query:', mesh_query)
#%% md
## Execute MeSH Query
#%%
# Filter and get unique results
result = query_dataset(pcg_grpm, mesh_query, 'mesh')
display(result)
#%% md
# QUERY Nutrigenetic Dataset
#%%
# Gene Query on Nutrigenetic ds
topic =  grpm_nutrigen_int.topic[0]

print(f'Displaying "{topic}" topic')
# Filter and get unique results
result = query_dataset(grpm_nutrigen_int, [topic], 'topic')
display(result)
#%% md
## Gene Query Example
#%%
# Gene Query on Nutrigenetic ds
my_genes = (
    'FTO',
    'APOB',
    'G6PD'
)
# Filter and get unique results
result = query_dataset(grpm_nutrigen_int, my_genes, 'gene')
display(result)
#%%
# Gene Query on Nutrigenetic-GWAS ds
result = query_dataset(grpm_nutrigen_int_gwas, my_genes, 'GRPM_GENE')
display(result)
#%%
