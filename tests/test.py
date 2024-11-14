import pandas as pd

from pygrpm import *

#%%
### GET Datasets ###
get_and_extract('grpm_dataset', record_id='14052302')
#%%
get_and_extract('nutrigenetic_dataset', record_id='14052302')

#%%
### LOAD Datasets ###
pcg_grpm, rna_grpm, pseudo_grpm = grpm_importer()
grpm_nutrigen, grpm_nutrigen_int, grpm_nutrigen_int_gwas = nutrig_importer()


#%%
### GET Stats ###
result = get_stats(pcg_grpm, group_by = 'gene')
print(result)
#%%
result = get_stats(grpm_nutrigen, group_by = 'gene')
print(result)

#%%
## Build MeSH Query ##
get_and_extract('ref-mesh', record_id='14052302')
get_topic_terms()
#%%
# LOAD MeSH
grpm_mesh = mesh_importer()

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
                      "kidney diseases",
                      ]
series1 = pd.Series(topic_terms_sample)

#%%
# Define MeSH Query
tab = create_corr_table(series1, series2, sentence_transformer, mesh_embeddings)
mesh_query = tab[tab.similarity >= 0.90].list2.to_list()
print(mesh_query)
#%%
### QUERY Datasets ###

# Filter and get unique results
result = query_dataset(pcg_grpm, mesh_query, 'mesh')
print(result)

#%%

# Gene Query on Nutrigenetic dataset
my_genes = (
    'FTO',
    'APOB',
    'G6PD'
)
# Filter and get unique results
result = query_dataset(grpm_nutrigen_int, my_genes, 'gene')
print(result)
#%%

# Gene Query on Nutrigenetic-GWAS ds
result = query_dataset(grpm_nutrigen_int_gwas, my_genes, 'GRPM_GENE')
print(result)

#%%

