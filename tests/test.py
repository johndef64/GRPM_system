from pygrpm import *

#%%
### GET Datasets ###
get_and_extract('grpm_dataset', record_id='8205724')
get_and_extract('nutrigenetic_dataset', record_id='8205724')

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
### QUERY Datasets ###

# MeSH Query on GRPM ds
my_mesh = (
    'Sodium-Coupled Vitamin C Transporters',
    'Micronutrients',
    'Biotinidase Deficiency',
    'Lipid Metabolism Disorders',
    'Vitamin D-Binding Protein',
    'Pyridoxal Phosphate',
    'Choline-Phosphate Cytidyltransferase',
    'Vitamin D3 24-Hydroxylase'
)

# Filter and get unique results
result = query_dataset(pcg_grpm, my_mesh, 'mesh')
print(result)

#%%

# Gene Query on Nutrigenetic ds
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