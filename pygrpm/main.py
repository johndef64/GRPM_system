import os
import io
import sys
import glob
import zipfile
import requests
import gdown
import pyarrow
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime
import pyperclip as pc

def simple_bool(message):
    choose = input(message+" (y/n): ").lower()
    your_bool = choose in ["y", "yes","yea","sure"]
    return your_bool

def get_gitfile(url, flag='', dir = os.getcwd()):
    url = url.replace('blob','raw')
    response = requests.get(url)
    file_name = flag + url.rsplit('/',1)[1]
    file_path = os.path.join(dir, file_name)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"File downloaded successfully. Saved as {file_name}")
    else:
        print("Unable to download the file.")


def get_from_zenodo(file, dir=os.getcwd(), record_id='8205724'):
    print('Requesting file to Zenodo...')
    # Construct the download URL
    url = f'https://zenodo.org/record/{record_id}/files/{file}?download=1'
    extracted_folder_name = dir
    # Define the output path for the file
    output_path = os.path.join(extracted_folder_name, file)
    # Use gdown to download the file
    gdown.download(url, output_path, quiet=False)
    return output_path

def get_and_extract(file, dir=os.getcwd(), record_id='8205724', ext='.zip', remove_zip=True):
    zip_file_name = file + ext
    extracted_folder_name = dir
    # Download the ZIP file using get_from_zenodo
    zip_file_path = get_from_zenodo(file + ext, dir, record_id)
    # Extract the ZIP contents
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        print('Extracting...')
        zip_ref.extractall(extracted_folder_name)
    print(f"ZIP file '{zip_file_name}' extracted to '{extracted_folder_name}' successfully.")
    # Remove the ZIP file after extraction
    if remove_zip: os.remove(zip_file_path)


##### Specific Functions ####
def grpm_importer():
    t1 = datetime.now()
    grpm_dataset = pd.read_parquet('grpm_dataset//grpm_dataset.parquet', engine='pyarrow')
    print('Importing time: ',datetime.now() - t1)
    def subset_df(df, subset):
        return df[df['type'] == subset].reset_index(drop=True)
    pcg_grpm    = subset_df(grpm_dataset, subset='PCG')
    rna_grpm    = subset_df(grpm_dataset, subset='RNA')
    pseudo_grpm = subset_df(grpm_dataset, subset='PSD')
    #pcg_grpm = pd.read_parquet('grpm_dataset/grpm_db_pcg/grpm_dataset.parquet', engine = 'pyarrow')
    def printer(df, tag):
        print(tag,   round(df.memory_usage().sum()    / 1024 / 1024, 2), 'MB')
    printer(pcg_grpm, 'pcg:')
    printer(rna_grpm, 'rna:')
    printer(pseudo_grpm, 'pseudo:')
    return pcg_grpm, rna_grpm, pseudo_grpm

def nutrig_importer():
    grpm_nutrigen = pd.read_parquet('nutrigenetic_dataset/grpm_nutrigen.parquet')
    grpm_nutrigen_int = pd.read_parquet('nutrigenetic_dataset/grpm_nutrigen_int.parquet')
    gwas_path = 'nutrigenetic_dataset/grpm_nutrigen_int_gwas.parquet'
    is_gwas = os.path.exists(gwas_path)
    grpm_nutrigen_int_gwas = pd.read_parquet(gwas_path) if is_gwas else pd.DataFrame()
    def printer(df, tag):
        print(tag,   round(df.memory_usage().sum()    / 1024 / 1024, 2), 'MB')
    printer(grpm_nutrigen, 'nutrigen dataset:')
    printer(grpm_nutrigen_int, 'nutrigen dataset filtered:')
    printer(grpm_nutrigen_int_gwas, 'nutrigen gwas dataset:')

    return   grpm_nutrigen, grpm_nutrigen_int, grpm_nutrigen_int_gwas

#%%
def get_stats(ds, group_by = 'gene', sublevel = 'unique'):
    print('Computing Stats...')
    t1 = datetime.now()
    ds = ds.groupby(group_by).describe(include=['object'])
    def norm_describe(df, sublevel):
        return df.loc[:, df.columns.get_level_values(1) == sublevel]
    print('runtime: ', datetime.now() - t1)
    return norm_describe(ds, sublevel)

#%%
def query_dataset(ds, my_tuple, field = 'mesh'):
    ds = ds[ds[field].isin(my_tuple)].drop_duplicates()
    return ds
