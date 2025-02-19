import os
import io
import sys
import glob
import zipfile
import requests
import importlib
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

def check_and_install_module(module_name):
    try:
        # Check if the module is already installed
        importlib.import_module(module_name)
        print(f"The module '{module_name}' is already installed.")
    except ImportError:
        # If the module is not installed, try installing it
        x = simple_bool(
            "\n" + module_name + "  module is not installed.\nwould you like to install it?")
        if x:
            import subprocess
            subprocess.check_call(["pip", "install", module_name])
            print(f"The module '{module_name}' was installed correctly.")
        else:
            pass

def get_gitfile(url: str, flag='', dir = os.getcwd(), print_=True):
    url = url.replace('blob','raw')
    response = requests.get(url)
    file_name = flag + url.rsplit('/',1)[1]
    file_path = os.path.join(dir, file_name)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
        if print_: print(f"File downloaded successfully. Saved as {file_name}")
    else:
        print("Unable to download the file.")

def get_file(url, file_name, dir = os.getcwd()):
    url = url
    file_name = file_name
    response = requests.get(url)
    if response.status_code == 200:
        content = response.content
        file_path = os.path.join(dir, file_name)
        with open(file_path, 'wb') as file:
            file.write(content)


releases = {'0.1':'8205724',
            '0.2':'14052302'}

def get_from_zenodo(file: str, record_id=releases['0.2'], dir=os.getcwd()): # 14052302
    # Construct the download URL
    url = f'https://zenodo.org/record/{record_id}/files/{file}?download=1'
    print(f'Downloading {file} from {url}')
    extracted_folder_name = dir
    # Define the output path for the file
    output_path = os.path.join(extracted_folder_name, file)
    # Use gdown to download the file
    gdown.download(url, output_path, quiet=False)
    return output_path

def get_and_extract(file: str, record_id=releases['0.2'], dir=os.getcwd(), ext='.zip', remove_zip=True):
    zip_file_name = file + ext
    extracted_folder_name = dir
    # Download the ZIP file using get_from_zenodo
    zip_file_path = get_from_zenodo(file + ext, record_id, dir)
    # Extract the ZIP contents
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        print('Extracting...')
        zip_ref.extractall(extracted_folder_name)
    print(f"ZIP file '{zip_file_name}' extracted in '{extracted_folder_name}' successfully.")
    # Remove the ZIP file after extraction
    if remove_zip: os.remove(zip_file_path)

def import_mesh_embeddings():
    file_path = 'ref-mesh/GrpmMeshEmbeddings_biobert-v1.1.pkl'
    url ="https://github.com/johndef64/GRPM_system/raw/refs/heads/main/ref-mesh/GrpmMeshEmbeddings_biobert-v1.1.pkl"
    if not os.path.exists(file_path):
        gdown.download(url, file_path, quiet=False)

    # Import the mesh_embeddings back from the pickle file
    with open(file_path, "rb") as file:
        grpm_mesh_embeddings = pickle.load(file)
    return grpm_mesh_embeddings

def get_topic_terms():
    files = ["topic_terms_nutri_old.csv",
             "topic_terms_ng_xeno.csv",
             "topic_terms_ng_vitam.csv",
             "topic_terms_ng_oxi_stress.csv",
             "topic_terms_ng_ob_bmi.csv",
             "topic_terms_ng_nutri.csv",
             "topic_terms_ng_intol.csv",
             "topic_terms_ng_eat_taste.csv",
             "topic_terms_ng_dmt2_ms.csv",
             "topic_terms_ng_cvd.csv",
             "topic_terms_ng_aller.csv",]
    for file in files :
        url = f"https://github.com/johndef64/GRPM_system/blob/main/ref-mesh/topic_terms/{file}"
        get_gitfile(url, flag='', dir = 'ref-mesh/topic_terms', print_=False)


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

def mesh_importer():
    grpm_mesh = pd.read_csv('ref-mesh/MESH_STY_LITVAR1.csv', index_col=0)
    print('GRPM MeSH count:', grpm_mesh['Preferred Label'].nunique())
    print('semantic types:', grpm_mesh['Semantic Types Label'].nunique())
    grpm_mesh = grpm_mesh[['Preferred Label', 'Semantic Types Label', 'Class ID', 'mesh_id','Semantic Types']]
    return grpm_mesh

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
    """
    :param ds: dataset
    :param my_tuple: query by element list
    :param field: query field
    :return: filtered dataset
    """
    ds = ds[ds[field].isin(my_tuple)].drop_duplicates()
    return ds
#%%

#%%
## SEMANTIC SIMILARITY

from sentence_transformers import SentenceTransformer, util

def load_language_model(language_model = 'dmis-lab/biobert-v1.1'):
    # Choose embedding model
    return SentenceTransformer(language_model)

def extract_embedding(input_text, model):
    # Encode the input text to get the embedding
    embedding = model.encode(input_text, show_progress_bar=True)
    return embedding

# Function to compute cosine similarity
def cosine_similarity(vec1, vec2):
    dot_product = np.dot(vec1, vec2)
    norm_vec1 = np.linalg.norm(vec1)
    norm_vec2 = np.linalg.norm(vec2)
    return dot_product / (norm_vec1 * norm_vec2)

def cosine_distance(array1, array2):
    # Calculate the cosine similarities
    similarities = np.zeros((array1.shape[0], array2.shape[0]))
    for i in tqdm(range(array1.shape[0])):
        for j in range(array2.shape[0]):
            similarities[i, j] = cosine_similarity(array1[i], array2[j])
    return similarities


def calculate_cos_sim(array1, array2):
    # Calculate the cosine similarities
    similarities = np.zeros((array1.shape[0], array2.shape[0]))
    for i in tqdm(range(array1.shape[0])):
        for j in range(array2.shape[0]):
            similarities[i, j] = cosine_similarity(array1[i], array2[j])
    return similarities

def calculate_cos_sim_tensor(tesor1, tensor2):
    # Calculate cosine similarities
    cosine_scores = util.pytorch_cos_sim(tesor1, tensor2)
    return cosine_scores

### METHOD 1 ###
def create_corr_table(series1, series2, model, series2_embeddings = []):
    """

    :param series1:
    :param series2:
    :param model:
    :param series2_embeddings:
    :return:
    """
    array1 = extract_embedding(series1.to_list(), model)
    if len(series2_embeddings) > 1:
        array2 = series2_embeddings
    else:
        array2 = extract_embedding(series2.to_list(), model)


    # Calculate the cosine similarities
    similarities = calculate_cos_sim(array1, array2)

    # Create a DataFrame with the similarities
    similarities_df = pd.DataFrame(similarities, index=[i for i in range(array1.shape[0])],
                                   columns=[j for j in range(array2.shape[0])])

    # Find the index of the maximum value in each row
    max_indices = similarities_df.idxmax(axis=1)
    # Find the maximum value in each row
    max_values = similarities_df.max(axis=1)

    # Create a DataFrame with the correspondence between row index and column values
    correspondence_table = pd.DataFrame({
        "list1_id": max_indices.index,
        "list2_id": max_indices.values,
        "Max Value": max_values.values
    })

    data_rows = []  # List to store each row as a dictionary

    for i in range(len(correspondence_table)):  # Loop over index range
        list1_id = correspondence_table['list1_id'][i]
        list2_id = correspondence_table['list2_id'][i]
        sim_value = correspondence_table['Max Value'][i]

        # Create a dictionary for the new row
        new_row = {
            'list1': series1[list1_id],
            'list2': series2[list2_id],
            'similarity': sim_value
        }

        data_rows.append(new_row)  # Append to list of rows

    df = pd.DataFrame(data_rows)
    df.sort_values(by=['similarity'], ascending=False)
    return df
#%%

### METHOD 2 ###
import torch

def get_mesh_rankings(user_query,  #"Enter your query here"
                      mesh_terms,
                      model,
                      mesh_embeddings=None,
                      ):
    """
    :param user_query: textual query e.c. "Diabetic Disorders, insulin resistance"
    :param mesh_terms: fixed mesh term list
    :param model: llm model
    :param mesh_embeddings: pre-made mesh embeddings (as numpy.array)
    :return: mesh ranking dataframe
    """
    # User query

    # Get embeddings for the query and MESH terms
    query_embedding = model.encode(user_query, convert_to_tensor=True, show_progress_bar=True)
    if mesh_embeddings:
        mesh_embeddings = torch.from_numpy(mesh_embeddings)
    else:
        mesh_embeddings = model.encode(mesh_terms, convert_to_tensor=True, show_progress_bar=True)

    # Get the default device (GPU if available, else CPU)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Move tensors to the same device
    query_embedding = query_embedding.to(device)
    mesh_embeddings = mesh_embeddings.to(device)

    # Calculate cosine similarities
    cosine_scores = util.pytorch_cos_sim(query_embedding, mesh_embeddings)

    df = pd.DataFrame({
        'mesh_terms': mesh_terms,
        'cosine_scores': cosine_scores.tolist()[0]
    })
    return df


# def filter_mesh_scores(mesh_terms, cosine_scores, threshold=0.88):
#     # Select relevant terms based on similarity threshold
#     cosine_threshold = threshold # Adjust threshold according to tests
#     pertinent_terms = [mesh_terms[i] for i in range(len(mesh_terms)) if cosine_scores[0][i] > cosine_threshold]
#
#     print("Relevant terms:", pertinent_terms)
#     return pertinent_terms


def filter_mesh_scores(df, threshold=0.88 ):
    """
    :param df: output of get_mesh_rankings()
    :param threshold: Adjust threshold according to tests
    :return: mesh list
    """

    # Select relevant terms based on similarity threshold
    pertinent_terms = df[df.cosine_scores > threshold].mesh_terms.to_list()

    print("Related MeSH:", pertinent_terms)
    return pertinent_terms



