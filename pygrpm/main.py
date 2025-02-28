import os
import io
import sys
import glob
import torch
import gdown
import pickle
import pyarrow
import zipfile
import requests
import importlib
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime
from sentence_transformers import SentenceTransformer, util

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

def import_mesh_embeddings(file_path = 'ref-mesh/GrpmMeshEmbeddings_biobert-v1.1.pkl'):
    url ="https://github.com/johndef64/GRPM_system/raw/refs/heads/main/ref-mesh/GrpmMeshEmbeddings_biobert-v1.1.pkl"
    if not os.path.exists(file_path):
        print('Downloading...')
        gdown.download(url, file_path, quiet=False)

    # Import the mesh_embeddings back from the pickle file
    with open(file_path, "rb") as file:
        grpm_mesh_embeddings = pickle.load(file)

    print("Done")
    return grpm_mesh_embeddings

####


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


def get_latest_mesh_data():
    gdown.download("https://data.bioontology.org/ontologies/MESH/download?apikey=8b5b7825-538d-40e0-9e9e-5ab9274a9aeb&download_format=csv", output="ref-mesh/MESH.gz")

def import_grpm_mesh(file_path = "ref-mesh/GrpmMesh.parquet"):
    os.makedirs('ref-mesh', exist_ok=True)
    if not os.path.isfile(file_path):
        gdown.download("https://github.com/johndef64/GRPM_system/raw/refs/heads/main/ref-mesh/GrpmMesh.parquet", output=file_path)
    grpm_mesh_data = pd.read_parquet(file_path)
    print('GRPM MeSH count:', grpm_mesh_data['Preferred Label'].nunique())
    return grpm_mesh_data

def import_grpm_mesh_zenodo(file_path='ref-mesh/MESH_STY_LITVAR1.csv'):
    if not os.path.exists(file_path):
        get_and_extract('ref-mesh', record_id='14052302')
    grpm_mesh = pd.read_csv(file_path, index_col=0)
    print('GRPM MeSH count:', grpm_mesh['Preferred Label'].nunique())
    print('semantic types:', grpm_mesh['Semantic Types Label'].nunique())
    grpm_mesh = grpm_mesh[['Preferred Label', 'Semantic Types Label', 'Class ID', 'mesh_id','Semantic Types']]
    return grpm_mesh


#%%
def get_stats(dataset, group_by = 'gene', sublevel = 'unique', gi_sort=False):
    """
    :param dataset: dataset to get stats for
    :param group_by: column to group by
    :param sublevel: sublevel to display
    :param gi_sort: sort results by GI
    :return: stats dataframe
    """
    print('Computing Stats...')
    t1 = datetime.now()
    group = dataset.groupby(group_by).describe(include=['object'])
    def norm_describe(df, sublevel):
        return df.loc[:, df.columns.get_level_values(1) == sublevel]

    print('runtime: ', datetime.now() - t1)

    stats = norm_describe(group, sublevel)

    if gi_sort and group_by in ['gene', 'GRPM_GENE']:
        # Create a sorted index by finding the order of elements in 'gene_sorted_gi'
        gene_sorted_gi = dataset[group_by].drop_duplicates().to_list()
        sorted_index = stats.index.map(lambda x: gene_sorted_gi.index(x) if x in gene_sorted_gi else float('inf'))

        # Sort 'stats' based on the created sorted index
        stats = stats.iloc[sorted_index.argsort()]

    return stats

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

def test_cuda():
    print("Torch version:",torch.__version__)
    print("Is CUDA enabled?",torch.cuda.is_available())
    # if torch.cuda.is_available():
    #     print(torch.randn(1).cuda())
    # else:
    #     print("""if CUDA not available and Windows system: uninstall torch and install "pip3 install torch --index-url https://download.pytorch.org/whl/cu118""")

def load_language_model(language_model = 'dmis-lab/biobert-v1.1'):
    # Choose embedding model
    return SentenceTransformer(language_model)

def extract_embedding(input_text, model):
    # Encode the input text to get the embedding
    embedding = model.encode(input_text, show_progress_bar=True)
    return embedding

def get_mesh_embeddings(grpm_mesh, model, file_path, retrain =False, definitions=False, synonyms=False):
    """

    :param grpm_mesh: GRPM Mesh DataFrame
    :param model: SentenceTransformer(language_model)
    :param file_path: mesh_embeddings.pkl file path
    :return: mesh embeddings dict
    """

    if torch.cuda.is_available() and not os.path.isfile(file_path) or retrain:
        print("Is CUDA enabled?",torch.cuda.is_available())
        print("No embeddings. Generating embeddings...")
        # Generate MeSH Embeddings
        grpm_meshes = list(grpm_mesh['Preferred Label'])
        if definitions:
            grpm_mesh['Definitions'] = grpm_mesh['Definitions'].fillna(" ")
            grpm_mesh['combined'] = grpm_mesh['Preferred Label'] + ": " + grpm_mesh['Definitions']
            grpm_meshes_to_encode = list(grpm_mesh['combined'])
        elif synonyms:
            grpm_mesh['Synonyms'] = grpm_mesh['Synonyms'].fillna(" ")
            grpm_mesh['combined'] = grpm_mesh['Preferred Label'] + "; " + grpm_mesh['Synonyms']
            grpm_meshes_to_encode = list(grpm_mesh['combined'])
        else:
            grpm_meshes_to_encode = grpm_meshes

        mesh_embeddings = extract_embedding(grpm_meshes_to_encode, model)
        grpm_mesh_embeddings = {"meshes":grpm_meshes, "embeddings":mesh_embeddings}

        # Open the file in write-binary mode to store the pickle
        with open(file_path, 'wb') as file:
            print("Saving embeddings in {}".format(file_path))
            # Use pickle to dump the dictionary into the file
            pickle.dump(grpm_mesh_embeddings, file)
    else:
        # Download and import MeSH Embeddings
        print("Importing pretrained embeddings...")
        grpm_mesh_embeddings = import_mesh_embeddings(file_path)


    return grpm_mesh_embeddings



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

def create_corr_table(series1, series2, model, series2_embeddings = [], n=3):
    """
    Returns only the top correspondence for each "series1" element

    :param series1:
    :param series2:
    :param model:
    :param series2_embeddings:
    :param n: number of top matches
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

    #####
    # List to store rows for the correspondence table
    data_rows = []

    # Loop through each row to find top N correspondences
    for row_index in similarities_df.index:
        # Get the top N similarities for each row
        top_n_similarities = similarities_df.loc[row_index].nlargest(n)
        for col_index, sim_value in top_n_similarities.items():
            # Create a new row entry
            new_row = {
                'list1': series1[row_index],
                'list2': series2[col_index],
                'similarity': sim_value
            }
            data_rows.append(new_row)  # Append new row to list

    # Create DataFrame from accumulated rows
    correspondence_df = pd.DataFrame(data_rows)

    # Sort by 'similarity' column
    return correspondence_df.sort_values(by=['similarity'], ascending=False)
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
    query_embedding = model.encode(user_query, convert_to_tensor=True, show_progress_bar=False)
    if mesh_embeddings is not None:
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


def get_mesh_query(user_query,
                   mesh_terms,
                   model,
                   mesh_embeddings=None,
                   threshold=0.88):
    df = get_mesh_rankings(user_query,  mesh_terms, model, mesh_embeddings=mesh_embeddings,)
    pertinent_terms = filter_mesh_scores(df, threshold=threshold)
    return pertinent_terms

