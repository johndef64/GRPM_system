import pandas as pd
from sentence_transformers import SentenceTransformer, util

# Load the BioBERT model via sentence-transformers
model = SentenceTransformer('dmis-lab/biobert-v1.1')

# List of MESH terms
grpm_nutrigen_int = pd.read_parquet('../nutrigenetic_dataset/grpm_nutrigen_int.parquet', engine='pyarrow')
grpm_meshs = grpm_nutrigen_int.mesh.drop_duplicates()
mesh_terms = grpm_meshs.to_list()  # sostituisci con termini MESH reali

# User query
user_query = """Diabetic Disorders, insulin resistance""" #"Enter your query here"


# Get embeddings for the query and MESH terms
query_embedding = model.encode(user_query, convert_to_tensor=True, show_progress_bar=True)
mesh_embeddings = model.encode(mesh_terms, convert_to_tensor=True, show_progress_bar=True)

# Calculate cosine similarities
cosine_scores = util.pytorch_cos_sim(query_embedding, mesh_embeddings)
print(cosine_scores)
#%%

# Select relevant terms based on similarity threshold
cosine_threshold = 0.85  # Adjust threshold according to tests
pertinent_terms = [mesh_terms[i] for i in range(len(mesh_terms)) if cosine_scores[0][i] > cosine_threshold]

print("Termini pertinenti:", pertinent_terms)
#%%
#%%