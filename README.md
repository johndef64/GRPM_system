# GRPM System

The GRPM (Gene-Rsid-Pmid-Mesh) system is a comprehensive tool designed to integrate and analyze genetic polymorphism data associated with nutrition. It comprises five modules that facilitate data retrieval, merging, analysis, and incorporation of GWAS data.

## Overview

- [Introduction](#introduction)
- [Modules](#modules)
- [Usage](#usage)
- [Requirements](#requirements)
- [Installation](#installation)

## Introduction

The GRPM System aims to build a comprehensive dataset of human genetic polymorphisms associated with nutrition. By combining data from multiple sources and utilizing MeSH terms as a framework, the system enables researchers and nutritionists to explore gene-diet interactions and personalized nutrition interventions.


## Modules

The GRPM System comprises five modules that perform various tasks. Module 01 retrieves data from LitVar and PubMed databases, merging them into a CSV format. Module 02 generates a coherent MeSH term list with the help of the ChatGPT language model and the OpenAI API. Module 03 integrates the MeSH term list into the database, extracting a survey for further analysis, while Module 04 assesses reports, GRPM association data, and uses `matplotlib`, and `seaborn` for data visualization. Lastly, Module 05 incorporates GWAS data from the complete catalog, associating GWAS phenotypes and potential risk/effect alleles with GRPM relationships.

The GRPM System consists of the following modules:

1. [Database Builder](https://github.com/johndef64/GRPM_playground/blob/main/GRPM_01_database_builder.ipynb): Retrieves and merges data from LitVar and PubMed databases.
-  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_playground/blob/main/GRPM_01_database_builder.ipynb)
2. [Reference Mesh List Builder](https://github.com/johndef64/GRPM_playground/blob/main/GRPM_02_ref-mesh_builder.ipynb): Generates a coherent MeSH term list for exploring the database.
-   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_playground/blob/main/GRPM_02_ref-mesh_builder.ipynb)
3. [Database and Survey Integration](https://github.com/johndef64/GRPM_playground/blob/main/GRPM_03_database_survey.ipynb): Incorporates the MeSH term list into the database and extracts a survey.
-   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_playground/blob/main/GRPM_03_database_survey.ipynb)
4. [Analysis of Reports and GRPM Association Data](https://github.com/johndef64/GRPM_playground/blob/main/GRPM_04_data-analyzer.ipynb): Analyzes reports and GRPM association data.
-  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_playground/blob/main/GRPM_04_data-analyzer.ipynb)
5. [Incorporation of GWAS Data](https://github.com/johndef64/GRPM_playground/blob/main/GRPM_05_gwas_data_analyzer.ipynb): Integrates GWAS data into the GRPM surveys.
-   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_playground/blob/main/GRPM_05_gwas_data_analyzer.ipynb)

## Usage

Detailed instructions on how to use each module of the GRPM System can be found inside the relative Jupyter Module provided in the repository. Make sure to follow the instructions and install the necessary Python modules specified for each module.

## Requirements

The GRPM System has the following requirements:

- `Python 3.9 or above`
- `pandas`
- `requests`
- `biopython`
- `nbib`
- `beautifulsoup`
- `openai`
- `matplotlib`
- `seaborn`
- `nltk`


## Installation

To install the GRPM System, clone the repository to your local machine:

```
git clone https://github.com/johndef64/GRPM_playground.git
```


