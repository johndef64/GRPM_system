# GRPM System

GRPM (Gene-Rsid-Pmid-Mesh) system is a comprehensive tool designed to integrate and analyze genetic polymorphism data associated with specific biomedical subjects. It comprises five modules that allow data retrieval, merging, analysis, and incorporation of GWAS data.

[![medrxiv Manuscript](https://img.shields.io/badge/medrxiv-10.1101/2023.08.04.23293659-blue.svg)](https://www.medrxiv.org/content/10.1101/2023.08.04.23293659v1.full.pdf+html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8205724.svg)](https://doi.org/10.5281/zenodo.8205724)

## Overview

- [Introduction](#introduction)
- [Modules](#modules)
- [Updates](#updates)
- [Installation](#installation)
- [Usage](#usage)
- [Requirements](#requirements)

## Introduction

GRPM System aims to build a comprehensive dataset of human genetic polymorphisms associated with nutrition. By combining data from multiple sources and utilizing MeSH terms as a framework, the system enables researchers to explore the vast genetic literature in search of variants significantly associated with a specific biomedical subject.
The main purpose of developing this resource was to assist nutritionists in investigating gene-diet interactions and implementing personalized nutrition interventions.

![Graphical Abstract](misc_data/graphical_abstract_s.png)

## Modules

The GRPM System comprises five modules that perform various tasks to facilitate the integration and analysis of genetic polymorphism data associated with nutrition. These modules are as follows:

To try out GRPM System. Run each module separately by clicking the "Open in Colab". Be careful to import all necessary dependencies and files. Google Drive folder synch option available.

| No. | Notebook | Module | Description |
| --- | --- | --- |-------------|
| 1. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_system/blob/main/GRPM_01_dataset_builder.ipynb) | [Dataset Builder](https://github.com/johndef64/GRPM_system/blob/main/GRPM_01_dataset_builder.ipynb)| Retrieves data from LitVar and PubMed databases, merging them into a CSV format.
| 2. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_system/blob/main/GRPM_02_ref-mesh_builder.ipynb) | [Reference Mesh List Builder](https://github.com/johndef64/GRPM_system/blob/main/GRPM_02_ref-mesh_builder.ipynb)| Generates a coherent MeSH term list using the ChatGPT language model and the OpenAI API for exploring the GRPM dataset.
| 3. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_system/blob/main/GRPM_03_dataset_survey.ipynb) | [GRPM Dataset Screening](https://github.com/johndef64/GRPM_system/blob/main/GRPM_03_dataset_survey.ipynb)| Integrates the MeSH term list into the GRPM dataset and extracts a survey for comprehensive analysis.
| 4. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_system/blob/main/GRPM_04_grpm-data_analyzer.ipynb) |  [GRPM Reports and Data Analyzer](https://github.com/johndef64/GRPM_system/blob/main/GRPM_04_grpm-data_analyzer.ipynb) |Analyzes reports and GRPM association data, utilizing `matplotlib` and `seaborn` for data visualization. 
| 5. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/johndef64/GRPM_system/blob/main/GRPM_05_gwas_grpm_merger.ipynb) | [Merge GWAS and GRPM Data](https://github.com/johndef64/GRPM_system/blob/main/GRPM_05_gwas_grpm_merger.ipynb) | Integrates GWAS data from the complete catalog, associating GWAS phenotypes and potential risk/effect alleles with the GRPM Dataset. 


![GRPM system: Integrating Genetic Polymorphism Data with PMIDs and MeSH Terms to Retrieve Genes and rsIDs for Biomedical Research Fields. GRPM Dataset: pcg, protein coding genes; rna, RNA genes; pseudo, presudogenes; in parentheses, dataset shape.](misc_data/grpm_system.png)

These modules provide a comprehensive framework for researchers and nutritionists to explore genetic polymorphism data and gain insights into gene-diet interactions and personalized nutrition interventions.



## Updates

The GRPM Dataset available on Zenodo is a snapshot of [LitVar1](https://www.ncbi.nlm.nih.gov/CBBresearch/Lu/Demo/LitVar/help.html). LitVar1 is now <u>**deprecated**</u> and has been fully replaced by [LitVar2](https://www.ncbi.nlm.nih.gov/research/litvar2/). Module 1 ([Dataset Builder](https://github.com/johndef64/GRPM_system/blob/main/GRPM_01_dataset_builder.ipynb)) has been updated to retrieve data from LitVar2. The subsequent modules in the pipeline remain functional and can be tested using the original version of the GRPM Dataset available on Zenodo.

## Installation

To install GRPM System, clone the repository to your local machine:

```
git clone https://github.com/johndef64/GRPM_system.git
```
Otherwise, run each module separately in Google Colab importing Google Drive to keep-up your progress.


## Usage

Detailed instructions on how to use each module of  GRPM System can be found inside the relative Jupyter Module provided in the repository. Make sure to follow the instructions and install the necessary Python packages specified for each module.


## Requirements

GRPM System has the following requirements:

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









