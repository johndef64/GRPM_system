## GRPM System

The GRPM (Gene-RsID-PubMedID-MeSH) system is a comprehensive tool designed to integrate and analyze genetic polymorphism data associated with nutrition. It consists of five modules that streamline the process of data retrieval, merging, analysis, and incorporation of GWAS (Genome-Wide Association Study) data. 

### Module 01: Database Builder

The Database Builder module retrieves data from the LitVar and PubMed databases and merges them into a CSV format. It uses Python modules such as `requests`, `pandas`, `biopython`, `nbib`, and `beautifulsoup` to extract, manipulate, and parse data.

### Module 02: Reference MeSH Term List Builder

The Reference MeSH Term List Builder module creates a coherent MeSH (Medical Subject Headings) term list to explore the database. It utilizes the ChatGPT language model and interacts with the OpenAI API. The Python modules used are `pandas` and `openai`.

### Module 03: Database and Survey Integration

The Database and Survey Integration module incorporates the MeSH term list generated in Module 02 into the database. It extracts a survey from the integrated data, providing a comprehensive dataset for further analysis. The required Python module for this module is `pandas`.

### Module 04: Analysis of Reports and GRPM Association Data

The Analysis of Reports and GRPM Association Data module analyzes the reports generated in Module 01 and Module 03, as well as the GRPM association data. It utilizes Python modules such as `pandas`, `matplotlib`, and `seaborn` for data manipulation and data visualization.

### Module 05: Incorporation of GWAS Data

The Incorporation of GWAS Data module integrates GWAS data from the complete GWAS catalog dataset into the GRPM surveys. It associates GWAS phenotypes and potential risk/effect alleles with the GRPM relationships. The required Python modules for this module are `pandas` and `Natural Language Toolkit (NLTK)`.

The GRPM System provides a comprehensive framework for researchers and nutritionists to explore genetic polymorphisms associated with nutrition. By integrating data from multiple sources and employing MeSH terms as a framework, it enables in-depth analysis of gene-diet interactions. The incorporation of GWAS data further enhances the system's functionality and potential for personalized nutrition interventions. 

For instructions on using the GRPM System and the specific Python modules required for each module, please refer to the documentation provided in the repository.
