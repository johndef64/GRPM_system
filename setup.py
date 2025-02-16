from setuptools import setup, find_packages


setup(
    name='pygrpm',
    version='0.1',
    packages=find_packages(),  # Automatically find packages in your project
    install_requires=[
        # 'os', 'io', 'sys', 'datetime', and 'zipfile' are standard libraries and do not need to be listed
        'requests',
        'pyarrow',
        'numpy',
        'pandas',
        'seaborn',
        'matplotlib',
        'tqdm',
        'gdown',
        'nltk',
        'pyperclip',
        'sentence-transformers',
        'bs4'
    ],

    #entry_points={             # Optional: command-line scripts
    #    'console_scripts': [
    #        'your-command=your_package.module:function',
    #    ],
    #},
    author='JohnDef64',
    author_email='giovannimaria.defilippis@unina.it',
    description='The GRPM System is a Python-based framework capable of constructing a comprehensive dataset of human genetic polymorphisms associated with nutrition. By integrating data from multiple sources and utilizing MeSH ontology as an organizational structure, this workflow enables researchers to investigate genetic variants with significant associations to specified biomedical subjects. The primary objective of developing this resource was to support nutritionists in exploring gene-diet interactions and implementing personalized nutrition strategies.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/johndef64/GRPM_system',
    classifiers=[  # Optional: supply classifiers for search indexing
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',  # Specify required Python versions
)