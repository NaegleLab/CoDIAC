from setuptools import setup, find_packages


# Run setup
setup(
    name="codiac",
    version="1.0.0",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/CoDIAC",
    install_requires=['pandas==2.1.*', 'numpy==1.26.*', 'scipy==1.11.*', 'matplotlib==3.7.*', 'seaborn==0.13.*', 'statsmodels==0.14.*', 'biopython==1.79.*','requests==2.31.*', 'cogent3==2023.7.18a1'],
    license='GNU General Public License v3',
    description='CoDIAC: Comprehensive Domain Interface Anlysis of Contacts',
    long_description='This is the source code for CoDIAC, an open source Python toolkit for harnessing InterPro, Uniprot, PDB, and AlphaFold for generating references of proteins containing domains of interest and analyzing the contacts that exist between that domain and other regions of interest.', 
    project_urls = {'Documentation': 'https://naeglelab.github.io/CoDIAC/index.html'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.6",
    package_data={'': ['CoDIAC/data/proteomescout_everything_20190701/*.tsv']},
    zip_safe = False,
    entry_points={
    'console_scripts': [
        'run-analysis = analysis.analysis:main'
    ]
}
)
