from setuptools import setup, find_packages


# Run setup
setup(
    name="codiac",
    version="0.0.1",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/CoDIAC",
    install_requires=['pandas', 'numpy', 'scipy', 'matplotlib', 'seaborn', 'statsmodels', 'biopython','requests'],
    license='GNU General Public License v3',
    description='CoDIAC: COmprehensive Domain Interface Anlysis of Contacts',
#    long_description=This is the source code for CoDIAC, an open source Python toolkit for harnessing InterPro, Uniprot, PDB, and AlphaFold for generating references of proteins containing domains of interest and analyzing the contacts that exist between that domain and other regions of interest."""   
    project_urls = {'Documentation': 'https://naeglelab.github.io/CoDIAC/index.html'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.6",
    zip_safe = False
)
