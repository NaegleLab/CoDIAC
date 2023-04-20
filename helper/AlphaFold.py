import os
import pandas as pd
import requests

def get_CifFile(inputFile, outputDir):
    '''Download .cif files for AlphaFold structures 
    
    Parameters
    ----------
        inputFile : csv file containing a list of UniProt IDs generated using InterPro.py
        outputDir : Location to store the .cif files
        
    Returns
    -------
        .cif files for AlphaFold structures
    
    '''
    inputFile = 'UniProt_IDs.csv'
    df = pd.read_csv(inputFile)

    for i in df['Uniprot ID']:
        Uniprot_id = i
        alphafold_id = 'AF-'+Uniprot_id+'-F1'
        path = outputDir + alphafold_id
        database_ver = 'v2'
        model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_{database_ver}.cif'
        os.system(f'curl {model_url} -o {path}.cif')