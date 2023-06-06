import requests
import pandas as pd
import numpy as np
import csv
#import InterProDomain

def makeRefFile(Uniprot_IDs, outputFile):
    '''Makes a Domain Reference file
    
    Parameters
    ----------
    Uniprot_IDs : list
        list of Uniprot Accession IDs generated by InterPro.py
    outputFile : string
        name of the output file
        
    Attributes
    ----------
        UniProt_ID : Uniprot Accession ID
        Gene : gene name
        Species : scientific name
        Domains : Reference domain names with boundary ranges
        Ref Sequence : Reference sequence
        PDB IDs : All PDB IDs linked to specific UniProt ID
        Domain Architecture : Domains found within the protein sequence arranged from N ter to C ter
        
    Returns
    -------
        Domain reference file in .csv format'''
    
    
    f = open(outputFile, 'w')
    writer = csv.writer(f)
    header = ['UniProt ID', 'Gene', 'Species', 'Uniprot Domains', 'Ref Sequence','PDB IDs','Uniprot Domain Architecture']
    writer.writerow(header)
    for ID in Uniprot_IDs:
        rowdata = []
        try:
            get_url = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{ID}')
            if get_url.status_code == 200:
                response = get_url.json()
                #Gene
                gene = response['gene'][0]['name']['value']

                #Species
                species = response['organism']['names'][0]['value']

                #Reference sequence
                if 'sequence' in response.keys():
                    seq = (response['sequence']['sequence'])
                else:
                    seq = 'None'  

                #Domain name along with its boundaries
                domainheader=[]
                if 'features' in response.keys():
                    for i in range(len(response['features'])):
                        s = response['features'][i]
                        for k, v in s.items():
                            if k == 'type':
                                if v == 'DOMAIN':
                                    start = s['begin']
                                    end = s['end']
                                    name = s['description']
                                    header = name+':'+start+':'+end
                                    domainheader.append(header)
                Domain = ';'.join(map(str, domainheader))

                #PDB IDs
                PDB_IDs=[]
                for i in range(len(response['dbReferences'])):

                    reftype = response['dbReferences'][i]['type']
                    PDB_ID = response['dbReferences'][i]['id']

                    if reftype == 'PDB':
                        PDB_IDs.append(PDB_ID)
                PDBID = ';'.join(map(str, PDB_IDs))

                #Domain Architecture
                domdict = {}
                if 'features' in response.keys():
                    for i in range(len(response['features'])):
                        s = response['features'][i]
                        for k, v in s.items():
                            if k == 'type':
                                if v == 'DOMAIN':
                                    start = int(s['begin'])
                                    end = s['end']
                                    name = s['description']
                                    domdict[start] = end, name

                sorted_dict = dict(sorted(domdict.items(),reverse=False))
                domain_arch = []
                for key, value in sorted_dict.items():
                    sort_start = key
                    sort_end = value[0]
                    domain = value[1]
                    domain_arch.append(domain)

               
                final_domarch = '|'.join(domain_arch)   

                    
        except requests.exceptions.RequestException as err:
            print('ERROR:',err)
            
        rowdata.append(ID)
        rowdata.append(gene)      
        rowdata.append(species)
        rowdata.append(Domain)
        rowdata.append(seq)
        rowdata.append(PDBID)
        rowdata.append(final_domarch)
        writer.writerow(rowdata)
    f.close()
    
    print('Domain Reference File successfully created!')
    #now add the InterPro Domains.
    #InterProDomain.appendRefFile(outputfile, outputfile)