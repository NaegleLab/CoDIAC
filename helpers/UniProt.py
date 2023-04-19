import requests
import pandas as pd
import numpy as np
import ast

class Uniprot:
    '''Given a Uniprot ID, this class can be used to fetch reference structure attributes such as Domains, Domain Architecture, Gene name, Reference sequence from Uniprot for specie of interest and returns a Domain reference csv file'''
    
    def __init__(self, UPROT_ID):
        self.UPROT_ID = UPROT_ID
    
    # Species function needs to be added

    def Domains(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json() 
        gene = res['id']
        get_dom = []
        if 'features' in res.keys():
            for i in range(len(res['features'])):
                s = res['features'][i]
                for k, v in s.items():
                    if k == 'type':
                        if v == 'DOMAIN':
                            start = s['begin']
                            end = s['end']
                            name = s['description']
                            header = self.UPROT_ID+':'+gene+':'+name+':'+start+':'+end
                            #print(header)
                            get_dom.append(header)                    
        else:
            get_dom.append('None')
        return(get_dom)
    
    def DomainArchitecture(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json() 
        gene = res['id']
        get_dom = []
        domdict = {}
        if 'features' in res.keys():
            for i in range(len(res['features'])):
                s = res['features'][i]
                for k, v in s.items():
                    if k == 'type':
                        if v == 'DOMAIN':
                            start = s['begin']
                            end = s['end']
                            name = s['description']
                            header = self.UPROT_ID+':'+gene+':'+name+':'+start+':'+end
                            get_dom.append(header)  
                            domdict[start] = end, name
                    
        sorted_dict = dict(sorted(domdict.items(),reverse=True))
        domain_arch = []
        for key, value in sorted_dict.items():
            sort_start = key
            sort_end = value[0]
            domain = value[1]
            domain_arch.append(domain)
            
        if len(domain_arch) > 1:
            final_domarch = '|'.join(domain_arch)   
        else:
            final_domarch = domain
        #print('{count} number of domains with Architecture {Arch}'.format(count=len(domain_arch),Arch = final_domarch))  
        return(final_domarch)
    
    def Gene(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json()
        gene = res['id']
        return(gene)
    
    def RefSeq(self):
        res = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{self.UPROT_ID}').json()
        if 'sequence' in res.keys():
            seq = (res['sequence']['sequence'])
        else:
            seq = 'None'                 
        return(seq)
        
def makeRefFile(UniProt_IDs, outputFile):
    '''Makes a Domain Reference file
    
    Parameters
    ----------
        UniProt_IDs : a list of Uniprot Accession IDs
        outputFile : name of the output file
        
    Attributes
    ----------
        UPROT_Accession : Uniprot Accession ID
        Gene : gene name
        Domain Boundaries : reference domain boundary ranges with its start and end positions
        Ref Sequence : Reference sequence
        Domain Architecture : Domains found within the protein sequence arranged from N ter to C ter
        
    Returns
    -------
        .csv Domain reference file with all attributes'''
        
    Domain_list = []
    Gene_list = []
    Refseq_list = []
    Domain_Arch = []
    print('Fetch successful for Uniprot IDs :\n')
    for ID in UniProt_IDs:
        instance = Uniprot(ID)
        DOMAIN = Uniprot.Domains(instance)
        GENE = Uniprot.Gene(instance)
        REFSEQ = Uniprot.RefSeq(instance)
        DOMAIN_ARCH = Uniprot.DomainArchitecture(instance)
        
        Domain_list.append(DOMAIN)
        Gene_list.append(GENE)
        Refseq_list.append(REFSEQ)
        Domain_Arch.append(DOMAIN_ARCH)
        print(ID)
        
    df = pd.DataFrame(columns=['UPROT Accession','Gene','Domain Boundaries','Ref Sequence','Domain Architecture'])
    df['UPROT Accession'] = UniProt_IDs
    df['Gene'] = Gene_list
    df['Domain Boundaries'] = Domain_list
    df['Ref Sequence'] = Refseq_list
    df['Domain Architecture'] = Domain_Arch
    df.index = np.arange(1, len(df) + 1)
    df.to_csv(outputFile)
    print('Domain Reference File successfully created!')
    
