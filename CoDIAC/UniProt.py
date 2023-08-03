import requests
import pandas as pd
import numpy as np
import csv
from CoDIAC import InterproDomain
import os
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
    
    outputfile_temp = outputFile+'.temp'
    f = open(outputfile_temp, 'w')
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
    print('Adding Interpro Domains')
    df = InterproDomain.appendRefFile(outputfile_temp,  outputFile)
    #print('Adding Interpro Domains')
    #delete the temp file
    os.remove(outputfile_temp)
    return df
    #now add the InterPro Domains.
    #InterProDomain.appendRefFile(outputfile, outputfile)

def print_domain_fasta_file(reference_csv, Interpro_ID, output_file, n_term_offset=0, c_term_offset=0, APPEND=False):
    """
    Given a Uniprot Reference File of proteins, which contain Interpro domain annotations, create a fasta 
    file of the domains of interest found in the reference file. The n_term and c_term offsets will 
    build a small padding in case of domain boundary issues (default are 0, but can be set to lengths up to 20)
    
    Parameters
    ----------
    reference_csv : string
        File location that contains the reference of interest (like produced from Uniprot.makeRefFile)
    Interpro_ID: string 
        Interpro ID - for example in a reference line such as 
        SH3_domain:IPR001452:82:143;SH2:IPR000980:147:246;Prot_kinase_dom:IPR000719:271:524
        the interpro ID for the SH3_domain is IPR001452; for the SH2 domain is IPR000980
    output_file: string
        location of output fasta. Fasta headers will be uniprot_ID|InterproID|domain_name|domain_number|domain_start|domain_end
        domain number will indicate in proteins with more than one domain of the same type, the occurrence of this domain
        from N-to-C numbering. Domain start and stop are relative to ones-based counting of the Uniprot
        reference protein
    n_term_offset: int
        Number of amino acids to extend in the n-term direction (up to start of protein)
    c_term_offset: int
        Number of amino acids to extend in the c-term direction (up to end of protein)
    APPEND: bool
        if APPEND=true then it will open an existing file and append, else it overwrites
    
    Returns
    -------
    No returns, prints a fasta file as described above. Fasta headers will be 
    '>uniprot_ID|domain_name|domain_number|InterproID|domain_start|domain_end

    """
    if n_term_offset < 0 or n_term_offset > 20:
        print("ERROR: resetting n_term_offset to a default number, should be between 0 and 20")
        n_term_offset = 5
    if c_term_offset < 0 or c_term_offset> 20:
        print("ERROR: resetting c_term_offset to a default number, should be between 0 and 20")
        c_term_offset = 5
    print("n offset is %d and c offset is %d"%(n_term_offset, c_term_offset))

    if not APPEND:
        f = open(output_file, "w")
    else:
        f= open(output_file, "a")

    domainRef = pd.read_csv(reference_csv)
    for row_index, row in domainRef.iterrows(): #for each protein in the reference
        domains = row['Interpro Domains']
        domainsArr = domains.split(';')
        uniprot_id = row['UniProt ID']
        gene_name = row['Gene']
        seq = row['Ref Sequence']
        #since there could be more than one value of array 
        domainNum = 0
        for index in range(0, len(domainsArr)):
            domainVals = domainsArr[index]
            domain_name, interpro, start, end = domainVals.split(':')

            if interpro == Interpro_ID:
                domainNum+=1
                new_start = int(start) - n_term_offset
                if new_start < 1:
                    new_start = 1 
                    print("NOTICE: %s n-term offset capped at start of protein"%(ID))
                new_end = int(end) + c_term_offset
                if new_end > len(seq):
                    new_end = len(seq)
                    print("NOTICE: %s c-term offset capped at end of protein"%(ID))
                domain_seq = seq[new_start-1:new_end-1] #0 based indexing for python access of protein indexes
                #uniprot_ID|gene|InterproID|domain_name|domain_number|domain_start|domain_end
                f.write(">%s|%s|%s|%d|%s|%d|%d\n"%(uniprot_id, gene_name, domain_name, domainNum, Interpro_ID, new_start, new_end))
                f.write(domain_seq+"\n")
        #if EDIT:
        #    domainRef.at[row_index, 'domains_edited'] = ";".join(domainsArr)

    f.close()
    #domainRef.to_csv(domainRefOutputFile)


