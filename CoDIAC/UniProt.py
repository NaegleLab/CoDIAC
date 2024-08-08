import requests
import pandas as pd
import numpy as np
import csv
from CoDIAC import InterPro
import os
from Bio import SeqIO
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


        protein_dict, error = return_uniprot_record(ID)
        if not error:   
            rowdata.append(protein_dict['uniprot_id'])
            rowdata.append(protein_dict['gene'])      
            rowdata.append(protein_dict['species'])
            rowdata.append(protein_dict['domains'])
            rowdata.append(protein_dict['sequence'])
            rowdata.append(protein_dict['PDB_IDs'])
            rowdata.append(protein_dict['domain_architecture'])
            writer.writerow(rowdata)
    f.close()
    
    print('Domain Reference File successfully created!')
    print('Adding Interpro Domains')
    df = InterPro.appendRefFile(outputfile_temp,  outputFile)
    #print('Adding Interpro Domains')
    #delete the temp file
    os.remove(outputfile_temp)
    return df
    #now add the InterPro Domains.
    #InterProDomain.appendRefFile(outputfile, outputfile)

def return_uniprot_record(uniprot_id):
    """
    For a uniprot_id, return information about the protein record in a dictionary format.

    Parameters
    ----------
    uniprot_id : string
        Uniprot ID of the protein of interest

    Returns
    -------
    protein_dict : dict
        Dictionary with keys of 'uniprot_id', 'gene', 'species', 'sequence', 'domains', 'PDB_IDs', 'domain_architecture'
        Values are the corresponding values for the protein record. Domains and PDB_IDs are string lists separated by ';' 
        domain_architecture is a string list separated by '|'
    error : int
        0 if no error, 1 if an error occurred in the request
    
    """
    protein_dict = {}
    error = 0
    protein_dict['uniprot_id'] = uniprot_id
    try:
        get_url = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}?canonical=true')
        if get_url.status_code == 200:
            response = get_url.json()
            #Gene
            if 'gene' in response.keys():
                if 'name' in response['gene'][0].keys():
                    gene = response['gene'][0]['name']['value']
                elif 'orfNames' in response['gene'][0].keys():
                    gene = response['gene'][0]['orfNames'][0]['value']
                else:
                    gene = 'None'
            else:
                gene = 'None'
            protein_dict['gene'] = gene

            #Species
            species = response['organism']['names'][0]['value']
            protein_dict['species'] = species
            #Reference sequence
            if 'sequence' in response.keys():
                seq = (response['sequence']['sequence'])
            else:
                seq = 'None'  
            protein_dict['sequence'] = seq

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
                                #replace the name values if they have ';' in them
                                name = name.replace(';', '')
                                header = name+':'+start+':'+end
                                domainheader.append(header)
            Domain = ';'.join(map(str, domainheader))
            protein_dict['domains'] = Domain

            #PDB IDs
            PDB_IDs=[]
            for i in range(len(response['dbReferences'])):

                reftype = response['dbReferences'][i]['type']
                PDB_ID = response['dbReferences'][i]['id']

                if reftype == 'PDB':
                    PDB_IDs.append(PDB_ID)
            PDBID = ';'.join(map(str, PDB_IDs))
            protein_dict['PDB_IDs'] = PDBID

            #Domain Architecture
            domdict = {}
            if 'features' in response.keys():
                for i in range(len(response['features'])):
                    s = response['features'][i]
                    for k, v in s.items():
                        if k == 'type':
                            if v == 'DOMAIN':
                                start_str = s['begin']
                                if not start_str[0].isdigit():
                                    start_str = start_str[1:]
                                start = int(start_str)

                                end_str = s['end']
                                if not end_str[0].isdigit():
                                    end_str = end_str[1:]
                                end = int(end_str)
                                name = s['description']
                                domdict[start] = end, name

            sorted_dict = dict(sorted(domdict.items(),reverse=False))
            domain_arch = []
            for key, value in sorted_dict.items():
                sort_start = key
                sort_end = value[0]
                domain = value[1]
                domain = domain.replace(';', '') #for the occasional issue of ';' in uniprot domain names
                domain_arch.append(domain)

            
            final_domarch = '|'.join(domain_arch) 
            protein_dict['domain_architecture'] = final_domarch  

                    
    except requests.exceptions.RequestException as err:
        print('ERROR:',err)
        error = 1


    return protein_dict, error


def make_domain_fasta_dict(reference_csv, Interpro_ID, n_term_offset=0, c_term_offset=0):
    """
    Given a Uniprot Reference File of proteins, which contain Interpro domain annotations, create a dictionary 
    with keys that are feasta headers for each domain of interest and the value is the domain sequence
    The n_term and c_term offsets will  build a small padding in case of domain boundary issues (default are 0, but can be set to lengths up to 20)
    This is a common function that can be used to call from different tools so that the domain boundaries, if extended, 
    are consistent across different aspects. 

    Parameters
    ----------
    reference_csv : string
        File location that contains the reference of interest (like produced from Uniprot.makeRefFile)
    Interpro_ID: string 
        Interpro ID - for example in a reference line such as 
        SH3_domain:IPR001452:82:143;SH2:IPR000980:147:246;Prot_kinase_dom:IPR000719:271:524
        the interpro ID for the SH3_domain is IPR001452; for the SH2 domain is IPR000980
    n_term_offset: int
        Number of amino acids to extend in the n-term direction (up to start of protein)
    c_term_offset: int
        Number of amino acids to extend in the c-term direction (up to end of protein)
    
    Returns
    -------
    domain_dict: dict
        key: fasta header in the format '>uniprot_ID|domain_name|domain_number|InterproID|domain_start|domain_end'
        value: sequence of the domain
    
    """
    if not isinstance(n_term_offset, int):
        print("ERROR: rounding n_term_offset to an integer")
        n_term_offset = int(n_term_offset)
    if not isinstance(c_term_offset, int):
        print("ERROR: rounding n_term_offset to an integer")
        c_term_offset = int(c_term_offset)
    # if n_term_offset < 0 or n_term_offset > 20:
    #     print("ERROR: resetting n_term_offset to a default number, should be between 0 and 20")
    #     n_term_offset = 5
    # if c_term_offset < 0 or c_term_offset> 20:
    #     print("ERROR: resetting c_term_offset to a default number, should be between 0 and 20")
    #     c_term_offset = 5
    print("n offset is %d and c offset is %d"%(n_term_offset, c_term_offset))

    domain_dict = {}
    domainRef = pd.read_csv(reference_csv)
    for row_index, row in domainRef.iterrows(): #for each protein in the reference
        domains = row['Interpro Domains']
        domainsArr = domains.split(';')
        uniprot_id = row['UniProt ID']
        gene_name = row['Gene']
        seq = row['Ref Sequence']
        species = row['Species']
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
                domain_seq = seq[new_start-1:new_end] #0 based indexing for python access of protein indexes for start, then the 1
                #uniprot_ID|gene|InterproID|domain_name|domain_number|domain_start|domain_end
                #print("DEBUG: >%s|%s|%s|%s|%d|%s|%d|%d\n"%(uniprot_id, gene_name, species, domain_name, domainNum, Interpro_ID, new_start, new_end))
                
                header = ">%s|%s|%s|%s|%d|%s|%d|%d"%(uniprot_id, gene_name, species, domain_name, domainNum, Interpro_ID, new_start, new_end)
                domain_dict[header] = domain_seq
    return domain_dict

def print_domain_fasta_file(reference_csv, Interpro_ID, output_file, n_term_offset=0, c_term_offset=0, APPEND=False):
    '''Given a Uniprot Reference File of proteins, which contain Interpro domain annotations, create a fasta 
    file of the domains of interest found in the reference file. The n_term and c_term offsets will 
    build a small padding in case of domain boundary issues (default are 0, but can be set to lengths up to 20)
    
    Parameters
    ----------
    reference_csv : string
        File location that contains the reference of interest (like produced from Uniprot.makeRefFile)
    Interpro_ID: string 
        Interpro ID - for example in a reference line such as 
        SH3_domain:IPR001452:82:143; SH2:IPR000980:147:246; Prot_kinase_dom:IPR000719:271:524
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
    domain_fasta_dict: dict
        As described in make_domain_fasta_dict - headers are fasta headers and values are the domain sequence
        Also writes to output_file

    '''

    domain_fasta_dict = make_domain_fasta_dict(reference_csv, Interpro_ID, n_term_offset=n_term_offset, c_term_offset=c_term_offset)

    if not APPEND:
        f = open(output_file, "w")
    else:
        f= open(output_file, "a")
    
    for header in domain_fasta_dict:
        f.write(header+"\n")
        f.write(domain_fasta_dict[header]+"\n")
        #header in format (">%s|%s|%s|%s|%d|%s|%d|%d\n"%(uniprot_id, gene_name, species, domain_name, domainNum, Interpro_ID, new_start, new_end))
    f.close()
    return domain_fasta_dict


def translate_fasta_to_new_headers(fasta_file, output_fasta, key_array_order):
    """
    UniProt.print_domain_fasta_file
    prints a highly informative fasta header in the production of printing fasta sequences
    for a given protein domain of interest. This structure is, in the following order, separated
    by '|'
    uniprot_id, gene_name, domain_name, domainNum, Interpro_ID, start, end
    
    If you wish to change the order, keeping track of this change, use this function to do so
    where you will use key_array_order as a list to indicate what items in what listing. 
    possible_values = ['uniprot', 'gene', 'domain_name', 'domain_num', 'Interpro_ID', 'start', 'end']

    For example, if you wanted to use uniprot ID first, gene, and the starting and ending position of the domain you would pass in for 
    key_array_order ['uniprot', 'gene', 'start', 'end']

    This will print the new fasta file at output_fasta and a mapping file, using the base of the output_fasta with _mapping.csv 

    Parameters
    ----------
    fasta_file: str
        location of input fasta file to translate
    output_fasta: str
        location to print the output fasta, using the longer headers
    key_array_order: list
        List of strings that includes the order and the values to keep from the possible values
        ['uniprot', 'gene', 'domain_name', 'domain_num', 'Interpro_ID', 'start', 'end']
    
    Returns
    -------
    output_fasta: str
        Confirmation of location of output file
    mapping_file: str
        Location of the mapping file created. This uses the same base name as the output_fasta and adds _mapping.csv

    """
    possible_values = ['uniprot', 'gene', 'species', 'domain_name', 'domain_num', 'Interpro_ID', 'start', 'end']
    values_to_keep = []
    for i in range(0, len(key_array_order)):
        item = key_array_order[i]
        if item not in possible_values:
            print("ERROR: %s not in possible values, skipping this item"%(item))
        else:
            values_to_keep.append(item)

    #now read in the fasta file and translate old headers to new. Keep this as dict too, 
    # fail and return a fatal error if naming convention produces non-unique fasta header
    header_trans_dict = {}
    with open(fasta_file, 'r') as handle, open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(handle, "fasta"):
            temp_dict = {}
            #record_id = record.id
            record_id = record.description #changed to description. record.id ends at first space, which happens in species
            record_id_vals = record_id.split('|')
            #record_id_vals = record_id.split('|(?=[^ ]))
            if len(record_id_vals) != 8:
                print("FATAL ERROR: Fasta file not formatted as expected coming from print_domain_fasta_file, has %d items"%(len(record_id_vals)))
                print(record_id)
                return
            temp_dict['uniprot'] = record_id_vals[0]
            temp_dict['gene'] = record_id_vals[1]
            temp_dict['species'] = record_id_vals[2]
            temp_dict['domain_name'] = record_id_vals[3]
            temp_dict['domain_num'] = record_id_vals[4] 
            temp_dict['Interpro_ID'] = record_id_vals[5] 
            temp_dict['start'] = record_id_vals[6] 
            temp_dict['end'] = record_id_vals[7]
            new_header_list = []
            for i in range(0, len(values_to_keep)):
                #print("DEBUG: adding %s, which is %s"%(values_to_keep[i],temp_dict[values_to_keep[i]]))
                new_header_list.append(temp_dict[values_to_keep[i]])
            #print("DEBUG: new header values are ")
            #print(new_header_list)
            new_header = "|".join(new_header_list)
            if new_header not in header_trans_dict:
                header_trans_dict[new_header] = record_id
            else:
                print("FATAL ERROR: Header request results in multiple fasta sequences with the same header. Suggest you consider adding identifying information for multiple domains from the same protein record")
                return
            modified_record = record
            modified_record.id = new_header
            modified_record.description = ''
            SeqIO.write(modified_record, out_handle, 'fasta')
        
        #now make an output file to keep track of the header matching
        mapping_file = output_fasta.strip(".fasta")
        mapping_file = mapping_file+"_mapping.csv"
        #df = pd.DataFrame.from_dict(header_trans_dict)
        #

        df = pd.DataFrame(header_trans_dict.items(), columns=['short', 'full'])
        df.to_csv(mapping_file)
        print("Created files: %s and %s"%(output_fasta, mapping_file))
        return output_fasta, mapping_file

def translate_fasta_back(fasta_file, header_trans_file, output_fasta):
    """
    Assuming that you have a fasta_file with shortened headers and would like to move those back 
    to the long form names, found in the mapping file (header_trans_file) created by translate_fasta_to_new_headers
    Use this function to print at the output_fasta location the fasta file with long headers. You would do this assuming 
    you wish to preserve a change, such as through alignment, of the shortened headers.

    Parameters
    ----------
    fasta_file: str
        location of input fasta file to translate
    header_trans_file: str
        location that stores the header translation (as written in translate_fasta_to_new_headers)
    output_fasta: str
        location to print the output fasta, using the longer headers
    
    Returns
    -------
        No returns, prints by non-append to output_fasta
    """
    
    #read in the header_translation
    trans_df = pd.read_csv(header_trans_file, index_col=0)
    if len(trans_df.columns)!=2 or ('full' not in trans_df.columns) or ('short' not in trans_df.columns):
        print("FATAL ERROR: translation file %s not formatted as expected with a short and full name"%(header_trans_file))
        return
    #a dict look up will be faster, so convert to dictionary
    trans_dict = {}
    count = 0
    count_total = 0
    for index, row in trans_df.iterrows():
        trans_dict[row['short']] = row['full']
    with open(fasta_file, 'r') as handle, open(output_fasta, 'w') as out_handle:
        for record in SeqIO.parse(handle, "fasta"):
            count_total+=1
            if record.id not in trans_dict:
                print("ERROR, cannot find %s in translation, skipping this record"%(record.id))
                break
            modified_record = record
            modified_record.id = trans_dict[record.id]
            modified_record.description = ''
            SeqIO.write(modified_record, out_handle, 'fasta')
            count+=1
    print("Complete, printed %d records out of %d to %s file"%(count, count_total, output_fasta))

        

