import pandas as pd
import requests
import logging
import os
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import csv
import shutil

def makeRefFile(refUniprot_file, cifFile_path, outputFile):
    '''make a structure reference file for alphafodl structures
    Parameters
    ----------
        refUniprot_file : str
            input uniprot reference file
        cifFile_path : str
            input path to find mmcif files downloaded for AlphaFold structures
    Returns
    -------
        outputFile : str
            file path for output reference csv file

    '''
        
    df = pd.read_csv(refUniprot_file)
    uniprot_ids = df['UniProt ID'].tolist()
    col_headers = ['PDB_ID','ENTITY_ID', 'CHAIN_ID', 'pdbx_description', 'rcsb_gene_name', 'pdbx_gene_src_scientific_name','pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can', 'rcsb_sample_sequence_length', 'modifications locations','modifications', 'rcsb_entity_polymer_type', 'macromolecular_type', 'database_name', 'database_accession','rcsb_uniprot_protein_sequence','entity_beg_seq_id', 'ref_beg_seq_id', 'aligned_regions_length', 'mutations exist (Y/N)','rcsb_mutation_count', 'mutations locations', 'pdbx_mutation joined', 'molecular_weight', 'experimental_method', 'resolution_combined','deposit_date', 'audit_author_name', 'title', 'pdbx_database_id_doi', 'pdbx_database_id_pub_med', 'ref:gene name', 'ref:struct/ref sequence','ref:reference range','ref:start position offset','Gaps ref:struct','Gaps struct:ref','ref:variants', 'ref:domains','ref:struct domain architecture','ref:full protein domain','ref:protein domain architecture',
 'domains_AF', 'struct_dom_arch_AF', 'pro_dom_arch_AF']
    print(col_headers)
    
    with open(outputFile,'w') as f:
        writer = csv.writer(f)
        header = col_headers
        writer.writerow(header)
        for uniprot_id in uniprot_ids:
            print(uniprot_id)
            data = []
            PDB_ID = 'AF-'+uniprot_id+'-F1'
            data.append('AF-'+uniprot_id+'-F1')
            # print(PDB_ID)
            
            cif_file = MMCIF2Dict(cifFile_path+PDB_ID+'/'+PDB_ID+'.cif')
            data.append(cif_file['_entity_poly.entity_id'][0])
            data.append(cif_file['_entity_poly.pdbx_strand_id'][0])
    
            get_url = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}')
            response = get_url.json()
            gene = response['gene'][0]['name']['value']
            print(PDB_ID,'-', gene)
            data.append(response['protein']['recommendedName']['fullName']['value'])
            data.append(response['gene'][0]['name']['value'])
            data.append('Homo sapiens')
            data.append(response['sequence']['sequence'])
            data.append(response['sequence']['sequence'])
            data.append(response['sequence']['length'])
            data.append('')
            data.append('')
            data.append('Protein')
            data.append(response['proteinExistence'])
            data.append('AlphaFold')
            data.append(response['accession'])
            data.append(response['sequence']['sequence'])
    
            for entry in range(len(response['features'])):
                if response['features'][entry]['type'] == 'CHAIN':
                    data.append(response['features'][entry]['begin'])
                    data.append(response['features'][entry]['begin'])
    
            data.append(response['sequence']['length'])
            data.append('N')
            data.append(0)
            data.append('')
            data.append('')
            data.append(response['sequence']['mass'])
            data.append('AlphaFold')
            data.append(response['sequence']['version'])
            data.append(response['sequence']['modified'])
            data.append('NA')
            data.append('NA')
            data.append('NA')
            data.append('NA')
            data.append(response['gene'][0]['name']['value'])
            data.append(response['sequence']['sequence'])
            data.append(str('1'+'-'+str(response['sequence']['length'])))
            data.append(0.0)
            data.append(0.0)
            data.append(0.0)
            data.append('')
            data.append(reference_arch(refUniprot_file,gene)[0])
            data.append(reference_arch(refUniprot_file,gene)[1])
            data.append(reference_arch(refUniprot_file,gene)[1])
            domain_list = []
            struct_arch = []
            for entry in range(len(response['features'])):          
                if response['features'][entry]['type'] == 'DOMAIN':
                    domain_entry = (response['features'][entry]) 
    
                    domain_name = domain_entry['description']
                    domain_start = domain_entry['begin']
                    domain_end = domain_entry['end']
    
                    domain_range = domain_name+':IPR000980:'+domain_start+','+domain_end+',0,0'
                    domain_list.append(domain_range)
                    struct_arch.append(domain_name)
            separator1 = ';'
            separator2 = '|'
            data.append(separator1.join(domain_list))
            data.append(separator2.join(struct_arch))
            data.append(separator2.join(struct_arch))
            # print(gene, reference_arch(gene)[0], reference_arch(gene)[1])
            writer.writerow(data)
            data.clear()



def reference_arch(refUniprot_file, gene):
    '''This retrieves the interpro domain architecture for a specific gene '''
        
    df = pd.read_csv(refUniprot_file)
    identical_domains=[]
    arch_with_index = []
    ref_arch = df.loc[df['Gene']==gene, ['Interpro Domain Architecture']].values[0][0]
    ref_ranges = df.loc[df['Gene']==gene, ['Interpro Domains']].values[0][0]

    parse1 = ref_ranges.split(';')
    update_ranges=[]
    for i in parse1:
        sh2, ipr, start, end = i.split(':')
        new_str = sh2+':'+ipr+':'+str(start)+','+str(end)+',0,0'
        update_ranges.append(new_str)
    return ';'.join(update_ranges), ref_arch


def get_cifFile(refUniprot_file, outputDir):
    '''Downloads .cif files for AlphaFold structures using database version 2
    
    Parameters
    ----------
        refUniprot_file : cstr
            uniprot reference file input
        
    Returns
    -------
        outputDir : str
            Location to store the .cif files
    
    '''
    # inputFile = 'UniProt_IDs.csv'
    # df = pd.read_csv(inputFile)
    df = pd.read_csv(refUniprot_file)
    uniprot_ids = df['UniProt ID'].tolist()
    
    for i in uniprot_ids:
        Uniprot_id = i
        alphafold_id = 'AF-'+Uniprot_id+'-F1'
        path = outputDir +alphafold_id
        database_ver = 'v2'
        model_url = f'https://alphafold.ebi.ac.uk/files/{alphafold_id}-model_{database_ver}.cif'
        os.system(f'curl {model_url} -o {alphafold_id}.cif')
        os.mkdir(path)
        # new_cifpath = path+'/'+alphafold_id+'.cif'
        shutil.move('./'+alphafold_id+'.cif', path)
