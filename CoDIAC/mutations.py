from selenium import webdriver
from selenium.webdriver import Chrome
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException
import shutil
import glob
import os
import re
import pandas as pd
from Bio import SeqIO
from pybiomart import Dataset

clinvar_sig = ['Benign', 'Likely benign','Likely pathogenic', 'Pathogenic','Likely pathogenic, low penetrance','Pathogenic, low penetrance',
               'Likely risk allele']

def makeMutationFeafile(fastafile, downloads_path, csvfiles_dir, output_feafile):
    '''make a feature file with mutations recorded as features. mutations extracted from GnomAD using Uniprot ID and their corresponding Ensemble ID.
    Parameters
    ----------
        fastafile : str
            the reference fasta file that is created using UniProt.py and key_array_order= ['uniprot', 'gene', 'domain_num', 'start', 'end']
        downloads_path : str
            the path to the downloads folder on the device where files downloaded are placed
        csvfiles_dir : str
            path of the directory where we would like to move and save the downloaded variant csv files from the downloads folder
        output_feafile : str
            the path to save the output feature file

    Returns
    -------
        output_feafile : str
            Mutation feature file (.fea) '''
    
    fetch_VariantCSV(fastafile, csvfiles_dir, downloads_path)
    
    listfiles = os.listdir(csvfiles_dir)
    
    with open(output_feafile,'w') as file:
        for csvfile in listfiles:
            csv_path = csvfiles_dir + csvfile
            if csvfile.startswith('gnomAD'):
                ensembleID = csvfile.split('_')[2]
                
                uniprot_dict = reference_Dict(fastafile)
                
                for key, values in uniprot_dict.items():
                    if ensembleID == values[1]:
                        header = values[0]
                        
                start = int(header.split('|')[3])
                end = int(header.split('|')[4])
        
                df = pd.read_csv(csv_path)
                df_missense = df.loc[df['VEP Annotation']=='missense_variant']
                df_sig = df_missense.loc[df_missense['ClinVar Clinical Significance'].isin(clinvar_sig)]
                mutation = df_sig['HGVS Consequence'].tolist()
                significance = set(df_missense['ClinVar Clinical Significance'].tolist())
                
                for m in mutation:
                    mut_resid = (re.findall(r'\d+', m))
              
                    if int(mut_resid[0]) in range(start, end+1):
                        resid_upd = int(mut_resid[0]) - start + 1
                        file.write('mut\t'+str(header)+'\t-1\t'+str(resid_upd)+'\t'+str(resid_upd)+'\t'+'mut\n')
    print('Created Feature file with mutations')


def get_transcriptID(uniprot_id):
    '''Fetch transcript ID that is associated to a specific UniProt ID
    
    Parameters
    ----------
        uniprot_id : str
            UniProt accession ID
    Returns
    -------
        ID : str
            Returns a transcript ID '''
    
    dataset = Dataset(name='hsapiens_gene_ensembl',host='http://www.ensembl.org')
    df = dataset.query(attributes=['ensembl_transcript_id','ensembl_gene_id', 'external_gene_name', 'description','uniprotswissprot'])
    df_filter = df.loc[df['UniProtKB/Swiss-Prot ID'] == uniprot_id]
    ID = (df_filter['Gene stable ID'].unique().tolist())[0]
    
    return ID
    
def reference_Dict(fastafile):
    '''Generates a dictionary with fasta headers and transcript IDs for every uniprot ID found in the input fasta file.

    Parameters
    ----------
        fastafile : str
            the reference fasta file that is created using UniProt.py and key_array_order= ['uniprot', 'gene', 'domain_num', 'start', 'end']
    Returns
    -------
        uniprot_dict : dict
            dictionary whose keys are uniprot IDs and values are lists that holds the corresponding fasta headers and transcript IDs. '''
    uniprot_dict = {}
    file = SeqIO.parse(open(fastafile), 'fasta')
    for fasta in file:
        name, sequence = fasta.id, str(fasta.seq)
        uniprot_id, gene, sh2_index, start, end = name.split('|')
        uniprot_dict[uniprot_id] = [name] 
    
    for uniprotID in uniprot_dict:
        transcriptID = get_transcriptID(uniprotID)
        uniprot_dict[uniprotID].append(transcriptID)
    
    return uniprot_dict

def check_exists_by_xpath(driver):
    '''Checks whether GnomAD has mutations recorded for specific transcript ID. For the ones that do not have any information available, it pops up 'Gene not found' error. And we use this to identify Uniprot IDs/genes for which data doesnt exist on GnomAD.
    Parameters
    ----------
        driver : str
            input the created Chrome browser object
    Returns
    -------
        boolean value '''
    try:
        driver.find_element(By.XPATH, "//*[contains(text(), 'Gene not found')]")
    except NoSuchElementException:
        return False
    return True
    
def fetch_VariantCSV(fastafile, csvfiles_dir, downloads_path):
    '''fetches variant CSV files from GnomAD for every UniProt ID/transcript ID found in the input fastafile.
    Parameters
    ----------
        fastafile : str
            the reference fasta file that is created using UniProt.py and key_array_order= ['uniprot', 'gene', 'domain_num', 'start', 'end']
        csvfiles_dir : str
            path of the directory where we would like to move and save the downloaded variant csv files from the downloads folder
        downloads_path : str
            the path to the downloads folder on the device where files downloaded are placed
    Returns
    -------
        .csv files downloaded from GnomAD '''
    
    uniprot_dict = reference_Dict(fastafile)
    print('Exported Variant CSV files for transcript ID ...')
    for entry in uniprot_dict:
        ensembleID = uniprot_dict[entry][1]
        options = webdriver.ChromeOptions()
        options.page_load_strategy = "none"
        driver = Chrome(options=options)
        driver.implicitly_wait(60)
        url = "https://gnomad.broadinstitute.org/gene/"+ensembleID+"?dataset=gnomad_r4"
        driver.get(url)
        myReturnValue = check_exists_by_xpath(driver)
        if (myReturnValue == False):
            button = driver.find_element(By.XPATH, "//*[contains(text(), 'Export variants to CSV')]")
            button.click()
        else:
            print('Gene not found for transcriptID',ensembleID)
            
        driver.close()
        print(ensembleID)

    changeDirectory(csvfiles_dir, downloads_path)
    print('Variant CSV files moved from %s to %s'%(downloads_path,csvfiles_dir))

def changeDirectory(csvfiles_dir, downloads_path):
    '''Change the directory of the downloaded files'''
        
    listfiles = os.listdir(downloads_path)
    variant_csvfiles = []
    for file in listfiles:
        if file.startswith('gnomAD'):
            variant_csvfiles.append(file)
            shutil.move(downloads_path+'/'+file, csvfiles_dir)
            