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
from xml.dom import minidom
import requests
from CoDIAC import featureTools


clinvar_sig = ['Benign', 'Likely benign','Likely pathogenic', 'Pathogenic','Likely pathogenic, low penetrance','Pathogenic, low penetrance',
               'Likely risk allele']

AA_dict = {"ALA":'A',
                  "ARG":'R',
                  "ASN":'N',
                  "ASP":'D',
                  "CYS":'C',
                  "GLU":'E',
                  "GLN":'Q',
                  "GLY":'G',
                  "HIS":'H',
                  "ILE":'I',
                  "LEU":'L',
                  "LYS":'K',
                  "MET":'M',
                  "PHE":'F',
                  "PRO":'P',
                  "SER":'S',
                  "THR":'T',
                  "TRP":'W',
                  "TYR":'Y',
                  "VAL":'V'}

def gnomAD_mutations(fastafile, downloads_path, csvfiles_dir, output_feafile, N_offset=0, C_offset=0):
    '''Mutations found on the domain of interest and are extracted from GnomAD using Uniprot ID and their corresponding Ensemble ID are stored in feature files
    
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
        N_offset : int
            number of residues to add or remove from the N terminal end of the domain boundary
        C_offset : int
            number of residues to add or remove from the C terminal end of the domain boundary
 
    '''
    
    fetch_VariantCSV(fastafile, csvfiles_dir, downloads_path)
    
    listfiles = os.listdir(csvfiles_dir)
    
    header_EnsmblID_dict = ensembleID_dict(fastafile)
    tmp=[]
    with open(output_feafile,'w') as file:
        file.write('gnomAD_Mutation\t882255\n')
        for csvfile in listfiles:
            csv_path = csvfiles_dir + csvfile
            if csvfile.startswith('gnomAD'):
                ensembleID = csvfile.split('_')[2]
        
                for key, values in header_EnsmblID_dict.items():
                    if ensembleID == key:
                        for fasta_header in values:
                            header = fasta_header
                            start = int(header.split('|')[3])
                            end = int(header.split('|')[4])
                            print(key)
                           
                            df = pd.read_csv(csv_path)
                            df_missense = df.loc[df['VEP Annotation']=='missense_variant']
                            df_sig = df_missense.loc[df_missense['ClinVar Clinical Significance'].isin(clinvar_sig)]
                            mutation = df_sig['HGVS Consequence'].tolist()
                            significance = set(df_missense['ClinVar Clinical Significance'].tolist())
                            
                            for m in mutation:
                                mut_resid = (re.findall(r'\d+', m))
                          
                                if int(mut_resid[0]) in range(start+N_offset, end+1+C_offset):
                                    resid_upd = int(mut_resid[0]) - start + 1
                                    if header+':'+str(resid_upd) not in tmp:
                                        tmp.append(header+':'+str(resid_upd))
                                        file.write('gnomAD_Mutation\t'+str(header)+'\t-1\t'+str(resid_upd)+'\t'+str(resid_upd)+'\t'+'gnomAD_Mutation\n')
    print('GnomAD mutation feature file created!')
    

def OMIM_mutations(uniprot_refFile, api_key, ref_fastaFile, output_featureFile, domain_of_interest,N_offset=0, C_offset=0): 
    """Generates a feature file whose features are mutations found on the domain of interest and are extracted from OMIM database

    Parameters
    ----------
        uniprot_refFile : str
            input a uniprot reference file to get a list of uniprot IDs
        api_key : str
            this key needs to be generated through OMIM to be able to access their programmatic interface 
        ref_fastaFile : str
            fasta file with reference sequences of the domains
        domain_of_interest : str
            mutations on domains of interest (SH2 domains)
        N_offset : int
            number of residues to add or remove from the N terminal end of the domain boundary
        C_offset : int
            number of residues to add or remove from the C terminal end of the domain boundary
            
    """
        
    df_ref = pd.read_csv(uniprot_refFile)
    uniprot_ids = df_ref['UniProt ID'].tolist() 
    mutation_dict = omim_mutation_dict(uniprot_ids, api_key)
    
    with open(output_featureFile,'w') as file:
        file.write('OMIM_Mutation\tCC6677\n')
        for entry in mutation_dict:
            
            uniprot_id = entry
            gene = mutation_dict[entry][0]
            mutations = mutation_dict[entry][1]
            
            sh2_mutations={} 
            if gene != 'NA' and mutations != '':
                mut_pos = mutations.split(',')
                
                interpro_domains = df_ref.loc[df_ref['Gene']==gene, 'Interpro Domains'].values[0].split(';')
        
                for i in interpro_domains:
                    domain = i.split(':')[0]
                    start = int(i.split(':')[2])
                    end = int(i.split(':')[3])
                    query_start = int(start)+N_offset
                    query_stop = int(end)+1+C_offset
                    tmp_fea=[]
                    if domain_of_interest in domain:
                        # print(i, start, end)
                        for mut in mut_pos:
                            # print(mut)
                            if int(mut) in range(query_start, query_stop):
                                feature_pos = int(mut) - int(start) + 1
                                tmp_fea.append(feature_pos)
                        sh2_mutations[i] = tmp_fea
        
            for entry in sh2_mutations:
                
                # target = gene+'|'+entry.split(':')[2]+'|'+entry.split(':')[3]
                ref_header= makeheader(gene, int(entry.split(':')[2]), int(entry.split(':')[3]), ref_fastaFile)
                for pos in sh2_mutations[entry]:
                    file.write('OMIM_Mutation\t'+ref_header+'\t-1\t'+str(pos)+'\t'+str(pos)+'\tOMIM_mutation\n')
            
            # print(ref_header, gene, sh2_mutations,'\n')
    print('OMIM mutation feature file created!')       

def PDB_mutations(PDB_refFile,ref_fastaFile, output_featureFile, domain_of_interest,N_offset=0, C_offset=0):
    '''Generates a feature file with mutations found on the domain of interest within PDB structures
    
    Parameters
    ----------
        PDB_refFile : str
            input PDB reference file path
        ref_fastaFile : str
            fasta file with reference sequences of the domains
        output_featureFile : str
            output feature file path
        domain_of_interest : str
            mutations on domains of interest (SH2 domains)
        N_offset : int
            number of residues to add or remove from the N terminal end of the domain boundary
        C_offset : int
            number of residues to add or remove from the C terminal end of the domain boundary
            
    Returns
    -------
        feature file with mutations/variants reported in PDB structures that are present within the domain of interest '''
    df = pd.read_csv(PDB_refFile)
    j = 1
    tmp_list = []
    with open(output_featureFile,'w') as file:
        file.write('PDB_Mutation\t44AA99\n')
        for i in range(len(df)):
            variants = df['ref:variants'][i]
            uniprotid = df['database_accession'][i]
            genename = df['ref:gene name'][i]
            pdb_mutation = df['pdbx_mutation joined'][i]
            if isinstance(variants, str) and isinstance(pdb_mutation,str):
                if 'SH2' in str(df['ref:struct domain architecture'][i]):
                    if '-' not in variants:
                        if len(pdb_mutation.split(',')) == len(variants.split(';')):
                            print(df['PDB_ID'][i])
                            variant_list = variants.split(';')
                            variant_pos = []
                            for var in variant_list:
                                posnum = re.findall('\d+', var)
                                variant_pos.append(int(posnum[0]))
                            # print(variant_pos, df['PDB_ID'][i])
                            domains = df['ref:domains'][i].split(';')
                            domain_dict = {}
                            # for dom in domains:
                            #     domname, iprid, ranges = dom.split(':')
                            #     start, stop, gap, mut = ranges.split(',')
                            #     if domain_of_interest in domname:
                            #         domain_dict[domname] = [start, stop]
                            dom_index = 1
                            for dom in domains:
                                domname, iprid, ranges = dom.split(':')
                                start, stop, gap, mut = ranges.split(',')
                                if domain_of_interest in domname:
                                    domain_dict[dom_index] = [domname, start, stop]
                                dom_index += 1
                                
                            for entry in domain_dict:
                                domain_id, (start), (stop) = domain_dict[entry]
                                for varpos in variant_pos:
                                    if varpos in range(int(start)+N_offset, int(stop)+1+C_offset):
                                        header = makeheader(genename, int(start), int(stop),ref_fastaFile)
                                        feature_pos = int(varpos) - int(start) + 1
                                        # print(j,df['PDB_ID'][i], varpos, start, stop, genename, header, feature_pos)
                                        check_entry = header+':'+str(feature_pos)
                                        if check_entry not in tmp_list:
                                            tmp_list.append(check_entry)
                                            file.write('PDB_Mutation\t'+str(header)+'\t-1\t'+str(feature_pos)+'\t'+str(feature_pos)+'\tPDB_Mutation\n')
                                        j+=1
    
    
                        if len(pdb_mutation.split(',')) != len(variants.split(';')):
                            print(df['PDB_ID'][i], variants)
                            variant_list = variants.split(';')
                            variant_pos = []
                            for var in variant_list:
                                posnum = re.findall('\d+', var)
                                variant_pos.append(int(posnum[0]))
                           
                            domains = df['ref:domains'][i].split(';')
                            domain_dict = {}
                            dom_index = 1
                            for dom in domains:
                                domname, iprid, ranges = dom.split(':')
                                start, stop, gap, mut = ranges.split(',')
                                if domain_of_interest in domname:
                                    domain_dict[dom_index] = [domname, start, stop]
                                dom_index += 1
                            print(domain_dict)
                            for entry in domain_dict:
                                domain_id, (start), (stop) = domain_dict[entry]
                                for varpos in variant_pos:
                                    if varpos in range(int(start)+N_offset, int(stop)+1+C_offset):
                                        header = makeheader(genename, int(start), int(stop),ref_fastaFile)
                                        feature_pos = int(varpos) - int(start) + 1
                                        # print(j,df['PDB_ID'][i], varpos, start, stop, genename, header, feature_pos,'\n')
                                        check_entry = header+':'+str(feature_pos)
                                        if check_entry not in tmp_list:
                                            tmp_list.append(check_entry)
                                            file.write('PDB_Mutation\t'+str(header)+'\t-1\t'+str(feature_pos)+'\t'+str(feature_pos)+'\tPDB_Mutation\n')
                                        # else:
                                        #     tmp_list.append(check_str)
                                        j+=1
    print('PDB Mutation feature file created!')

def makeheader(gene, start, stop, ref_fastaFile):
    '''Using the reference sequence fasta file, we retrive specific headers of SH2 domains. Here to handle the differences in start and stop positions between experiments and references, we use +- 5 amino acid cut off of start/end domain boundary positions and identify the right domain header
    
    Parameters
    ----------
        gene : str
            gene of interest
        start : int
            start position of the sH2 domain
        stop : int
            end position of SH2 domain
        ref_fastaFile : str
            fasta file with unaligned reference sequences
    Returns
    -------
        header : str
            based on the input, thsi function identifies the fasta header from reference sequence fasta file '''
        
    ref_seq = SeqIO.parse(open(ref_fastaFile), 'fasta')
    ref_headers = []
    for fasta in ref_seq:
        name, sequence = fasta.id, str(fasta.seq)
        ref_headers.append(name)
        
    for entry in ref_headers:
        uid, refgene, index, ref_start, ref_end = entry.split('|')
        if gene == refgene:
            if start in range(int(ref_start)-5, int(ref_start)+5) and stop in range(int(ref_end)-5, int(ref_end)+5):
                header = entry
            
    return header
    
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

def ensembleID_dict(fasta_refFile):
    ''' Using an input fasta reference file, this function generates a dictionary with keys as Ensemble ID (associated to Uniprot ID) and values are list of headers. For tandem SH2 domain proteins, we will find more than one header linked to a single Uniprot/Ensemble ID. '''
    df_id = pd.DataFrame()
    ensmbl_id = []
    headerlist = []
    file = SeqIO.parse(open(fasta_refFile), 'fasta')
    for fasta in file:
        name, sequence = fasta.id, str(fasta.seq)
        uniprot_id = name.split('|')[0]
        ensmbl_id.append(get_transcriptID(uniprot_id))
        headerlist.append(name)
    df_id['header'] = headerlist
    df_id['ensmbl_id'] = ensmbl_id

    header_EnsmblID_dict ={}
    for name, group in df_id.groupby('ensmbl_id'):
        header_EnsmblID_dict[name]= group['header'].unique().tolist()

    return header_EnsmblID_dict

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

def get_MIMID(uniprot_id):
    """ Returns the MIM ID (OMIM accession number) for a given UniProt ID """
    get_url = requests.get(f'http://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}')
    if get_url.status_code == 200:
        response = get_url.json()
    mim_id = ''  
    for i in range(len(response['dbReferences'])):
        for j in response['dbReferences']:
            if 'MIM' in j.values():
                if response['dbReferences'][i]['type'] == 'MIM':
                    if response['dbReferences'][i]['properties']['type'] == 'gene':
                        mim_id = int(response['dbReferences'][i]['id'])   
            
    return mim_id
    
def separateNumbersAlphabets(str):
    ''' separates numeric and characters in a string'''
    numbers = []
    alphabets = []
    res = re.split('(\d+)', str)
    
    for i in res:
        if i >= '0' and i <= '9':
            numbers.append(i)
        else:
            alphabets.append(i)
            
    return numbers, alphabets

def get_omim_content(api_key, mim_id):
    ''' downloading content from OMIM '''
    api_url = f'https://api.omim.org/api/entry/allelicVariantList?mimNumber={mim_id}'
    headers = {'apiKey':f'{api_key}'}
    response = requests.get(api_url, headers=headers)
    content = response.__dict__['_content']
    return content

def get_omim_mutations(content):
    ''' extracts gene specific mutation from OMIM database
    
    Parameters
    ----------
        content : bytes
            The response object returned for a MIM ID is parsed to extract mutations associated to the MIM ID 
    Returns
    -------
        Outputs a list of mutations for specific gene'''
    
    dat = minidom.parseString(content)
    mutation_list = []
    mutation_positions =[]
    if (dat.getElementsByTagName('mutations')) == []:
        gene = ''
    for i in range(len(dat.getElementsByTagName('mutations'))):
        tagname = dat.getElementsByTagName('mutations')[i].firstChild.nodeValue
        values = tagname.split(',')
        gene = values[0] 
        mutations = values[1:]
        # print(tagname)
        mut = [] 
        for j in mutations:
            
            j1 = j.replace(" ", "")
            if j1[:3] in AA_dict.keys() and j1[-3:] in AA_dict.keys():
                #print(j)
                numbers, chars = separateNumbersAlphabets(j)
                
                short = []
                for c in chars:
                    if c.replace(" ", "") in AA_dict.keys():
                        #print('true')
                        short.append(AA_dict[c.replace(" ", "")])
                        if j not in mutation_list:
                            mutation_list.append(j)
         
                for n in numbers:
                    if n not in mutation_positions:
                        mutation_positions.append(int(n))
        
    return(gene, mutation_positions)

def listToString(s):
 
    str1 = " "
 
    # return string
    return (str1.join(str(s)))

def omim_mutation_dict(uniprot_ids, api_key):
    ''' Creates a dictionary to identify what mutations are extracted for a UniProt ID from OMIM 
    
    Parameters 
    ----------
        uniprot_ids : list
            list of uniprot IDS to assess
        api_key : str
            this key needs to be generated through OMIM to be able to access their programmatic interface 
    Returns
    -------
        mutation_dict : dict
            dictionary whose keys are uniprot IDS and values are a list of gene name and mutation list in the form of a string
            '''
    mutation_dict = {}
    for id in uniprot_ids:
        print(id)
        mim_id = get_MIMID(id)
        # print(mim_id)
        
        if mim_id != '': 
            api_url = f'https://api.omim.org/api/entry/allelicVariantList?mimNumber={mim_id}' 
            headers = {'apiKey':f'{api_key}'}  
            response = requests.get(api_url, headers=headers) 
            content = response.__dict__['_content'] 
            gene, mutation_list = get_omim_mutations(content) 
            mut_list_str = ','.join(map(str, set(mutation_list))) 
            # print(id, gene, mutation_list, mut_list_str) 
            mutation_dict[id] = [gene, mut_list_str]
            
    return mutation_dict

def get_mutations_in_domains(uniprot_ref_file, mutations_file, Interpro_ID):
    """
    Maps mutations from a mutations file to the domains in a UniProt reference file based on InterPro ID.


    Parameters
    ----------
        uniprot_ref_file : str
            Path to the UniProt reference file.
        mutations_file : str
            Path to the mutations file.
        Interpro_ID : str
            The InterPro ID to filter domains.

    Returns
    -------
        pd.DataFrame
            DataFrame containing mapped mutations to domains.
    """
    uniprot_reference_df = pd.read_csv(uniprot_ref_file)
    df_mut = pd.read_csv(mutations_file)

    # Check that the uniprot and mutation dataframes have the necessary columns
    if 'UniProt ID' not in uniprot_reference_df.columns:
        raise ValueError("UniProt reference DataFrame must contain 'UniProt ID' column.")
    if 'Ref Sequence' not in uniprot_reference_df.columns:
        raise ValueError("UniProt reference DataFrame must contain 'Ref Sequence' column.")
    if 'UniProt ID' not in df_mut.columns:
        raise ValueError("Mutations DataFrame must contain 'UniProt ID' column.")
    if 'mutation' not in df_mut.columns:
        raise ValueError("Mutations DataFrame must contain 'mutation' column.")


    uniprot_mapping_dict, IDs_not_mapped = map_report_uniprot_ids(uniprot_reference_df, df_mut)
    df_mutations_in_domain, error_dict = map_mutations_to_domains(uniprot_mapping_dict, df_mut, uniprot_reference_df, Interpro_ID)

    #Let's write a report of errors
    print(f"Uniprot IDs not mapped between the mutation file and the reference file: {len(IDs_not_mapped)} \n\t")
    print(IDs_not_mapped)

    for uniprot_id, mutations in error_dict.items():
        if uniprot_id not in IDs_not_mapped:
            print(f"Uniprot ID: {uniprot_id} has the following mutation sites that do not match the reference sequence: \n\t {mutations}")

    return df_mutations_in_domain


    
def map_report_uniprot_ids(uniprot_reference_df, mutations_df):
    """
    Maps Uniprot IDs from the mutations DataFrame to the reference DataFrame.
    This function checks if the Uniprot IDs in the mutations DataFrame are present in the reference DataFrame.
    If a Uniprot ID is not found, it checks if it is a canonical isoform (ending with -1) and maps it to the canonical ID.
    If the Uniprot ID is not found in the reference DataFrame, it is added to a list of IDs not mapped.
    Parameters
    ----------
        uniprot_reference_df : pd.DataFrame
            DataFrame containing UniProt reference sequences with a column 'UniProt ID'.
        mutations_df : pd.DataFrame
            DataFrame containing mutations with a column 'Uniprot ID'.  
    Returns 
    -------
        uniprot_mapping_dict : dict
            Dictionary mapping Uniprot IDs to their canonical forms.
        IDs_not_mapped : list
            List of Uniprot IDs that could not be mapped to the reference DataFrame.
    """
    uniprot_IDs = uniprot_reference_df['UniProt ID'].to_list()

    if 'UniProt ID' not in mutations_df.columns:
        print("No Uniprot ID column in mutations, please make sure the mutations DataFrame has a 'Uniprot ID' column.")
        return {}, []

    uniprot_mutations = list(set(mutations_df['UniProt ID'].to_list()))
    uniprot_mapping_dict = {}
    IDs_not_mapped = []
    for ID in uniprot_mutations:
        if ID not in uniprot_IDs:
            #check to see if it's a canonical isoform with a -1 at the end
            if ID.endswith('-1'):
                canonical_ID = ID[:-2]
                if canonical_ID in uniprot_IDs:
                    uniprot_mapping_dict[ID] = canonical_ID
                else:
                    # If the canonical ID is not found, add to IDs_not_mapped
                    IDs_not_mapped.append(ID)
            else:
                IDs_not_mapped.append(ID)
        else:
            uniprot_mapping_dict[ID] = ID
    return uniprot_mapping_dict, IDs_not_mapped

def map_mutations_to_domains(uniprot_mapping_dict, df_mut, uniprot_reference_df, Interpro_ID):
    """
    Maps mutations to domains in a UniProt reference file based on InterPro ID.
    This assumes you have already run uniprot_mapping_dict, IDs_not_mapped = map_report_uniprot_ids(uniprot_reference_file, mutations_file)
    Here, we walk through each mutation, check that it is correctly the same amino acid in the
    uniprot reference sequence, then see if it falls within a domain of interest and print this to a feature file.
    
    Parameters
    ----------
        uniprot_mapping_dict : dict
            Dictionary mapping Uniprot IDs to their canonical forms.
        df_mut : pd.DataFrame
            DataFrame containing mutations with 'Uniprot ID' and 'mutation' columns.
        uniprot_reference_df : pd.DataFrame
            DataFrame containing UniProt reference sequences and domain information.
        Interpro_ID : str
            The InterPro ID to filter domains.      

    Returns
    -------
        df_mutations_in_domain : pd.DataFrame
            DataFrame containing mutations mapped to domains. Has columns:
            - 'Uniprot ID'
            - 'mutation Uniprot ID' - the Uniprot ID of the mutation
            - 'mutation' - original mutation string
            - 'Original Sequence Number' - original sequence number of the mutation
            - 'fasta header' - fasta header of the domain this mutation maps into 
            - 'Domain Sequence Number' - the sequence number of the mutation within the domain
        error_dict : dict
            Dictionary containing errors for mutations that do not match the reference sequence.    

    """
    df_mutations_in_domain = pd.DataFrame(columns=['Uniprot ID', 'mutation Uniprot ID', 'mutation', 'Original Sequence Number', 'fasta header', 'Domain Sequence Number'])
    row_list = []
    error_dict = {}
    for uniprot_id, group in df_mut.groupby('UniProt ID'):
        if uniprot_id not in uniprot_mapping_dict:
            continue
        seq = uniprot_reference_df[uniprot_reference_df['UniProt ID'] == uniprot_mapping_dict[uniprot_id]]['Ref Sequence'].values[0]    
        domains = uniprot_reference_df[uniprot_reference_df['UniProt ID'] == uniprot_mapping_dict[uniprot_id]]['Interpro Domains'].values[0]  # Assuming this column contains domain information
        domain_arr = domains.split(';')
        domain_dict = {}
        num = 1
        aa_list_in_domains = []
        domain_fasta_headers = return_long_fasta_header(uniprot_reference_df.loc[uniprot_reference_df['UniProt ID'] == uniprot_mapping_dict[uniprot_id]], Interpro_ID)

        for domain in domain_arr:
            name, domain_intepro_ID, start, end = domain.split(':')
            if domain_intepro_ID == Interpro_ID:
                domain_dict[num] = {'start':start, 'end':end, 'name':name, 'Interpro_ID':domain_intepro_ID}
                aa_list_in_domains.extend(range(int(start), int(end) + 1))  # +1 to include the end site
                num+=1
        # now I have start and end areas I care about, let's go ahead and create a list of sequeneces that matter to me and if the mutation is in that list, then I will process it
        for index, row in group.iterrows():
            mutation = row['mutation']
            # check that the mutation is in the correct format, that it starts and ends with a character and has a number in the middle
            if not mutation[0].isalpha() or not mutation[-1].isalpha() or not mutation[1:-1].isdigit():
                #print(f"Mutation {mutation} at index {index} is not in the correct format, skipping.")
                continue
            aa = mutation[0]
            to_aa = mutation[-1]
            site_num = int(mutation[1:-1])
            if site_num in aa_list_in_domains:
                # This mutation is within the domain of interest, process it
                #print(f"Mutation {mutation} at site {site_num} is within the domain of interest for Uniprot ID {uniprot_id}.")
                if seq[site_num - 1] != aa:
                    #print(f"Warning: Mutation {mutation} at site {site_num} does not match the reference sequence for Uniprot ID {uniprot_id}. Reference sequence has {seq[site_num - 1]} at this position.")
                    if uniprot_id not in error_dict:
                        error_dict[uniprot_id] = []
                    error_dict[uniprot_id].append(mutation)
                    continue
                if len(domain_dict) == 1:
                    # If there is only one domain, we can set the domain sequence number to based on the start site
                    domain_start = int(domain_dict[1]['start'])
                    domain_sequence_number = site_num - domain_start + 1
                    header = domain_fasta_headers[1]
                else:
                    for domain_num, domain_info in domain_dict.items():
                        if int(domain_info['start']) <= site_num <= int(domain_info['end']):
                            domain_sequence_number = site_num - int(domain_info['start']) + 1
                            header = domain_fasta_headers[domain_num]
                            break
                temp_dict = {
                    'Uniprot ID': uniprot_mapping_dict[uniprot_id],
                    'mutation Uniprot ID': uniprot_id,
                    'mutation': mutation,
                    'Original Sequence Number': site_num,
                    'fasta header': header,
                    'Domain Sequence Number': domain_sequence_number #int(site_num) - int(domain_dict[0]['start']) + 1  # Adjusting to domain sequence number
                }
                row_list.append(temp_dict)
    df_mutations_in_domain = pd.DataFrame(row_list)
    return df_mutations_in_domain, error_dict

def return_long_fasta_header(uniprot_reference_row, Interpro_ID):
    """ 
    Given a row from the uniprot reference file, return a long fasta header for each domain of interest in a
    domain dict, with keys equal to the domain number and values equal to the header.

    Parameters
    ----------
        uniprot_reference_row : pd.DataFrame
            A single row of a DataFrame containing UniProt reference information.
        Interpro_ID : str
            The InterPro ID to filter domains.
    Returns
    -------
        header_dict : dict
            A dictionary where keys are domain numbers and values are formatted fasta headers.
            The format is:
            ">UniProt ID|Gene Name|Species|Domain Name|Domain Number|InterPro ID|Start|End"
    """
    if(len(uniprot_reference_row) != 1):
        print("Expecting a single row of a dataframe, this will return the last row only")
    for index, row in uniprot_reference_row.iterrows():
        header_dict = {}
        domainNum = 1
        domains = row['Interpro Domains']
        domainsArr = domains.split(';')
        uniprot_id = row['UniProt ID']
        gene_name = row['Gene']
        species = row['Species']
        domainNum = 1
        for domain in domainsArr:
            domain_name, interpro, start, end = domain.split(':')
            if interpro == Interpro_ID:
                #>O95696|BRD1|Homo sapiens|Bromodomain|1|IPR001487|560|668
                header_dict[domainNum] = ">%s|%s|%s|%s|%d|%s|%d|%d"%(uniprot_id, gene_name, species, domain_name, domainNum, Interpro_ID, int(start), int(end))
                domainNum += 1
    return header_dict

def print_domain_mutation_features(df_mutations_in_domain, mapping_header_dict, color_dict, output_feature_file):
    """
    Prints a feature file for mutations in domains, where each mutation is represented as a feature in the format required by Jalview.
    This function takes a DataFrame of mutations in domains and generates a feature file that can be used in Jalview for visualization.
    Parameters
    ----------
        df_mutations_in_domain : pd.DataFrame
            DataFrame containing mutations mapped to domains such as returned by `get_mutations_in_domains`.
        color_dict : dict
            Dictionary mapping mutation types to colors for visualization in Jalview.
        output_feature_file : str
            Path to the output feature file where the mutations will be written.    

    """
    mutation_feature_dict = {}
    for name, group in df_mutations_in_domain.groupby('fasta header'):
        # get the unique mutations for this fasta header
        #make a mapping to the short header
        name = name.split('>')[1]
        if name in mapping_header_dict:
            name_short = mapping_header_dict[name]['short']
        else:
            print("ERROR: Fasta header not found in mapping dictionary, using full header.")
            name_short = name
        unique_mutations = group['Domain Sequence Number'].unique()
        mutation_feature_dict[name_short] ={}
        for residue in unique_mutations:
            mutation_feature_dict[name_short][residue] = ['mutation']
    featureTools.print_jalview_feature_file(mutation_feature_dict, color_dict, output_feature_file)
    print("Feature file written to %s"%output_feature_file)  