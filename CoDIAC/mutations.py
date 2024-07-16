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
        N_offset : int
            number of residues to add or remove from the N terminal end of the domain boundary
        C_offset : int
            number of residues to add or remove from the C terminal end of the domain boundary

    Returns
    -------
        output_feafile : str
            Mutation feature file (.fea) '''
    
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
    """Generates a feature file with mutations extracted from OMIM database.

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
    Returns
    -------
        output_featureFile : str """
        
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
    '''Generates a feature file for mutations extracted from PDB structures on the domain of interest.
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
        feature file with mutations/variants reported in PDB structures that are present within the SH2 domain boundary '''
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

    
            