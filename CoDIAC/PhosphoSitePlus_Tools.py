import pandas as pd
import os 
from Bio import SeqIO
import codecs


#Given licensing terms of PhosphoSitePlus, you must create a download and point to your 
#local directory of Phosphosite data and point that here. 
package_directory = os.path.dirname(os.path.abspath(__file__))
PHOSPHOSITE_FILE = package_directory + '/data/phosphositeplus_data.csv'


def convert_pSiteDataFiles(PHOSPHOSITEPLUS_DATA_DIR):
    """
    First, download all the data from PhosphositePlus and place it in where you will 
    reference as PHOSPHOSTIEPLUS_DATA_DIR. 
    
    Given the files as they are downloaded from PhosphositePlusrearrange this rearranges 
    data from those files to create a dataframe, written to the CoDIAC data location
    as a dataframe that can tehn be used to query sequence and PTMs by a Uniprot ID. 

    This code will need to be updated in the PhosphositePlus website changes their data format.
    Currently assumes the files are the following and each of a 3 line header preamble:
    Phosphorylation_site_dataset
    Ubiquitination_site_dataset
    Sumoylation_site_dataset
    O-GalNAc_site_dataset
    Phosphosite_PTM_seq.fasta

    Returns
    -------
    df : pandas dataframe
        A dataframe with the UniprotID, species, protein name, sequence, and PTMs for each protein

    """

    phospho = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Phosphorylation_site_dataset')
    ubiq = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Ubiquitination_site_dataset')
    sumo = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Sumoylation_site_dataset')
    glyco = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'O-GalNAc_site_dataset')
    sequence_file = PHOSPHOSITEPLUS_DATA_DIR+'Phosphosite_PTM_seq.fasta'

    #check that all the files exist
    if not os.path.exists(sequence_file):
        print("ERROR: %s does not exist"%(sequence_file))
        return None
    if not os.path.exists(phospho):
        print("ERROR: %s does not exist"%(phospho))
        return None
    if not os.path.exists(ubiq):
        print("ERROR: %s does not exist"%(ubiq))
        return None
    if not os.path.exists(sumo):
        print("ERROR: %s does not exist"%(sumo))
        return None
    if not os.path.exists(glyco):
        print("ERROR: %s does not exist"%(glyco))
        return None
    

    #first remove the non-fasta lines at the preamble of the file, then load the rest of the data
    with codecs.open(sequence_file, 'r', encoding='utf-8',
                 errors='ignore') as f:
        lines = f.readlines()
    for line_number in range(0, len(lines)):
        if lines[line_number][0] != '>': #means we haven't yet hit the first fasta line
            lines.pop(line_number)
        else:
            break 
    file_temp = PHOSPHOSITEPLUS_DATA_DIR+'Phosphosite_PTM_seq_temp.fasta'
    with open(file_temp, 'w') as f:
        f.writelines(lines)

    DEBUG = 0
    count = 0
    df = pd.DataFrame(columns=['Uniprot_ID', 'species', 'name', 'sequence'])
    with open(file_temp, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if DEBUG and count > 100:
                #return df
                break
            temp_dict = {}
            #record_id = record.id
            record_id = record.description #changed to description. record.id ends at first space, which happens in species
            seq = record.seq
            gn, temp_dict['name'], temp_dict['species'], temp_dict['Uniprot_ID'] = record_id.split('|')
            df.loc[len(df)] = [temp_dict['Uniprot_ID'], temp_dict['species'], temp_dict['name'], str(seq)]
            #record_id_vals = record_id.split('|(?=[^ 
            count+=1 
    df.set_index('Uniprot_ID', inplace=True) #this is the parent sequence and info df. 

    #next, let's handle the PTM files, put them as a dataframe then, merge their PTMs into a string 
    #for appending to the larger parent df.
    #merge the PTMs into a string
    for uniprot_id, row in df.iterrows():
        print("DEBUG: working on %s"%(uniprot_id))
        PTM_list = []
        if uniprot_id in phospho:
            PTM_list += phospho[uniprot_id]
        if uniprot_id in ubiq:
            PTM_list += ubiq[uniprot_id]
        if uniprot_id in sumo:
            PTM_list += sumo[uniprot_id]
        if uniprot_id in glyco:
            PTM_list += glyco[uniprot_id]
        PTM_string = ';'.join(PTM_list)
        df.loc[uniprot_id, 'modifications'] = PTM_string
        print("DEBUG, adding PTMs %s"%(PTM_string))

    print("Writing New Phosphosite Data to %s"%(PHOSPHOSITE_FILE))
    df.to_csv(PHOSPHOSITE_FILE)

    return df

def read_PTM_file_to_df(file):
    """
    Have to skip the first 3 lines of the PTM files, then read the rest of the data. 
    
    """
    df = pd.read_csv(file, sep='\t', skiprows=3)
    PTM_dict = {}
    DEBUG = 0
    counter = 0
    for index, row in df.iterrows():
        if DEBUG and counter > 100:
            return PTM_dict
        uniprot_id = row['ACC_ID']
        mod_rsd = row['MOD_RSD']
        pos, type = mod_rsd.split('-')
        if type == 'p':
            if 'S' in pos:
                type_name = 'Phosphoserine'
            elif 'T' in pos:
                type_name = 'Phosphothreonine'
            elif 'Y' in pos:
                type_name = 'Phosphotyrosine'
        elif type == 'ub':
            type_name = 'Ubiquitination'
        elif type == 'sm':
            type_name = 'Sumoylation'
        elif type == 'ga':
            type_name = 'N-Glycosylation'
        else: 
            print("ERROR: don't recognize PTM type %s"%(type))
        if uniprot_id not in PTM_dict:
            PTM_dict[uniprot_id] = []
        PTM_dict[uniprot_id].append(pos+'-'+type_name)
        counter +=1 
    return PTM_dict


