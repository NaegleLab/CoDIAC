import pandas as pd
import os 
from Bio import SeqIO
import codecs


#Given licensing terms of PhosphoSitePlus, you must create a download and point to your 
#local directory of Phosphosite data and point that here. 
package_directory = os.path.dirname(os.path.abspath(__file__))
PHOSPHOSITE_FILE = package_directory + '/data/phosphositeplus_data.csv'
if os.path.exists(PHOSPHOSITE_FILE):
    PHOSPHOSITE = pd.read_csv(PHOSPHOSITE_FILE, index_col=0)
    PSITE_INIT = True
else:
    PSITE_INIT = False



def get_PTMs(uniprot_ID):
    """
    Get PTMs for a given Uniprot ID from the PhosphoSitePlus data. You must have created
    the data file from the PhosphoSitePlus using convert_pSiteDataFiles.

    Parameters
    ----------
    uniprot_ID : str
        Uniprot ID for the protein of interest
    
    Returns
    -------
    PTMs : tuples
        Returns a list of tuples of modifications
        [(position, residue, modification-type),...,]
        
        Returns -1 if unable to find the ID

        Returns [] (empty list) if no modifications  
    Uniprot_ID: string
        The Uniprot ID of the protein of interest, this may change if mods are not found on the root protein
    
    """
    USE_ISO = False
    if PSITE_INIT:
        if uniprot_ID in PHOSPHOSITE.index:
            mods = PHOSPHOSITE.loc[uniprot_ID, 'modifications']
            if mods == 'nan' or pd.isna(mods):
                NO_MODS = True
                print("Found no mods for %s"%(uniprot_ID))
                # There are some cases where the root has no data, but an isoform does, so let's test that
                iso_num = 1
                uniprot_ID_iso = uniprot_ID + '-'+str(iso_num)
                while uniprot_ID_iso in PHOSPHOSITE.index and NO_MODS:
                    mods = PHOSPHOSITE.loc[uniprot_ID_iso, 'modifications']
                    if mods == 'nan' or pd.isna(mods):
                        iso_num += 1
                        uniprot_ID_iso = uniprot_ID + '-'+str(iso_num)
                    else:
                        USE_ISO = True
                        print("Using an isoform for PTMs %s, found mods"%(uniprot_ID_iso))
                        print(mods)
                        NO_MODS = False
                        break
                
                #out of the while loop, see if we have to return empty
                if NO_MODS:
                    return [], uniprot_ID #hit the end of the isoforms and found no mods in while, so leaving with empty

            mods_raw=mods.split(";")
            mods_clean =[]
            for i in mods_raw:
                tmp = i.strip()
                tmp = tmp.split("-")
            
            # append a tuple of (position, residue, type)
                mods_clean.append((tmp[0][1:], tmp[0][0], "-".join(tmp[1:])))
            return_ID = uniprot_ID
            if USE_ISO:
                return_ID = uniprot_ID_iso
            return mods_clean, return_ID
        else:
            print("ERROR: %s not found in PhosphositePlus data"%(uniprot_ID))
            return '-1', uniprot_ID
    else:
        print("ERROR: PhosphositePlus data not found. Run convert_pSiteDataFiles to create it")
        return None, uniprot_ID

def get_sequence(uniprot_ID):
    """
    Get the sequence for a given Uniprot ID from the PhosphoSitePlus data. You must have created
    the data file from the PhosphoSitePlus using convert_pSiteDataFiles.

    Parameters
    ----------
    uniprot_ID : str
        Uniprot ID for the protein of interest
    
    Returns
    -------
    sequence : str
        The sequence of the protein of interest. 
        Returns '-1' if sequence not found
    
    """
    if PSITE_INIT:
        if uniprot_ID in PHOSPHOSITE.index:
            return PHOSPHOSITE.loc[uniprot_ID, 'sequence']
        else:
            #print("ERROR: %s not found in PhosphositePlus data"%(uniprot_ID))
            return '-1'
    else:
        print("ERROR: PhosphositePlus data not found. Run convert_pSiteDataFiles to create it")
        return None

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

    #check that all the files exist
    sequence_file = PHOSPHOSITEPLUS_DATA_DIR+'Phosphosite_seq.fasta'

    if not os.path.exists(sequence_file):
        print("ERROR: %s does not exist"%(sequence_file))
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'Phosphorylation_site_dataset'):
        print("ERROR: Phosphorylation file Phosphorylation_site_dataset does not exist")
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'Ubiquitination_site_dataset'):
        print("ERROR: Ubiquitination file Ubiquitination_site_dataset does not exist")
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'Sumoylation_site_dataset'):
        print("ERROR: Sumoylation file Sumoylation_site_dataset does not exist")
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'O-GalNAc_site_dataset'):
        print("ERROR: O-GalNAc file O-GalNAc_site_dataset does not exist")
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'O-GlcNAc_site_dataset'):
        print("ERROR: O-GlcNAc file O-GlcNAc_site_dataset does not exist")
        return None
    if not os.path.exists(PHOSPHOSITEPLUS_DATA_DIR+'Methylation_site_dataset'):
        print("ERROR: Methylation file Methylation_site_dataset does not exist")
        return None

    phospho = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Phosphorylation_site_dataset')
    ubiq = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Ubiquitination_site_dataset')
    sumo = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Sumoylation_site_dataset')
    galnac = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'O-GalNAc_site_dataset')
    glcnac = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'O-GlcNAc_site_dataset')
    meth = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Methylation_site_dataset')
    acetyl = read_PTM_file_to_df(PHOSPHOSITEPLUS_DATA_DIR+'Acetylation_site_dataset')



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
        #print("DEBUG: working on %s"%(uniprot_id))
        PTM_list = []
        if uniprot_id in phospho:
            PTM_list += phospho[uniprot_id]
        if uniprot_id in ubiq:
            PTM_list += ubiq[uniprot_id]
        if uniprot_id in sumo:
            PTM_list += sumo[uniprot_id]
        if uniprot_id in galnac:
            PTM_list += galnac[uniprot_id]
        if uniprot_id in glcnac:
            PTM_list += glcnac[uniprot_id]
        if uniprot_id in meth:
            PTM_list += meth[uniprot_id]
        if uniprot_id in acetyl:
            PTM_list += acetyl[uniprot_id]
        
        PTM_string = ';'.join(PTM_list)
        df.loc[uniprot_id, 'modifications'] = PTM_string
        #print("DEBUG, adding PTMs %s"%(PTM_string))

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
            if 'N' in pos:
                type_name = 'N-Glycosylation'
            elif 'S' in pos: #if you want a more general classification, this and T can be either just O-GalNAc or O-glycosylation
                type_name = 'O-GalNAc Serine'
            elif 'T' in pos:
                type_name = 'O-GalNAc Threonine'
        elif type == 'gl':
            if 'S' in pos: #if you want a more general classification, this can be either just O-GlcNAc or O-glycosylation
                type_name = 'O-GlcNAc Serine'
            elif 'T' in pos:
                type_name = 'O-GlcNAc Threonine'
        elif type == 'm1' or type == 'me':
            type_name = 'Methylation'
        elif type == 'm2':
            type_name = 'Dimethylation'
        elif type == 'm3':
            type_name = 'Trimethylation'
        elif type == 'ac':
            type_name = 'Acetylation'
        else: 
            print("ERROR: don't recognize PTM type %s"%(type))

        if uniprot_id not in PTM_dict:
            PTM_dict[uniprot_id] = []
        PTM_dict[uniprot_id].append(pos+'-'+type_name)
        counter +=1 
    return PTM_dict


